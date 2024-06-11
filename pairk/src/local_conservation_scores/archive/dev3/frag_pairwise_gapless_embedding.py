# %%
import json
import time
from collections import defaultdict
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
from Bio.SeqRecord import SeqRecord

import local_seqtools.esm_tools as esm_tools
from local_conservation_scores.tools import capra_singh_2007_scores, general
from local_env_variables import matrices
from local_seqtools import alignment_tools as aln_tools
from local_seqtools import general_utils as tools
from local_seqtools import jch_alignment as jch_aln


def euclidean_distance(tensor1, tensor2):
    squared_diff = np.square(tensor1 - tensor2)
    sum_squared_diff = np.sum(squared_diff)
    distance = np.sqrt(sum_squared_diff)
    return distance


def make_empty_kmer_ortho_df(positions, ortholog_ids: list[str]):
    cols = ["reference_kmer"] + ortholog_ids
    df = pd.DataFrame(
        index=positions,
        columns=cols,
    )
    return df


def find_best_ortho_match_embeddings(hit_tensor, ortho_str, ortho_tensor):
    pos, dist, seq = [], [], []
    matches = defaultdict(list)
    for i in range(len(ortho_str) - (len(hit_tensor) - 1)):
        target_tensor = ortho_tensor[i : i + len(hit_tensor), :]
        dist_i = euclidean_distance(np.array(hit_tensor), np.array(target_tensor))
        matches[i].append(dist_i)
        matches[i].append(ortho_str[i : i + len(hit_tensor)])
        pos.append(i)
        dist.append(dist_i)
        seq.append(ortho_str[i : i + len(hit_tensor)])
    best_match_ind = np.argmin(dist)
    best_subseq = seq[best_match_ind]
    best_pos = pos[best_match_ind]
    return np.min(dist), best_subseq, best_pos


def generate_kmer_embeddings(seq: str, tensor, k: int):
    """Generates kmers and their corresponding embeddings from a sequence and
    its embedding tensor
    Assumes that the sequence and tensor are of the same length (i.e. no start
    and end tokens are included in the tensor)"""
    k2 = k - 1
    kmers = []
    kmer_tensors = []
    positions = []
    for i in range(len(seq) - k2):
        kmers.append(seq[i : i + k])
        kmer_tensors.append(tensor[i : i + k, :])
        positions.append(i)
    return positions, kmers, kmer_tensors


def run_pairwise_kmer_emb_aln(
    reference_id: str,
    embedding_dict: dict,
    k: int,
):
    ref_idr_str, ref_idr_embedding = embedding_dict.pop(reference_id)
    positions, kmers, kmer_tensors = generate_kmer_embeddings(
        seq=ref_idr_str,
        tensor=ref_idr_embedding,
        k=k,
    )
    score_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    subseq_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    pos_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    rbm_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    for position, kmer, kmer_tensor in zip(positions, kmers, kmer_tensors):
        score_df.loc[position, "reference_kmer"] = kmer
        subseq_df.loc[position, "reference_kmer"] = kmer
        pos_df.loc[position, "reference_kmer"] = kmer
        rbm_df.loc[position, "reference_kmer"] = kmer
        for ortholog_id, v in embedding_dict.items():
            ortholog_idr = v[0]
            ortholog_embedding = v[1]
            if ortholog_idr is None or ortholog_embedding is None:
                continue
            if len(ortholog_idr) < k:
                continue
            try:
                best_score, best_subseq, best_pos = find_best_ortho_match_embeddings(
                    kmer_tensor, ortholog_idr, ortholog_embedding
                )
            except ValueError as e:
                print(e)
                continue
            score_df.loc[position, ortholog_id] = best_score
            subseq_df.loc[position, ortholog_id] = best_subseq
            pos_df.loc[position, ortholog_id] = best_pos
            try:
                reci_best_score, reci_best_subseq, reci_best_pos = find_best_ortho_match_embeddings(
                    ortholog_embedding[best_pos : best_pos + k, :],
                    ref_idr_str,
                    ref_idr_embedding,
                )
            except ValueError as e:
                print(e)
                continue
            # for if the reciprocal best match is the same sequence as the original
            rbm_df.loc[position, ortholog_id] = reci_best_subseq == kmer
            # for if the reciprocal best match is the exact same as the original
            # rbm_df.loc[position, ortholog_id] = reci_best_pos == position
    return score_df, subseq_df, pos_df, rbm_df


def get_idr_embeddings(aln_seq: str, idr_aln_start: int, idr_aln_end: int, mod):
    seq_str, index = tools.reindex_alignment_str(aln_seq)
    idr_aln_str = aln_seq[idr_aln_start : idr_aln_end + 1]
    idr_str, _ = tools.reindex_alignment_str(idr_aln_str)
    if len(idr_str) == 0:
        print("no idr")
        return None, None
    idr_orth_st = seq_str.find(idr_str)
    idr_orth_end = idr_orth_st + len(idr_str) - 1
    if idr_orth_st == -1 or idr_orth_end == -1:
        print('couldn"t find idr in fl ortholog seq')
        return None, None
    orth_tensor = mod.encode(seq_str, device="cuda")
    # orth_tensor = mod.encode(seq_str, device='cpu', threads=60)
    idr_ortho_tensor = orth_tensor[idr_orth_st + 1 : idr_orth_end + 2, :]
    return idr_str, idr_ortho_tensor


def get_idr_embedding_dict(
    aln: Align.MultipleSeqAlignment, idr_aln_st: int, idr_aln_end: int, mod: esm_tools.ESM_Model
):
    embedding_dict = defaultdict(list)
    for seq in aln:
        idr_str, idr_tensor = get_idr_embeddings(
            str(seq.seq), idr_aln_st, idr_aln_end, mod
        )
        embedding_dict[seq.id].append(idr_str)
        embedding_dict[seq.id].append(idr_tensor)
    return embedding_dict


def main(
    input_alignment_file: str | Path,
    output_file: str | Path,
    reference_id: str,
    k: int,
    idr_aln_st: int,
    idr_aln_end: int,
    mod: esm_tools.ESM_Model,
    overwrite: bool = False,
    **kwargs,
):
    # check if output file exists
    if Path(output_file).exists() and not overwrite:
        # print(f'{output_file} exists and overwrite is False')
        # print('exiting...')
        return
    fasta_importer = tools.FastaImporter(input_alignment_file)
    aln = fasta_importer.import_as_alignment()
    embedding_dict = get_idr_embedding_dict(aln, idr_aln_st, idr_aln_end, mod)
    score_df, subseq_df, pos_df, rbm_df = run_pairwise_kmer_emb_aln(
        reference_id,
        embedding_dict,
        k,
    )
    output_dict = {
        "score_dataframe": score_df.to_dict(orient="split"),
        "subseq_dataframe": subseq_df.to_dict(orient="split"),
        "position_dataframe": pos_df.to_dict(orient="split"),
        "reciprocal_best_match_dataframe": rbm_df.to_dict(orient="split"),
    }
    with open(output_file, "w") as json_file:
        json.dump(output_dict, json_file)


def main_no_output_file(
    input_alignment_file: str | Path,
    reference_id: str,
    k: int,
    idr_aln_st: int,
    idr_aln_end: int,
    mod: esm_tools.ESM_Model,
    **kwargs,
):
    fasta_importer = tools.FastaImporter(input_alignment_file)
    aln = fasta_importer.import_as_alignment()
    embedding_dict = get_idr_embedding_dict(aln, idr_aln_st, idr_aln_end, mod)
    score_df, subseq_df, pos_df, rbm_df = run_pairwise_kmer_emb_aln(
        reference_id,
        embedding_dict,
        k,
    )
    output_dict = {
        "score_dataframe": score_df.to_dict(orient="split"),
        "subseq_dataframe": subseq_df.to_dict(orient="split"),
        "position_dataframe": pos_df.to_dict(orient="split"),
        "reciprocal_best_match_dataframe": rbm_df.to_dict(orient="split"),
    }
    return output_dict
    




# %%
# ==============================================================================
# // testing
# ==============================================================================
# import local_conservation_analysis_pipeline.group_conservation_objects as group_tools

# json_file = "../../benchmark/benchmark_multi_species/p3_conservation_analysis_pipeline/conservation_analysis/2-9606_0_004caa/2-9606_0_004caa.json"
# og = group_tools.ConserGene(json_file)
# lvlo=og.get_level_obj('Vertebrata')

# def pad_hit(seq: str, st_pos: int, end_pos: int, l_flank: int = 0, r_flank: int = 0):
#     st = max(0, st_pos - l_flank)
#     end = min(len(seq)-1, end_pos + r_flank)
#     return st, end, seq[st : end + 1]

# _,_,flanked_hit = pad_hit(og.query_sequence, og.hit_start_position, og.hit_end_position, 5, 5)
# k = len(flanked_hit)
# og.hit_sequence
# main(
#     input_alignment_file=lvlo.alignment_file,
#     reference_id=og.query_gene_id,
#     idr_aln_st=lvlo.idr_aln_start,
#     idr_aln_end=lvlo.idr_aln_end,
#     k=k,
#     output_file='test.json',
#     model_name='esm2_t33_650M_UR50D',
#     overwrite=True,
# )
