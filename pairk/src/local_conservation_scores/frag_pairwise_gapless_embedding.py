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
import torch
from attrs import asdict, define, field, validators


@define
class PairwiseKmerEmbAlnParams:
    input_alignment_file: str | Path
    reference_id: str = field(validator=validators.instance_of(str))
    k: int = field(converter=int)
    idr_aln_st: int = field(converter=int)
    idr_aln_end: int = field(converter=int)
    mod: esm_tools.ESM_Model



def make_empty_kmer_ortho_df(positions, ortholog_ids: list[str]):
    cols = ["reference_kmer"] + ortholog_ids
    df = pd.DataFrame(
        index=positions,
        columns=cols,
    )
    return df


def get_kmers(seq: str, k: int):
    """Generates kmers from a sequence"""
    k2 = k - 1
    kmers = []
    for i in range(len(seq) - k2):
        kmers.append(seq[i : i + k])
    return kmers


def get_subsequences(indices, input_string, k):
    """Get subsequences of ortho string based on best indices"""
    char_list = list(input_string)
    start_indices = indices.view(-1, 1)
    end_indices = start_indices + k
    max_length = len(char_list)
    subsequences = [
        "".join(char_list[start:end]) for start, end in zip(start_indices, end_indices)
    ]
    return subsequences


def run_pairwise_kmer_emb_aln(
    reference_id: str,
    embedding_dict: dict,
    k: int,
):
    # get the reference sequence and remove it from the embedding_dict
    ref_seq_str, ref_seq_embedding = embedding_dict.pop(reference_id)
    kmers = get_kmers(ref_seq_str, k)
    positions = list(range(len(ref_seq_str) - (k - 1)))
    score_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    subseq_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    pos_df = make_empty_kmer_ortho_df(positions, list(embedding_dict.keys()))
    score_df.loc[positions, "reference_kmer"] = kmers
    subseq_df.loc[positions, "reference_kmer"] = kmers
    pos_df.loc[positions, "reference_kmer"] = kmers

    # Make expanded ref tensor
    expand_inds_ref = torch.arange(k).view(1, -1) + torch.arange(
        ref_seq_embedding.shape[0]
    ).view(-1, 1)
    expand_inds_ref[ref_seq_embedding.shape[0] - (k - 1) :] = 0
    expand_inds_ref = (
        expand_inds_ref.unsqueeze(-1)
        .expand(-1, -1, ref_seq_embedding.shape[1])
        .to(dtype=torch.int64)
    )
    expand_ref = ref_seq_embedding.unsqueeze(1).expand(-1, expand_inds_ref.shape[1], -1)
    expand_ref = torch.gather(expand_ref, 0, expand_inds_ref)
    expand_ref = expand_ref[: ref_seq_embedding.shape[0] - (k - 1)].reshape(
        -1, k * expand_ref.shape[2]
    )

    # for each ortholog sequence
    for ortholog_id, v in embedding_dict.items():
        ortholog_seq = v[0]
        ortholog_embedding = v[1]
        if ortholog_seq is None or ortholog_embedding is None:
            continue
        if len(ortholog_seq) < k:
            continue

        # Make expanded ortholog tensor
        expand_inds_ortho = torch.arange(k).view(1, -1) + torch.arange(
            ortholog_embedding.shape[0]
        ).view(-1, 1)
        expand_inds_ortho[ortholog_embedding.shape[0] - (k - 1) :] = 0
        expand_inds_ortho = (
            expand_inds_ortho.unsqueeze(-1)
            .expand(-1, -1, ortholog_embedding.shape[1])
            .to(dtype=torch.int64)
        )
        expand_ortho = ortholog_embedding.unsqueeze(1).expand(
            -1, expand_inds_ortho.shape[1], -1
        )
        expand_ortho = torch.gather(expand_ortho, 0, expand_inds_ortho)
        expand_ortho = expand_ortho[: ortholog_embedding.shape[0] - (k - 1)].reshape(
            -1, k * expand_ortho.shape[2]
        )

        # Calculate pairwise distances and get stats
        pairwise_dists = torch.cdist(
            expand_ref, expand_ortho, p=2
        )  # Optional: compute_mode='donot_use_mm_for_euclid_dist'
        min_dists, min_dists_pos = torch.min(pairwise_dists, dim=-1)
        score_df.loc[positions, ortholog_id] = min_dists.cpu().numpy()
        pos_df.loc[positions, ortholog_id] = min_dists_pos.cpu().numpy()
        subseq_df.loc[positions, ortholog_id] = get_subsequences(
            min_dists_pos, ortholog_seq, k
        )
    return score_df, subseq_df, pos_df


def get_idr_embeddings(aln_seq: str, idr_aln_start: int, idr_aln_end: int, mod: esm_tools.ESM_Model, device='cuda', threads=1):
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
    orth_tensor = mod.encode(seq_str, device=device, threads=threads)
    # orth_tensor = mod.encode(seq_str, device='cpu', threads=60)
    idr_ortho_tensor = orth_tensor[idr_orth_st + 1 : idr_orth_end + 2, :]
    return idr_str, idr_ortho_tensor


def get_idr_embedding_dict(
    aln: Align.MultipleSeqAlignment, idr_aln_st: int, idr_aln_end: int, mod: esm_tools.ESM_Model, **kwargs
):
    embedding_dict = defaultdict(list)
    for seq in aln:
        idr_str, idr_tensor = get_idr_embeddings(
            str(seq.seq), idr_aln_st, idr_aln_end, mod, **kwargs
        )
        embedding_dict[seq.id].append(idr_str)
        embedding_dict[seq.id].append(idr_tensor)
    return embedding_dict


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
    embedding_dict = get_idr_embedding_dict(aln, idr_aln_st, idr_aln_end, mod, **kwargs)
    score_df, subseq_df, pos_df = run_pairwise_kmer_emb_aln(
        reference_id,
        embedding_dict,
        k,
    )
    output_dict = {
        "score_dataframe": score_df.to_dict(orient="split"),
        "subseq_dataframe": subseq_df.to_dict(orient="split"),
        "position_dataframe": pos_df.to_dict(orient="split"),
        # "reciprocal_best_match_dataframe": rbm_df.to_dict(orient="split"),
    }
    return output_dict


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
    output_dict = main_no_output_file(
        input_alignment_file,
        reference_id,
        k,
        idr_aln_st,
        idr_aln_end,
        mod,
        **kwargs,
    )
    with open(output_file, "w") as json_file:
        json.dump(output_dict, json_file)



    




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
