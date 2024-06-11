# %%
import json
import time
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
from Bio.SeqRecord import SeqRecord

from local_conservation_scores.tools import capra_singh_2007_scores, general
from local_env_variables import matrices
from local_seqtools import alignment_tools as aln_tools
from local_seqtools import general_utils as tools
from local_seqtools import jch_alignment as jch_aln


# %%
def score_alignment(seq1: str, seq2: str, subs_mat_df: pd.DataFrame) -> float:
    assert len(seq1) == len(seq2)
    score = 0
    for s1, s2 in zip(seq1, seq2):
        score += float(subs_mat_df.loc[s1, s2]) # type: ignore
    return score


def score_kmer_2_seq_no_gaps_old(
    kmer: str,
    sequence: str,
    substitution_matrix_df: pd.DataFrame,
):
    '''
    best_score = np.max(scores)
    best_position = np.argmax(scores)
    best_subseq = sequence[best_position : best_position + len(kmer)]
    '''
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")
    scores = []
    for i in range(len(sequence) - len(kmer) + 1):
        subseq = sequence[i : i + len(kmer)]
        # score = aln_tools.score_alignment(kmer, subseq, substitution_matrix_df)
        score = score_alignment(kmer, subseq, substitution_matrix_df)
        scores.append(score)
    return scores


def kmer_best_matches_2_ortho(
    kmer: str,
    sequence: str,
    substitution_matrix_df: pd.DataFrame,
):
    scores = score_kmer_2_seq_no_gaps_old(kmer, sequence, substitution_matrix_df)
    best_score = np.max(scores)
    best_position = np.argmax(scores)
    best_subseq = sequence[best_position : best_position + len(kmer)]
    return best_score, best_subseq, best_position

                    
def make_empty_kmer_ortho_df(kmers: list[str], ortholog_ids: list[str]):
    return pd.DataFrame(
        index=kmers,
        columns=ortholog_ids,
    )

            
def run_pairwise_kmer_alignment(
    ref_idr: str,
    ortholog_idrs: dict[str, str],
    k: int,
    scoring_matrix_name: str="grantham_similarity_norm",
):
    substitution_matrix_df = matrices.load_precomputed_matrix_df(scoring_matrix_name)
    kmers = tools.gen_kmers(ref_idr, k)
    kmers = list(set(kmers))
    score_df = make_empty_kmer_ortho_df(kmers, list(ortholog_idrs.keys()))
    subseq_df = make_empty_kmer_ortho_df(kmers, list(ortholog_idrs.keys()))
    rbm_df = make_empty_kmer_ortho_df(kmers, list(ortholog_idrs.keys()))
    for kmer in kmers:
        for ortholog_id, ortholog_idr in ortholog_idrs.items():
            if len(ortholog_idr) < k:
                continue
            try:
                best_score, best_subseq, _ = kmer_best_matches_2_ortho(
                    kmer,
                    ortholog_idr,
                    substitution_matrix_df,
                )
            except ValueError as e:
                print(e)
                continue
            score_df.loc[kmer, ortholog_id] = best_score
            subseq_df.loc[kmer, ortholog_id] = best_subseq
            try:
                reci_best_score, reci_best_subseq, _ = kmer_best_matches_2_ortho(
                    best_subseq,
                    ref_idr,
                    substitution_matrix_df,
                )
            except ValueError as e:
                print(e)
                continue
            rbm_df.loc[kmer, ortholog_id] = (reci_best_subseq==kmer)
    return score_df, subseq_df, rbm_df
        

def main(
    input_alignment_file: str|Path,
    output_file: str|Path,
    reference_id: str,
    k: int,
    idr_aln_st: int,
    idr_aln_end: int,
    overwrite: bool = False,
    scoring_matrix_name: str="grantham_similarity_norm",
    **kwargs,
):
    # check if output file exists
    if Path(output_file).exists() and not overwrite:
        # print(f'{output_file} exists and overwrite is False')
        # print('exiting...')
        return
    fasta_importer = tools.FastaImporter(input_alignment_file)
    aln = fasta_importer.import_as_alignment()
    idr_aln = aln[:, idr_aln_st : idr_aln_end+1]
    idrs = tools.strip_dashes_from_sequences(list(idr_aln)) # type: ignore
    idrs = {i.id: str(i.seq) for i in idrs}
    ref_idr = idrs.pop(reference_id)
    score_df, subseq_df, rbm_df = run_pairwise_kmer_alignment(
        ref_idr,
        idrs, # type: ignore
        k,
        scoring_matrix_name,
    )
    output_dict = {
        "score_dataframe": score_df.to_dict(orient="split"),
        "subseq_dataframe": subseq_df.to_dict(orient="split"),
        "reciprocal_best_match_dataframe": rbm_df.to_dict(orient="split"),
    }
    with open(output_file, "w") as json_file:
        json.dump(output_dict, json_file, indent=4)




# %%
# ==============================================================================
# // testing
# ==============================================================================    
# import local_conservation_analysis_pipeline.group_conservation_objects as group_tools

# json_file = "../../../../05-conservation_pipeline/examples/table_annotation/conservation_analysis/3-9606_0_002f40/3-9606_0_002f40.json"
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
#     output_file='test.json'
# )



# %%
# ==============================================================================
# // old
# ==============================================================================    





def score_kmer_2_seq_no_gaps_old2(
    kmer: str,
    sequence: str,
    substitution_matrix_df: pd.DataFrame,
):
    '''
    best_score = np.max(scores)
    best_position = np.argmax(scores)
    best_subseq = sequence[best_position : best_position + len(kmer)]
    '''
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")
    bscore = 0
    for i in range(len(sequence) - len(kmer) + 1):
        subseq = sequence[i : i + len(kmer)]
        # score = aln_tools.score_alignment(kmer, subseq, substitution_matrix_df)
        score = score_alignment(kmer, subseq, substitution_matrix_df)
        if score > bscore:
            bscore = score
            bsubseq = subseq
            bpos = i
    return bscore, bsubseq, bpos


def score_2_kmers(kmer1: str, kmer2: str, substitution_matrix_df: pd.DataFrame):
    kmer1_arr = np.array(list(kmer1))
    kmer2_arr = np.array(list(kmer2))
    return np.trace(substitution_matrix_df.loc[kmer1_arr, kmer2_arr].values)


def score_2_kmers_arrays(kmer1_arr: np.ndarray, kmer2_arr: np.ndarray, substitution_matrix_df: pd.DataFrame):
    return np.trace(substitution_matrix_df.loc[kmer1_arr, kmer2_arr].values)

def score_kmer_2_seq_no_gaps_v1(
    kmer: str,
    sequence: str,
    substitution_matrix_df: pd.DataFrame,
):
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")

    kmer_len = len(kmer)
    sequence_len = len(sequence)
    kmerarr = np.array(list(kmer))
    seqarr = np.array(list(sequence))
    # maxes = [max]
    score_list = []
    # Calculate scores for all possible subsequences
    for i in range(sequence_len - kmer_len + 1):
        seqarr_i = seqarr[i:i + kmer_len]
        scores = np.array([float(substitution_matrix_df.loc[k, s]) for k, s in zip(kmerarr, seqarr_i)])
        score_list.append(np.sum(scores))
    return score_list

def score_kmer_2_seq_no_gaps_v2(
    kmer: str,
    sequence: str,
    substitution_matrix_df: pd.DataFrame,
):
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")

    kmer_len = len(kmer)
    sequence_len = len(sequence)
    kmerarr = np.array(list(kmer))
    seqarr = np.array(list(sequence))
    # maxes = [max]
    scores = [np.trace(substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]].values) for i in range(sequence_len - kmer_len + 1)] # type: ignore
    return scores


def score_kmer_2_seq_no_gaps_v3(
    kmer: str,
    sequence: str,
    substitution_matrix_df: pd.DataFrame,
):
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")

    kmer_len = len(kmer)
    sequence_len = len(sequence)
    kmerarr = np.array(list(kmer))
    seqarr = np.array(list(sequence))
    # maxes = [max]
    # counts = 0
    # scores = []
    # for i in range(sequence_len - kmer_len + 1):
    #     scores.append(np.trace(substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]].values))
    #     counts += 1
    # print(counts)
    # return scores
    scores = np.sum([np.diag(substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]].values) for i in range(sequence_len - kmer_len + 1)], axis=1) # type: ignore
    return scores






# %%



def main2(
    input_alignment_file: str|Path,
    reference_id: str,
    idr_aln_st: int,
    idr_aln_end: int,
    k: int,
    output_file: str|Path,
    scoring_matrix_name: str="grantham_similarity_norm",
):
    fasta_importer = tools.FastaImporter(input_alignment_file)
    aln = fasta_importer.import_as_alignment()
    idr_aln = aln[:, idr_aln_st : idr_aln_end+1]
    idrs = tools.strip_dashes_from_sequences(list(idr_aln)) # type: ignore
    idrs = {i.id: i for i in idrs}
    ref_idr = idrs.pop(reference_id)
    idrs = list(idrs.values())
    score_df, subseq_df, rbm_df = run_pairwise_kmer_alignmentv2(
        str(ref_idr.seq),
        idrs, # type: ignore
        k,
        scoring_matrix_name,
    )
    # output_dict = {
    #     "score_dataframe": score_df.to_dict(orient="split"),
    #     "subseq_dataframe": subseq_df.to_dict(orient="split"),
    #     "reciprocal_best_match_dataframe": rbm_df.to_dict(orient="split"),
    # }
    # with open(output_file, "w") as json_file:
    #     json.dump(output_dict, json_file, indent=4)
# ==============================================================================



