'''
exhaustive_pairk_align:
    - reference_id: str    
    - idr_dictionary: dict[str, str]
    - k: int
    - matrix: dict[str, dict[str, float]]


needleman_pairk_align:
    - reference_id: str    
    - idr_dictionary: dict[str, str]
    - k: int
    - aligner: Bio.Align.PairwiseAligner


input set 1:
    - reference_id: str
    - MSA: list[str]
    - idr_aln_st: int
    - idr_aln_end: int
    - k: int
    - alignment_type: str = ("exhaustive", "dynamic_programming")
    - matrix_name: str
        - or-
    - matrix_file: str = "grantham_similarity_normx100_aligner_compatible"
'''

from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
import copy

from local_env_variables import matrices
from local_seqtools import general_utils as tools

def get_aligner(
    matrix: Align.substitution_matrices.Array, 
) -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -100000
    aligner.extend_gap_score = -100000
    aligner.query_end_gap_score = 0.0
    aligner.mode = "global"
    aligner.substitution_matrix = matrix
    # print('getting aligner')
    return aligner


def make_empty_kmer_ortho_df(positions, ortholog_ids: list[str]):
    cols = ["reference_kmer"] + ortholog_ids
    df = pd.DataFrame(
        index=positions,
        columns=cols,
    )
    return df


def kmer_align_aligner(
    kmer: str,
    sequence: str,
    aligner: Align.PairwiseAligner, 
) -> Align.Alignment:
    alignment = aligner.align(Seq.Seq(sequence), Seq.Seq(kmer))[0]
    return alignment


def first_non_dash_index(s):
    index = 0
    while index < len(s) and s[index] == '-':
        index += 1
    return index


def score_kmer_2_seq_no_gaps_needleman(
    kmer: str,
    sequence: str,
    aligner: Align.PairwiseAligner,
):
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")
    aln = kmer_align_aligner(kmer, sequence, aligner=aligner)
    # bseq, bpos = get_aas_aligned_2_query_nongap_positions(aln)
    # bpos = get_position_of_query_start_in_target(aln)
    s = first_non_dash_index(aln[1])
    best_subseq = aln[0][s:s+len(kmer)]
    best_pos = s
    return aln.score, best_subseq, best_pos




def needleman_pairwise_kmer_alignment(
    ref_idr: str,
    ortholog_idrs: dict[str, str],
    k: int,
    aligner: Align.PairwiseAligner,
):
    # substitution_matrix_df = matrices.load_precomputed_matrix_df(scoring_matrix_name)
    kmers = tools.gen_kmers(ref_idr, k)
    positions = list(range(len(kmers)))

    score_df = make_empty_kmer_ortho_df(positions, list(ortholog_idrs.keys()))
    subseq_df = make_empty_kmer_ortho_df(positions, list(ortholog_idrs.keys()))
    pos_df = make_empty_kmer_ortho_df(positions, list(ortholog_idrs.keys()))
    rbm_df = make_empty_kmer_ortho_df(positions, list(ortholog_idrs.keys()))
    score_df.loc[positions, "reference_kmer"] = kmers
    subseq_df.loc[positions, "reference_kmer"] = kmers
    pos_df.loc[positions, "reference_kmer"] = kmers
    rbm_df.loc[positions, "reference_kmer"] = kmers
    for position, kmer in zip(positions, kmers):
        for ortholog_id, ortholog_idr in ortholog_idrs.items():
            if len(ortholog_idr) < k:
                subseq_df.loc[position, ortholog_id] = "-" * k
                continue
            try:
                best_score, best_subseq, best_pos = score_kmer_2_seq_no_gaps_needleman(kmer, ortholog_idr, aligner)
            except ValueError as e:
                subseq_df.loc[position, ortholog_id] = "-" * k
                print(e)
                continue
            except KeyError as e:
                subseq_df.loc[position, ortholog_id] = "-" * k
                print(e)
                continue
            score_df.loc[position, ortholog_id] = best_score
            subseq_df.loc[position, ortholog_id] = best_subseq
            pos_df.loc[position, ortholog_id] = best_pos
            try:
                reci_best_score, reci_best_subseq, _ = score_kmer_2_seq_no_gaps_needleman(best_subseq, ref_idr, aligner)
            except ValueError as e:
                print(e)
                continue
            except KeyError as e:
                print(e)
                continue
            rbm_df.loc[position, ortholog_id] = (reci_best_subseq==kmer)
    return score_df, subseq_df, pos_df, rbm_df


def main(
    idr_dict_in: dict[str, str],
    reference_id: str,
    k: int,
    matrix_name: str="grantham_similarity_norm",
):
    idr_dict = copy.deepcopy(idr_dict_in)
    ref_idr = idr_dict.pop(reference_id)
    matrix=matrices.load_matrix_for_aligner(matrix_name)
    aligner = get_aligner(matrix)
    score_df, subseq_df, pos_df, rbm_df = needleman_pairwise_kmer_alignment(
        ref_idr,
        idr_dict,
        k,
        aligner=aligner,
    )
    output_dict = {
        "score_dataframe": score_df.to_dict(orient="split"),
        "subseq_dataframe": subseq_df.to_dict(orient="split"),
        "position_dataframe": pos_df.to_dict(orient="split"),
        "reciprocal_best_match_dataframe": rbm_df.to_dict(orient="split"),
    }
    return output_dict


