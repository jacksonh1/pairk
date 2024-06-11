import json
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, SeqIO, SeqRecord

import local_conservation_scores.tools.capra_singh_2007_scores as cs
from local_seqtools import alignment_tools as aln_tools
from local_seqtools import general_utils as tools


def matrix_dict_2_df(matrix_dict, mat_key):
    df = pd.DataFrame(
        matrix_dict[mat_key]["data"],
        columns=matrix_dict[mat_key]["columns"],
        index=matrix_dict[mat_key]["index"],
    )
    return df


def import_pairwise_matrices(
    filepath,
    matrix_keys: list[str] = [
        "score_dataframe",
        "subseq_dataframe",
        "reciprocal_best_match_dataframe",
        "position_dataframe",
    ],
):
    with open(filepath, "r") as json_file:
        data = json.load(json_file)
    matrices = {}
    for k in matrix_keys:
        if k in data:
            matrices[k] = matrix_dict_2_df(data, k)
    if len(matrices) == 0:
        raise ValueError("No matrices found in the json file")
    return matrices


def pad_hit_target_length_approximate(query_sequence: str, st: int, end: int, target_hit_length: int = 36):
    hit_sequence = query_sequence[st : end + 1]
    if len(hit_sequence) >= target_hit_length:
        return st, end, hit_sequence
    flank = round((target_hit_length - len(hit_sequence)) / 2)
    s_i_new = max(0, st - flank)
    s_e_new = min(len(query_sequence), end + flank)
    return s_i_new, s_e_new, query_sequence[s_i_new : s_e_new + 1]


def pad_hit_target_length_exact(query_sequence: str, st: int, end: int, target_hit_length: int = 36):
    hit_sequence = query_sequence[st : end + 1]
    if len(hit_sequence) >= target_hit_length:
        return st, end, hit_sequence
    st_i=st
    end_i=end
    lside=True
    while len(hit_sequence) < target_hit_length:
        if lside:
            st_i = max(0, st_i-1)
            hit_sequence = query_sequence[st_i:end_i+1]
            lside = False
        else:
            end_i = min(len(query_sequence)-1, end_i+1)
            hit_sequence = query_sequence[st_i:end_i+1]
            lside = True
    return st_i, end_i, query_sequence[st_i : end_i + 1]








# def split_query_seq_for_kmer_gen(query_seq: str, hit_start: int, hit_end: int):
#     pt1 = query_seq[:hit_start]
#     pt2 = query_seq[hit_end + 1 :]
#     return pt1, pt2


# def get_kmers(query_seq, hit_st, hit_end):
#     k = hit_end - hit_st + 1
#     pt1, pt2 = split_query_seq_for_kmer_gen(query_seq, hit_st, hit_end)
#     kmers = []
#     kmers.extend(tools.gen_kmers(pt1, k))
#     kmers.extend(tools.gen_kmers(pt2, k))
#     kmers.append(query_seq[hit_st : hit_end + 1])
#     return list(set(kmers))



# def aln_kmer_2_seq_no_gaps_normed(
#     kmer: str,
#     sequence: str,
#     substitution_matrix_df: pd.DataFrame,
# ):
#     positions = []
#     subsequences = []
#     scores = []
#     if len(sequence) < len(kmer):
#         return 0, "-" * len(kmer), -1
#     for i in range(len(sequence) - len(kmer) + 1):
#         subseq = sequence[i : i + len(kmer)]
#         max_score = max(
#             aln_tools.score_alignment(subseq, subseq, substitution_matrix_df),
#             aln_tools.score_alignment(kmer, kmer, substitution_matrix_df),
#         )
#         score = (
#             aln_tools.score_alignment(kmer, subseq, substitution_matrix_df) / max_score
#         )
#         positions.append(i)
#         scores.append(score)
#         subsequences.append(subseq)
#     best_score = np.max(scores)
#     best_subseq = subsequences[np.argmax(scores)]
#     best_pos = positions[np.argmax(scores)]
#     # or
#     # best_subseq = sequence[best_pos : best_pos + len(kmer)]
#     return np.max(scores), subsequences[np.argmax(scores)], positions[np.argmax(scores)]





# def score_alignment_columnwise(
#     aln: Align.MultipleSeqAlignment, score_func: Callable = cs.property_entropy
# ):
#     scores = []
#     for c in range(aln.get_alignment_length()):
#         col = aln[:, c]
#         score = score_func(col)
#         scores.append(score)
#     return scores


# def valdar_column_score(q_res, alignment_column, scoring_matrix_df):
#     score = 0
#     for i in alignment_column:
#         score += cons_tools.mut_valdar_df_mat(q_res, i, scoring_matrix_df)
#     return score / len(alignment_column)


# def score_alignment_columnwise_valdar(
#     aln: Align.MultipleSeqAlignment, query_id, scoring_matrix_df
# ):
#     aln_dict = {i.id: i for i in aln}
#     query_seq = str(aln_dict.pop(query_id).seq)
#     aln2 = AlignIO.MultipleSeqAlignment(list(aln_dict.values()))
#     scores = []
#     for c, r in enumerate(query_seq):
#         col = aln2[:, c]
#         score = valdar_column_score(r, col, scoring_matrix_df)
#         scores.append(score)
#     return scores
