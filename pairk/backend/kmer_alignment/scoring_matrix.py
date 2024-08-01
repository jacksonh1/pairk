import json
from pathlib import Path
import numpy as np
import pandas as pd
import pairk.backend.tools.sequence_utils as tools
import pairk.backend.tools.matrices as matrices
import pairk.backend.tools.pairwise_tools as pairwise_tools
import pairk.backend.exceptions as _exceptions
import copy


def score_alignment(
    seq1: str, seq2: str, substitution_matrix_dict: dict[str, dict[str, float | int]]
) -> float:
    assert len(seq1) == len(seq2)
    score = 0
    for s1, s2 in zip(seq1, seq2):
        score += float(substitution_matrix_dict[s1][s2])  # type: ignore
    return score


def score_kmer_2_seq_dict_dict(
    kmer: str,
    sequence: str,
    substitution_matrix_dict: dict[str, dict[str, float | int]],
):
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")
    scores = []
    for i in range(len(sequence) - len(kmer) + 1):
        subseq = sequence[i : i + len(kmer)]
        score = score_alignment(kmer, subseq, substitution_matrix_dict)
        scores.append(score)
    return scores


def best_from_scores(sequence: str, scores: list[float], k):
    best_match_pos = np.argmax(scores)
    bestscore = scores[best_match_pos]
    best_subseq = sequence[best_match_pos : best_match_pos + k]
    return bestscore, best_subseq, int(best_match_pos)


def run_pairwise_kmer_alignment(
    query_idr: str,
    ortholog_idrs: dict[str, str],
    k: int,
    score_matrix_dict: dict[str, dict[str, float | int]],
):
    kmers = tools.gen_kmers(query_idr, k)
    positions = list(range(len(kmers)))

    score_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    orthokmer_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    pos_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    score_df.loc[positions, "query_kmer"] = kmers
    orthokmer_df.loc[positions, "query_kmer"] = kmers
    pos_df.loc[positions, "query_kmer"] = kmers
    for position, kmer in zip(positions, kmers):
        for ortholog_id, ortholog_idr in ortholog_idrs.items():
            if len(ortholog_idr) < k:
                orthokmer_df.loc[position, ortholog_id] = "-" * k
                continue
            query_k_scores = score_kmer_2_seq_dict_dict(
                kmer, ortholog_idr, score_matrix_dict
            )
            best_score, best_subseq, best_pos = best_from_scores(
                ortholog_idr, query_k_scores, k
            )
            score_df.loc[position, ortholog_id] = best_score
            orthokmer_df.loc[position, ortholog_id] = best_subseq
            pos_df.loc[position, ortholog_id] = best_pos
    return score_df, orthokmer_df, pos_df


def pairk_alignment(
    idr_dict: dict[str, str],
    query_id: str,
    k: int,
    matrix_name: str = "EDSSMat50",
) -> pairwise_tools.PairkAln:
    """run pairwise k-mer alignment method using an exhaustive comparison of k-mers.
    Each query k-mer is scored against each ortholog k-mer to find the best matching
    ortholog k-mer in each ortholog. If a ortholog IDR is shorter than the k-mer, a
    string of "-" characters ("-"\\*k) is assigned as the best matching ortholog k-mer for that
    ortholog.

    **Note**: if there are multiple top-scoring matches, only one is returned.

    Parameters
    ----------
    idr_dict : dict[str, str]
        input sequences in dictionary format with the key being the sequence id and
        the value being the sequence as a string
    query_id : str
        the id of the query sequence within the `idr_dict` dictionary
    k : int
        the length of the k-mers to use for the alignment
    matrix_name : str, optional
        The name of the scoring matrix to use in the algorithm, by default "EDSSMat50".
        The available matrices can be viewed with the function `pairk.print_available_matrices()`.

    Returns
    -------
    pairwise_tools.PairkAln
        an object containing the alignment results. See the `pairk.PairkAln` class for more information.
    """
    _exceptions.validate_matrix_name(matrix_name)
    _exceptions.check_queryid_in_idr_dict(idr_dict, query_id)
    scoremat_df = matrices.load_matrix_as_df(matrix_name)
    scoremat_dict = matrices.matrix_df_to_dict(scoremat_df)
    idr_str_dict = copy.deepcopy(idr_dict)
    _exceptions.check_sequence_characters_dict(idr_str_dict, matrix_name)
    query_idr = idr_str_dict.pop(query_id)
    score_df, orthokmer_df, pos_df = run_pairwise_kmer_alignment(
        query_idr, idr_str_dict, k, score_matrix_dict=scoremat_dict
    )
    return pairwise_tools.PairkAln(
        orthokmer_df=orthokmer_df,
        pos_df=pos_df,
        score_df=score_df,
    )


def pairk_alignment_single_kmer(
    kmer: str,
    ortholog_idrs: dict[str, str],
    matrix_name: str = "EDSSMat50",
):
    """
    Align a single kmer to a dictionary of sequences and return the best scoring subsequence for each sequence in the dictionary.

    Parameters
    ----------
    kmer : str
        the kmer to align
    ortholog_idrs : dict[str, str]
        a dictionary of sequences to align the kmer to, with the key being the sequence id and the value being the sequence as a string
    matrix_name : str, optional
        The name of the scoring matrix to use in the algorithm, by default "EDSSMat50".
        The available matrices can be viewed with the function `print_available_matrices()`
        in `pairk.backend.tools.matrices`.

    Returns
    -------
    dict[str, float], dict[str, str], dict[str, int]
        the best scores, best subsequences, and best positions for the kmer in each sequence in the input `ortholog_idrs` dictionary
    """
    _exceptions.validate_matrix_name(matrix_name)
    # I don't like how innefficient this is, but I find the normal error for
    # invalid characters to be confusing
    check_dict = copy.deepcopy(ortholog_idrs)
    check_dict["kmer"] = kmer
    _exceptions.check_sequence_characters_dict(check_dict, matrix_name)
    scoremat_df = matrices.load_matrix_as_df(matrix_name)
    scoremat_dict = matrices.matrix_df_to_dict(scoremat_df)
    best_scores = {}
    best_subseqs = {}
    best_positions = {}
    for ortholog_id, ortholog_idr in ortholog_idrs.items():
        if len(ortholog_idr) < len(kmer):
            best_subseqs[ortholog_id] = "-" * len(kmer)
            continue
        query_k_scores = score_kmer_2_seq_dict_dict(kmer, ortholog_idr, scoremat_dict)
        best_score, best_subseq, best_pos = best_from_scores(
            ortholog_idr, query_k_scores, len(kmer)
        )
        if query_k_scores.count(best_score) > 1:
            print(f"Warning: multiple best scores found for {ortholog_id}")
        best_scores[ortholog_id] = best_score
        best_subseqs[ortholog_id] = best_subseq
        best_positions[ortholog_id] = best_pos
    return best_scores, best_subseqs, best_positions


# def reciprocal_best_match(query_idr, query_kmer, ortho_kmers, matrix_name):
#     """
#     Calculate the reciprocal best match (RBM) of a query kmer to a set of ortholog kmers.
#     The RBM is calculated by aligning the best scoring ortholog kmer to the query sequence and checking if the best scoring query kmer is the same as the original query kmer.

#     Parameters
#     ----------
#     query_idr : str
#         the query sequence
#     query_kmer : str
#         the query kmer
#     ortho_kmers : dict[str, str]
#         a dictionary of ortholog kmers with the key being the ortholog id and the value being the ortholog kmer that originally aligned to the query kmer
#     matrix_name : str
#         the name of the scoring matrix to use in the alignment

#     Returns
#     -------
#     dict[str, bool]
#         a dictionary of boolean values indicating whether the query kmer ortholog kmer match is reciprocal
#     """
#     _exceptions.validate_matrix_name(matrix_name)
#     _exceptions.check_sequence_characters_list(
#         [query_idr] + list(ortho_kmers.values()), matrix_name
#     )
#     scoremat_df = matrices.load_matrix_as_df(matrix_name)
#     scoremat_dict = matrices.matrix_df_to_dict(scoremat_df)
#     rbm_dict = {}
#     for ortholog_id, ortholog_kmer in ortho_kmers.items():
#         try:
#             rec_scores = score_kmer_2_seq_dict_dict(
#                 ortholog_kmer, query_idr, scoremat_dict
#             )
#             reci_best_score, reci_best_subseq, _ = best_from_scores(
#                 query_idr, rec_scores, len(query_kmer)
#             )
#         except ValueError as e:
#             print(e)
#             continue
#         except KeyError as e:
#             print(e)
#             continue
#         rbm_dict[ortholog_id] = reci_best_subseq == query_kmer
#     return rbm_dict
