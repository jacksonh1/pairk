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
    """
    best_score = np.max(scores)
    best_position = np.argmax(scores)
    best_subseq = sequence[best_position : best_position + len(kmer)]
    """
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


def pairk_alignment_single_kmer(
    kmer: str,
    ortholog_idrs: dict[str, str],
    score_matrix_dict: dict[str, dict[str, float | int]],
):
    """
    Align a kmer to a dictionary of sequences and return the best scoring subsequence for each sequence in the dictionary.

    Parameters
    ----------
    kmer : str
        the kmer to align
    ortholog_idrs : dict[str, str]
        a dictionary of sequences to align the kmer to, with the key being the sequence id and the value being the sequence as a string
    score_matrix_dict : dict[str, dict[str, float | int]]
        the scoring matrix to use for the alignment

    Returns
    -------
    dict[str, float], dict[str, str], dict[str, int]
        the best scores, best subsequences, and best positions for the kmer in each sequence in the input `ortholog_idrs` dictionary
    """
    best_scores = {}
    best_subseqs = {}
    best_positions = {}
    for ortholog_id, ortholog_idr in ortholog_idrs.items():
        if len(ortholog_idr) < len(kmer):
            best_subseqs[ortholog_id] = "-" * len(kmer)
            continue
        try:
            query_k_scores = score_kmer_2_seq_dict_dict(
                kmer, ortholog_idr, score_matrix_dict
            )
            best_score, best_subseq, best_pos = best_from_scores(
                ortholog_idr, query_k_scores, len(kmer)
            )
        except ValueError as e:
            best_subseqs[ortholog_id] = "-" * len(kmer)
            print(e)
            continue
        except KeyError as e:
            best_subseqs[ortholog_id] = "-" * len(kmer)
            print(e)
            continue
        best_scores[ortholog_id] = best_score
        best_subseqs[ortholog_id] = best_subseq
        best_positions[ortholog_id] = best_pos
    return best_scores, best_subseqs, best_positions


def run_pairwise_kmer_alignment(
    query_idr: str,
    ortholog_idrs: dict[str, str],
    k: int,
    score_matrix_dict: dict[str, dict[str, float | int]],
    rbm: bool = False,
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
    rbm_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    score_df.loc[positions, "query_kmer"] = kmers
    orthokmer_df.loc[positions, "query_kmer"] = kmers
    pos_df.loc[positions, "query_kmer"] = kmers
    if rbm:
        rbm_df.loc[positions, "query_kmer"] = kmers
    else:
        rbm_df = None
    for position, kmer in zip(positions, kmers):
        for ortholog_id, ortholog_idr in ortholog_idrs.items():
            if len(ortholog_idr) < k:
                orthokmer_df.loc[position, ortholog_id] = "-" * k
                continue
            try:
                query_k_scores = score_kmer_2_seq_dict_dict(
                    kmer, ortholog_idr, score_matrix_dict
                )
                best_score, best_subseq, best_pos = best_from_scores(
                    ortholog_idr, query_k_scores, k
                )
            except ValueError as e:
                orthokmer_df.loc[position, ortholog_id] = "-" * k
                print(e)
                continue
            except KeyError as e:
                orthokmer_df.loc[position, ortholog_id] = "-" * k
                print(e)
                continue
            score_df.loc[position, ortholog_id] = best_score
            orthokmer_df.loc[position, ortholog_id] = best_subseq
            pos_df.loc[position, ortholog_id] = best_pos
            if rbm:
                try:
                    rec_scores = score_kmer_2_seq_dict_dict(
                        best_subseq, query_idr, score_matrix_dict
                    )
                    reci_best_score, reci_best_subseq, _ = best_from_scores(
                        query_idr, rec_scores, k
                    )
                except ValueError as e:
                    print(e)
                    continue
                except KeyError as e:
                    print(e)
                    continue
                rbm_df.loc[position, ortholog_id] = reci_best_subseq == kmer  # type: ignore
    return score_df, orthokmer_df, pos_df, rbm_df


def pairk_alignment(
    idr_dict_in: dict[str, str],
    query_id: str,
    k: int,
    matrix_name: str = "EDSSMat50",
    rbm: bool = False,
) -> pairwise_tools.PairkAln:
    """run pairwise k-mer alignment method using an exhaustive comparison of k-mers. Each query k-mer is scored against each ortholog k-mer to find the best matching ortholog k-mer in each ortholog.

    Parameters
    ----------
    idr_dict_in : dict[str, str]
        input sequences in dictionary format with the key being the sequence id and the value being the sequence as a string
    query_id : str
        the id of the query sequence within the `idr_dict_in` dictionary
    k : int
        the length of the k-mers to use for the alignment
    matrix_name : str, optional
        The name of the scoring matrix to use in the algorithm, by default "EDSSMat50".
        The available matrices can be viewed with the function `print_available_matrices()`
        in `pairk.backend.tools.matrices`.
    rbm : bool, optional
        whether to calculate the reciprocal best match (RBM) for each k-mer, by default False.
        It true, the RBM matrix will be calculated by aligning the best scoring
        ortholog k-mer to the query sequence and checking if the best
        scoring query k-mer is the same as the original query k-mer.
        This is currently unused in the pairk package. The RBM matrix will be
        included in the output `PairkAln` object, it is a boolean dataframe
        indicating whether the query k-mer ortholog k-mer match is reciprocal.

    Returns
    -------
    pairwise_tools.PairkAln
        an object containing the alignment results. See the `pairk.PairkAln` class for more information.
    """
    _exceptions.validate_matrix_name(matrix_name)
    _exceptions.check_queryid_in_idr_dict(idr_dict_in, query_id)
    idr_dict = copy.deepcopy(idr_dict_in)
    query_idr = idr_dict.pop(query_id)
    scoremat_df = matrices.load_matrix_as_df(matrix_name)
    scoremat_dict = matrices.matrix_df_to_dict(scoremat_df)
    score_df, orthokmer_df, pos_df, rbm_df = run_pairwise_kmer_alignment(
        query_idr, idr_dict, k, score_matrix_dict=scoremat_dict, rbm=rbm
    )
    return pairwise_tools.PairkAln(
        orthokmer_df=orthokmer_df,
        pos_df=pos_df,
        score_df=score_df,
        rbm_df=rbm_df,
    )
