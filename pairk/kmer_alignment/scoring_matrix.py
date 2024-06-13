import json
from pathlib import Path
import numpy as np
import pandas as pd
import pairk.tools.sequence_utils as tools
import pairk.tools.matrices as matrices
import pairk.tools.pairwise_tools as pairwise_tools
import copy
import pairk.exceptions as _exceptions


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
    '''
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
    '''
    best_scores = {}
    best_subseqs = {}
    best_positions = {}
    for ortholog_id, ortholog_idr in ortholog_idrs.items():
        if len(ortholog_idr) < len(kmer):
            best_subseqs[ortholog_id] = "-" * len(kmer)
            continue
        try:
            ref_k_scores = score_kmer_2_seq_dict_dict(
                kmer, ortholog_idr, score_matrix_dict
            )
            best_score, best_subseq, best_pos = best_from_scores(
                ortholog_idr, ref_k_scores, len(kmer)
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
    ref_idr: str,
    ortholog_idrs: dict[str, str],
    k: int,
    score_matrix_dict: dict[str, dict[str, float | int]],
):
    kmers = tools.gen_kmers(ref_idr, k)
    positions = list(range(len(kmers)))

    score_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    subseq_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    pos_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    rbm_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
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
                ref_k_scores = score_kmer_2_seq_dict_dict(
                    kmer, ortholog_idr, score_matrix_dict
                )
                best_score, best_subseq, best_pos = best_from_scores(
                    ortholog_idr, ref_k_scores, k
                )
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
                rec_scores = score_kmer_2_seq_dict_dict(
                    best_subseq, ref_idr, score_matrix_dict
                )
                reci_best_score, reci_best_subseq, _ = best_from_scores(
                    ref_idr, rec_scores, k
                )
            except ValueError as e:
                print(e)
                continue
            except KeyError as e:
                print(e)
                continue
            rbm_df.loc[position, ortholog_id] = reci_best_subseq == kmer
    return score_df, subseq_df, pos_df, rbm_df


def pairk_alignment(
    idr_dict_in: dict[str, str],
    reference_id: str,
    k: int,
    matrix_name: str = "EDSSMat50",
) -> dict[str, pd.DataFrame]:
    """run pairwise k-mer alignment method using an exhaustive comparison of k-mers. Each reference k-mer is scored against each ortholog k-mer to find the best matching ortholog k-mer in each ortholog.

    Parameters
    ----------
    idr_dict_in : dict[str, str]
        input sequences in dictionary format with the key being the sequence id and the value being the sequence as a string
    reference_id : str
        the id of the reference sequence within the `idr_dict_in` dictionary
    k : int
        the length of the k-mers to use for the alignment
    matrix_name : str, optional
        The name of the scoring matrix to use in the algorithm, by default "EDSSMat50". The available matrices can be viewed with the function `print_available_matrices()` in `pairk.tools.matrices`.

    Returns
    -------
    dict[str, pd.DataFrame]
        results of the pairwise alignment in dictionary format, where the keys are the names of the dataframes and the values are the dataframes. All dataframes have the same structure. One column is the reference k-mer sequence ('reference_kmer'). The other columns are named as the ortholog sequence ids. The dataframe indexes are the reference k-mer start position in the reference sequence. The returned dataframes are:\n
        - 'score_dataframe': the alignment scores for each k-mer in the reference sequence against the corresponding best matching ortholog k-mer.
        - 'subseq_dataframe': the best scoring k-mer from each ortholog for each reference k-mer.
        - 'position_dataframe': the start position of the best scoring k-mer from each ortholog for each reference k-mer.
        - 'reciprocal_best_match_dataframe': a boolean dataframe indicating whether the reference k-mer is the reciprocal best scoring k-mer to the ortholog k-mer.

    """
    _exceptions.validate_matrix_name(matrix_name)
    _exceptions.check_refid_in_idr_dict(idr_dict_in, reference_id)
    idr_dict = copy.deepcopy(idr_dict_in)
    ref_idr = idr_dict.pop(reference_id)
    scoremat_df = matrices.load_matrix_as_df(matrix_name)
    scoremat_dict = matrices.matrix_df_to_dict(scoremat_df)
    score_df, subseq_df, pos_df, rbm_df = run_pairwise_kmer_alignment(
        ref_idr, idr_dict, k, score_matrix_dict=scoremat_dict
    )
    output_dict = {
        "score_dataframe": score_df,
        "subseq_dataframe": subseq_df,
        "position_dataframe": pos_df,
        "reciprocal_best_match_dataframe": rbm_df,
    }
    return output_dict
