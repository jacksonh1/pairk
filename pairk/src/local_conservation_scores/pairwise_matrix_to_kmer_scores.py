# %%
import concurrent.futures
import json
import multiprocessing
import time
from pathlib import Path
from typing import Callable

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_conservation_scores.tools.capra_singh_2007_scores as cs
import local_seqtools.general_utils as tools
import numpy as np
import pandas as pd
from local_conservation_scores.tools import pairwise_tools
from collections import namedtuple
from attrs import asdict, define, field, validators
from alfpy import word_distance, word_pattern, word_vector
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords as alf_seqrecords



@define
class PairwiseScoreResults:
    hit_sequence: str
    hit_scores: list[float]
    hit_z_scores: list[float]
    flank_hit_sequence: str | None = None
    flank_hit_scores: list[float] | None = None
    flank_hit_z_scores: list[float] | None = None
    background_scores: list[float] | None = None


def kmer_distance_matrix(kmer_list, word_size=2):
    alf_seq_records = alf_seqrecords.SeqRecords(id_list=kmer_list, seq_list=kmer_list)
    p = word_pattern.create(alf_seq_records.seq_list, word_size=word_size)
    # counts = word_vector.Counts(alf_seq_records.length_list, p)
    freqs = word_vector.Freqs(alf_seq_records.length_list, p)
    dist = word_distance.Distance(freqs, 'google')
    matrix = distmatrix.create(alf_seq_records.id_list, dist)
    return matrix

def _alfpy_query_matrix(
    refid: str, matrix
) -> tuple[list[str], list[float]]:
    """
    get the row from the alfpy matrix that corresponds to the query gene
    """
    query_row = [c for c, i in enumerate(matrix.id_list) if refid in i][0]
    query_row_distance = matrix.data[query_row]
    query_row_similarity = 1 - query_row_distance
    return matrix.id_list, query_row_similarity


def filter_similar_kmers(kmer_list: list[str], hit_seq:str, similarity_threshold: float=0.5):
    distmat=kmer_distance_matrix(kmer_list, word_size=2)
    id_list, sim = _alfpy_query_matrix(hit_seq, distmat)
    passing_kmers = [i for i, s in zip(id_list, sim) if s <= similarity_threshold]
    if hit_seq not in passing_kmers:
        passing_kmers.append(hit_seq)
    return passing_kmers


def list_of_strings_to_list_of_columns(seqlist: list[str]):
    seqlist_filtered = [s for s in seqlist if isinstance(s, str)]
    col_strings = []
    for c in range(len(seqlist_filtered[0])):
        col = "".join([seqlist_filtered[i][c] for i in range(len(seqlist_filtered))])
        col_strings.append(col)
    return col_strings


def score_pseudo_aln_columns(
    seqlist: list[str],
    score_func: Callable = cs.property_entropy,
):
    col_strings = list_of_strings_to_list_of_columns(seqlist)
    scores = []
    for c in col_strings:
        scores.append(score_func(c))
    return scores


def subseq_df_2_background_scores(
    subseq_df: pd.DataFrame,
    score_func: Callable = cs.property_entropy,
):
    background_scores = []
    for position in subseq_df.index:
        background_scores.extend(
            score_pseudo_aln_columns(
                subseq_df.loc[position, :].to_list(), score_func=score_func
            )
        )
    return background_scores


def matrix_json_2_pairwise_scores(
    matrix_json: str | Path,
    hit_position: int,
    columnwise_score_func: Callable = cs.property_entropy,
    reciprocal_best_match: bool = False,
    bg_cutoff: int = 50,
    bg_kmer_cutoff: int = 10,
    similarity_threshold: float = 0.5,
) -> PairwiseScoreResults:
    '''
    This is the function that you would want to run if you did it standalone
    '''
    matrix_dict = pairwise_tools.import_pairwise_matrices(matrix_json)
    subseq_df = matrix_dict["subseq_dataframe"].copy()
    k=len(subseq_df.loc[hit_position, "reference_kmer"])
    subseq_df = subseq_df.fillna("-"*k)
    if 'reciprocal_best_match_dataframe' in matrix_dict:
        rbm_df = matrix_dict["reciprocal_best_match_dataframe"].copy()
    if reciprocal_best_match:
        rbm_orthoids = list(
            rbm_df.loc[hit_position][rbm_df.loc[hit_position] == True].index
        )
        if len(rbm_orthoids) == 0:
            raise ValueError("no reciprocal best match found")
        subseq_df = subseq_df[["reference_kmer"] + rbm_orthoids]
    if similarity_threshold < 1:
        # remove kmers from background that are too similar to the hit
        # function input: kmer_list (will be querykmers=seqdf['reference_kmer'].to_list())
        querykmers=subseq_df['reference_kmer'].to_list()
        bg_kmers = filter_similar_kmers(querykmers, subseq_df.loc[hit_position, "reference_kmer"], similarity_threshold=similarity_threshold)
        # if len(bg_kmers) < bg_kmer_cutoff:
        #     raise ValueError(f"not enough background kmers: {len(bg_kmers)}")
        subseq_df=subseq_df[subseq_df['reference_kmer'].isin(bg_kmers)]
    if len(subseq_df) < bg_kmer_cutoff:
        raise ValueError(f"not enough kmers in the matrices: {len(subseq_df)}. Need at least {bg_kmer_cutoff}")
    try:
        background_scores = subseq_df_2_background_scores(subseq_df, score_func=columnwise_score_func)
    except ZeroDivisionError as e:
        print(f"ZeroDivisionError: {matrix_json}")
        raise ValueError(f"ZeroDivisionError: {e}")
    hit_scores = score_pseudo_aln_columns(
        subseq_df.loc[hit_position, :].to_list(), score_func=columnwise_score_func
    )
    if len(background_scores) < bg_cutoff:
        raise ValueError(f"not enough background scores: {len(background_scores)}")
    hit_z_scores = tools.z_score_comparison(hit_scores, background_scores)
    scores = PairwiseScoreResults(
        hit_sequence=subseq_df.loc[hit_position, "reference_kmer"],
        hit_scores=hit_scores,
        hit_z_scores=hit_z_scores,
        background_scores=background_scores,
    )
    return scores



