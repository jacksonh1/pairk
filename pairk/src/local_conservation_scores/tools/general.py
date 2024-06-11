import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO

import local_env_variables.matrices as mats


def aln_gaps_mask(aln, gap_frac_cutoff=0.20):
    """returns a mask of positions in alignment with less than gap_perc_cutoff gaps"""
    gap_frac_list = []
    alignment_width = aln.get_alignment_length()
    n_seqs_in_aln = len(aln)
    for i in range(alignment_width):
        gap_count = aln[:, i].count("-")
        gap_frac = gap_count / n_seqs_in_aln
        gap_frac_list.append(gap_frac)
    return np.array(gap_frac_list) < gap_frac_cutoff


def add_seq_gap_mask_to_mask(
    score_list_seq_str: str,
    gappy_mask: np.ndarray,
):
    assert len(score_list_seq_str) == len(gappy_mask)
    master_mask = []
    for s, m in zip(score_list_seq_str, gappy_mask):
        if s == "-":
            master_mask.append(False)
        else:
            master_mask.append(m)
    return master_mask


def make_score_mask(alignment: Align.MultipleSeqAlignment, query_seq_str_in_aln:str, gap_frac_cutoff=0.20):
    """make gap and score mask for an alignment

    Parameters
    ----------
    alignment : Align.MultipleSeqAlignment
        multiple sequence alignment
    query_seq_str_in_aln : str
        query sequence string in the alignment (including gaps)
    gap_frac_cutoff : float, optional
        if a position in the msa has < `gap_frac_cutoff`, then the corresponding position in the gap mask will be True, otherwise it will be set to False. by default 0.20

    Returns
    -------
    list
        [gap_mask, score_mask]. Both are boolean np.arrays. The gap mask is True for positions in the alignment with less than `gap_frac_cutoff` gaps. The score mask is True for positions in the alignment with less than `gap_frac_cutoff` gaps, and not a gap in the query sequence
    """    
    gap_mask = aln_gaps_mask(
        alignment,
        gap_frac_cutoff=gap_frac_cutoff,
    )
    score_mask = add_seq_gap_mask_to_mask(
        query_seq_str_in_aln,
        gap_mask,
    )
    return np.array(gap_mask), np.array(score_mask)


def mut_valdar_df_mat(a, b, mat_df):
    if a == "-" or b == "-":
        return 0
    if a not in mat_df.index or b not in mat_df.index:
        return 0
    return mat_df.loc[a, b]


def mut_pid(a, b):
    if a == "-" or b == "-":
        return 0
    if a == b:
        return 1
    return 0


def asymmetric_valdar_score_df_unweighted(
    query_alignment_record, alignment, mat_df, score_function=mut_valdar_df_mat
):
    scores = []
    for pos_i in range(len(query_alignment_record.seq)):
        score_i = 0
        for seq_j in alignment:
            if seq_j.id == query_alignment_record.id:
                continue
            score_i += score_function(
                query_alignment_record.seq[pos_i], seq_j.seq[pos_i], mat_df
            )
        scores.append(score_i / (len(alignment) - 1))
    return scores


def asymmetric_valdar_score_df_unweighted_strings(
    query_alignment_str,
    query_id,
    alignment,
    mat_df=mats.load_precomputed_matrix_df(matrix_name='EDSSMat50_max_off_diagonal_norm'),
):
    scores = []
    found_query = False
    for pos_i in range(len(query_alignment_str)):
        score_i = 0
        for seq_j in alignment:
            if seq_j.id == query_id:
                found_query = True
                continue
            score_i += mut_valdar_df_mat(
                query_alignment_str[pos_i], seq_j.seq[pos_i], mat_df
            )
        scores.append(score_i / (len(alignment) - 1))
    if not found_query:
        raise ValueError("query_id not found in alignment")
    return scores
    
    
def symmetric_valdar_score_df_unweighted(alignment, mat_df):
    scores = []
    for pos_i in range(len(alignment[0])):
        score_i = 0
        l = 0
        for j in range(len(alignment)):
            for k in range(j + 1, len(alignment)):
                score_i += mut_valdar_df_mat(
                    alignment[j][pos_i], alignment[k][pos_i], mat_df
                )
                l += 1
        scores.append(score_i / l)
    return scores


def conservation_string(scores, seq, z_score_cutoff: float = 1.0):
    assert len(scores) == len(seq)
    cons_str = ''
    for score, aa in zip(scores, seq):
        if score > z_score_cutoff:
            cons_str += aa
        else:
            cons_str += '_'
    return cons_str