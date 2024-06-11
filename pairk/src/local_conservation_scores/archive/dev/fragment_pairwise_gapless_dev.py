# %%
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
from Bio.SeqRecord import SeqRecord

from local_conservation_scores.tools import capra_singh_2007_scores, general
from local_env_variables import matrices
from local_seqtools import alignment_tools as aln_tools
from local_seqtools import general_utils as tools
from local_seqtools import jch_alignment as jch_aln
from typing import Callable
import time


# %%

def aln_small_str_2_seqs_no_gaps(
    query_seq: str,
    seqrecord_list: list[SeqRecord],
    substitution_matrix_path: str|Path = matrices.MATRIX_DIR / "grantham_similarity_normx100_aligner_compatible", 
) -> dict[str, Align.Alignment]:
    """return list of sequences aligning to query_seq"""
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -100000
    aligner.extend_gap_score = -100000
    aligner.query_end_gap_score = 0.0
    aligner.mode = "global"
    aligner.substitution_matrix = Align.substitution_matrices.read(
        substitution_matrix_path
    )
    # scoring_matrix_name = 'BLOSUM62'
    # aligner.substitution_matrix = Align.substitution_matrices.load(scoring_matrix_name)
    alignments_dict = {}
    for seq in seqrecord_list:
        if len(seq.seq) == 0:
            alignment = Align.Alignment(["-" * len(query_seq), query_seq])
        else:
            alignment = aligner.align(seq.seq, Seq.Seq(query_seq))[0]
        alignments_dict[seq.id] = alignment
    return alignments_dict

def kmer_align_aligner(
    kmer: str,
    sequence: str,
    substitution_matrix_path: str|Path = matrices.MATRIX_DIR / "grantham_similarity_normx100_aligner_compatible", 
) -> Align.Alignment:
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -100000
    aligner.extend_gap_score = -100000
    aligner.query_end_gap_score = 0.0
    aligner.mode = "global"
    aligner.substitution_matrix = Align.substitution_matrices.read(
        substitution_matrix_path
    )
    alignment = aligner.align(Seq.Seq(sequence), Seq.Seq(kmer))[0]
    return alignment


def get_aas_aligned_2_query_nongap_positions(alignment: Align.Alignment):
    query_arr = np.array(list(alignment[1])) # type: ignore
    target_arr = np.array(list(alignment[0])) # type: ignore
    nongap_indices = [c for c, i in enumerate(query_arr) if i != "-"]
    result = target_arr[nongap_indices]
    return ''.join(result)


def kmer_align_aligner_allseqs(
    query_seq: str,
    seqrecord_list: list[SeqRecord],
    substitution_matrix_path: str|Path = matrices.MATRIX_DIR / "grantham_similarity_normx100_aligner_compatible", 
):
    """return list of sequences aligning to query_seq"""
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -100000
    aligner.extend_gap_score = -100000
    aligner.query_end_gap_score = 0.0
    aligner.mode = "global"
    aligner.substitution_matrix = Align.substitution_matrices.read(
        substitution_matrix_path
    )
    # scoring_matrix_name = 'BLOSUM62'
    # aligner.substitution_matrix = Align.substitution_matrices.load(scoring_matrix_name)
    alignments_dict = {}
    score_dict = {}
    seq_dict = {}
    for seq in seqrecord_list:
        alignment = aligner.align(seq.seq, Seq.Seq(query_seq))[0]
        alignments_dict[seq.id] = alignment
        score_dict[seq.id] = alignment.score
        seq_dict[seq.id] = get_aas_aligned_2_query_nongap_positions(alignment)
    return alignments_dict, score_dict, seq_dict


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


def score_kmer_2_seq_no_gaps_v4(
    kmer: str,
    sequence: str,
    substitution_matrix_file = matrices.MATRIX_DIR / "grantham_similarity_normx100_aligner_compatible",
):
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")
    aln = kmer_align_aligner(kmer, sequence, substitution_matrix_path=substitution_matrix_file)
    bseq = get_aas_aligned_2_query_nongap_positions(aln)
    return aln.score, bseq

kmer = 'KQK'
# sequence = 'KKKKKKKK'
sequence = 'ACDEFGHI'
substitution_matrix_df = matrices.load_precomputed_matrix_df('grantham_similarity_norm')
substitution_matrix_df = round(substitution_matrix_df*100)
substitution_matrix_df = substitution_matrix_df.astype(int)

# %%

def time_func(func: Callable, n_iters, **args):
    a = time.time()
    for i in range(n_iters):
        func(**args)
    b = time.time()
    print(f'time: {b-a} seconds')
args ={
    'kmer':kmer,
    'sequence':sequence,
    'substitution_matrix_df':substitution_matrix_df
}
time_func(score_kmer_2_seq_no_gaps_old, 1000, **args)
time_func(score_kmer_2_seq_no_gaps_old2, 1000, **args)
time_func(score_kmer_2_seq_no_gaps_v1, 1000, **args)
time_func(score_kmer_2_seq_no_gaps_v2, 1000, **args)
time_func(score_kmer_2_seq_no_gaps_v3, 1000, **args)
args.pop('substitution_matrix_df')
time_func(score_kmer_2_seq_no_gaps_v4, 1000, **args)

score_kmer_2_seq_no_gaps_v1(kmer, sequence, substitution_matrix_df)
score_kmer_2_seq_no_gaps_v2(kmer, sequence, substitution_matrix_df)
score_kmer_2_seq_no_gaps_v3(kmer, sequence, substitution_matrix_df)
score_kmer_2_seq_no_gaps_v4(kmer, sequence)



# %%


def kmer_best_matches_2_ortho(
    kmer: str,
    sequence: str,
    substitution_matrix_df: pd.DataFrame,
):
    scores = score_kmer_2_seq_no_gaps_v3(kmer, sequence, substitution_matrix_df)
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
    kmer_count = 0
    for kmer in kmers:
        a = time.time()
        for ortholog_id, ortholog_idr in ortholog_idrs.items():
            if len(ortholog_idr) < k:
                continue
            best_score, best_subseq = score_kmer_2_seq_no_gaps_v4(kmer, ortholog_idr)
            # best_score, best_subseq, _ = kmer_best_matches_2_ortho(
            #     kmer,
            #     ortholog_idr,
            #     substitution_matrix_df,
            # )
            score_df.loc[kmer, ortholog_id] = best_score
            subseq_df.loc[kmer, ortholog_id] = best_subseq
            reci_best_score, reci_best_subseq= score_kmer_2_seq_no_gaps_v4(
                best_subseq,
                ref_idr,
            )
            rbm_df.loc[kmer, ortholog_id] = (reci_best_subseq==kmer)
        kmer_count += 1
        b = time.time()
        print(f"Finished kmer: {kmer}\n{kmer_count} of {len(kmers)}")
        print(f"Time: {b-a} seconds")
    return score_df, subseq_df, rbm_df

    
def main(
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
import local_conservation_analysis_pipeline.group_conservation_objects as group_tools

json_file = "../../../../05-conservation_pipeline/examples/table_annotation/conservation_analysis/3-9606_0_002f40/3-9606_0_002f40.json"
og = group_tools.ConserGene(json_file)
lvlo=og.get_level_obj('Vertebrata')

def pad_hit(seq: str, st_pos: int, end_pos: int, l_flank: int = 0, r_flank: int = 0):
    st = max(0, st_pos - l_flank)
    end = min(len(seq)-1, end_pos + r_flank)
    return st, end, seq[st : end + 1]

_,_,flanked_hit = pad_hit(og.query_sequence, og.hit_start_position, og.hit_end_position, 5, 5)
k = len(flanked_hit)
og.hit_sequence
main(
    input_alignment_file=lvlo.alignment_file,
    reference_id=og.query_gene_id,
    idr_aln_st=lvlo.idr_aln_start,
    idr_aln_end=lvlo.idr_aln_end,
    k=k,
    output_file='test.json'
)































# %%

def run_pairwise_kmer_alignmentv2(
    ref_idr: str,
    ortholog_idrs: list[SeqRecord],
    k: int,
    scoring_matrix_name: str="grantham_similarity_norm",
):
    substitution_matrix_df = matrices.load_precomputed_matrix_df(scoring_matrix_name)
    kmers = tools.gen_kmers(ref_idr, k)
    kmers = list(set(kmers))
    kmer_count = 0
    for kmer in kmers:
        a = time.time()
        _,_,_ = kmer_align_aligner_allseqs(ref_idr, ortholog_idrs)
        kmer_count += 1
        b = time.time()
        print(f"Finished kmer: {kmer}\n{kmer_count} of {len(kmers)}")
        print(f"Time: {b-a} seconds")
    return score_df, subseq_df, rbm_df

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








# %%


import numpy as np


def score_kmer_2_seq_no_gaps(
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

    score_list = []
    # Calculate scores for all possible subsequences
    for i in range(kmer_len):
        seqarr_i = seqarr[i:i + kmer_len]
        scores = np.array([float(substitution_matrix_df.loc[k, s]) for k, s in zip(kmerarr, seqarr_i)])
        score_list.append(np.sum(scores))
    return score_lista

substitution_matrix_df = matrices.load_precomputed_matrix_df('grantham_similarity_norm')
%timeit score_kmer_2_seq_no_gaps('KQK', 'KKKKKKKK', substitution_matrix_df)

# %%

kmer = 'KQK'
# sequence = 'KKKKKKKK'
sequence = 'ACDEFGHI'
substitution_matrix_df = matrices.load_precomputed_matrix_df('grantham_similarity_norm')


kmer_len = len(kmer)
seq_len = len(sequence)

# Convert kmer and sequence to NumPy arrays
kmerarr = np.array(list(kmer))
seqarr = np.array(list(sequence))

# Create a 2D array of substitution scores
for i in range(seq_len - kmer_len + 1):
    print(i, sequence[i:i + kmer_len])
    substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]]


sub_slices = [substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]] for i in range(seq_len - kmer_len + 1)]
np.sum([np.diag(substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]].values) for i in range(seq_len - kmer_len + 1)], axis=1)
scores = [np.trace(substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]].values) for i in range(seq_len - kmer_len + 1)]

# np.trace(np.eye(3))

# a = np.arange(8).reshape((2,2,2))


scores_matrix = np.array([substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]].values for i in range(seq_len - kmer_len + 1)])
scores_matrix = np.array([substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]] for i in range(seq_len - kmer_len + 1)])

scores_matrix = np.array([substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]].values for i in range(seq_len - kmer_len + 1)])
matrix = np.array([substitution_matrix_df.loc[kmerarr, seqarr[i:i + kmer_len]].values for i in range(seq_len - kmer_len + 1)])
matrix.shape
np.sum(matrix, axis=1).tolist()


# Sum the scores along the K-mer axis
scores = np.sum(scores_matrix, axis=1)

score_kmer_2_seq_no_gaps(kmer, sequence, substitution_matrix_df).shape


# %%
def score_kmer_2_seq_no_gaps(
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
    return score

%timeit score_kmer_2_seq_no_gaps('KQK', 'KKKKKKKK', substitution_matrix_df)
# 'KKKKKKKK'
#      'KQK'   
#     'KQK'   
#    'KQK'   
#   'KQK'    
#  'KQK'     
# 'KQK'     




# %%
# with open("matrices_with_info.json", "r") as json_file:
#     data = json.load(json_file)

# df1_reconstructed = pd.DataFrame(
#     data["matrix1"]["data"],
#     columns=data["matrix1"]["columns"],
#     index=data["matrix1"]["index"],
# )
# df2_reconstructed = pd.DataFrame(
#     data["matrix2"]["data"],
#     columns=data["matrix2"]["columns"],
#     index=data["matrix2"]["index"],
# )
# additional_info = data["additional_info"]



















    















