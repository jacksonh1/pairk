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


def get_aligner(
    matrix_file_for_aligner: str|Path = matrices.MATRIX_DIR / "grantham_similarity_normx100_aligner_compatible", 
) -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -100000
    aligner.extend_gap_score = -100000
    aligner.query_end_gap_score = 0.0
    aligner.mode = "global"
    aligner.substitution_matrix = Align.substitution_matrices.read(
        matrix_file_for_aligner
    )
    # print('getting aligner')
    return aligner


def kmer_align_aligner(
    kmer: str,
    sequence: str,
    aligner: Align.PairwiseAligner, 
) -> Align.Alignment:
    alignment = aligner.align(Seq.Seq(sequence), Seq.Seq(kmer))[0]
    return alignment


def get_aas_aligned_2_query_nongap_positions(alignment: Align.Alignment):
    query_arr = np.array(list(alignment[1])) # type: ignore
    target_arr = np.array(list(alignment[0])) # type: ignore
    nongap_indices = [c for c, i in enumerate(query_arr) if i != "-"]
    result = target_arr[nongap_indices]
    return ''.join(result), nongap_indices[0]


def get_position_of_query_start_in_target(alignment: Align.Alignment):
    query_arr = np.array(list(alignment[1])) # type: ignore
    target_arr = np.array(list(alignment[0])) # type: ignore
    for c, i in enumerate(query_arr):
        if i != "-":
            return c
    # nongap_indices = [c for c, i in enumerate(query_arr) if i != "-"]
    # return nongap_indices[0]


def first_non_dash_index(s):
    index = 0
    while index < len(s) and s[index] == '-':
        index += 1
    return index


def score_kmer_2_seq_no_gaps_v4(
    kmer: str,
    sequence: str,
    aligner: Align.PairwiseAligner,
):
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")
    aln = kmer_align_aligner(kmer, sequence, aligner=aligner)
    bseq = get_aas_aligned_2_query_nongap_positions(aln)
    bpos = get_position_of_query_start_in_target(aln)
    print(bseq, bpos)
    # s = first_non_dash_index(aln[1])
    # bseq = aln[0][s:s+len(kmer)]
    # bpos = s
    return aln.score, bseq, bpos


def make_empty_kmer_ortho_df(kmers: list[str], ortholog_ids: list[str]):
    return pd.DataFrame(
        index=kmers,
        columns=ortholog_ids,
    )


def run_pairwise_kmer_alignment(
    ref_idr: str,
    ortholog_idrs: dict[str, str],
    k: int,
    matrix_file_for_aligner: str|Path=matrices.ALIGNER_MATRIX_FILE_DICT["grantham_similarity"],
):
    # substitution_matrix_df = matrices.load_precomputed_matrix_df(scoring_matrix_name)
    kmers = tools.gen_kmers(ref_idr, k)
    kmers = list(set(kmers))
    score_df = make_empty_kmer_ortho_df(kmers, list(ortholog_idrs.keys()))
    subseq_df = make_empty_kmer_ortho_df(kmers, list(ortholog_idrs.keys()))
    pos_df = make_empty_kmer_ortho_df(kmers, list(ortholog_idrs.keys()))
    rbm_df = make_empty_kmer_ortho_df(kmers, list(ortholog_idrs.keys()))
    aligner = get_aligner(matrix_file_for_aligner=matrix_file_for_aligner)
    for kmer in kmers:
        for ortholog_id, ortholog_idr in ortholog_idrs.items():
            if len(ortholog_idr) < k:
                continue
            try:
                best_score, best_subseq, best_pos = score_kmer_2_seq_no_gaps_v4(kmer, ortholog_idr, aligner)
            except ValueError as e:
                print(e)
                continue
            score_df.loc[kmer, ortholog_id] = best_score
            subseq_df.loc[kmer, ortholog_id] = best_subseq
            pos_df.loc[kmer, ortholog_id] = best_pos
            try:
                reci_best_score, reci_best_subseq, _= score_kmer_2_seq_no_gaps_v4(
                    best_subseq,
                    ref_idr,
                    aligner,
                )
            except ValueError as e:
                print(e)
                continue
            rbm_df.loc[kmer, ortholog_id] = (reci_best_subseq==kmer)
    return score_df, subseq_df, pos_df, rbm_df


def main(
    input_alignment_file: str|Path,
    output_file: str|Path,
    reference_id: str,
    k: int,
    idr_aln_st: int,
    idr_aln_end: int,
    overwrite: bool = False,
    aligner_matrix_name: str="grantham_similarity",
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
    if aligner_matrix_name not in matrices.ALIGNER_MATRIX_FILE_DICT:
        raise ValueError(f"aligner_matrix_name {aligner_matrix_name} not valid. Must be one of {list(matrices.ALIGNER_MATRIX_FILE_DICT.keys())}")
    matrix_file_for_aligner = matrices.ALIGNER_MATRIX_FILE_DICT[aligner_matrix_name]
    score_df, subseq_df, pos_df, rbm_df = run_pairwise_kmer_alignment(
        ref_idr,
        idrs, # type: ignore
        k,
        matrix_file_for_aligner=matrix_file_for_aligner
    )
    output_dict = {
        "score_dataframe": score_df.to_dict(orient="split"),
        "subseq_dataframe": subseq_df.to_dict(orient="split"),
        "position_dataframe": pos_df.to_dict(orient="split"),
        "reciprocal_best_match_dataframe": rbm_df.to_dict(orient="split"),
    }
    with open(output_file, "w") as json_file:
        json.dump(output_dict, json_file, indent=4)


# ==============================================================================
# // testing
# ==============================================================================    
import local_conservation_analysis_pipeline.group_conservation_objects as group_tools

json_file = "../../../benchmark/benchmark_multi_species/p3_conservation_analysis_pipeline/conservation_analysis/10-9606_0_00051d/10-9606_0_00051d.json"
og = group_tools.ConserGene(json_file)
lvlo=og.get_level_obj('Vertebrata')

def pad_hit(seq: str, st_pos: int, end_pos: int, l_flank: int = 0, r_flank: int = 0):
    st = max(0, st_pos - l_flank)
    end = min(len(seq)-1, end_pos + r_flank)
    return st, end, seq[st : end + 1]

_,_,flanked_hit = pad_hit(og.query_sequence, og.hit_start_position, og.hit_end_position, 5, 5)
k = len(flanked_hit)
a = time.time()
main(
    lvlo.alignment_file,
    output_file="./test.json",
    reference_id=og.query_gene_id,
    k=k,
    idr_aln_st=lvlo.idr_aln_start,
    idr_aln_end=lvlo.idr_aln_end,
    overwrite=True,
)
b = time.time()
print(f"Time: {b-a} seconds")