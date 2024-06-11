#!/usr/bin/env python


import argparse
# %%
import json
from argparse import RawTextHelpFormatter
from pathlib import Path

import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO, SeqRecord

import local_env_variables.matrices as submats
from local_conservation_scores.tools import general as cons_tools


def mask_alignment(
    alignment: Align.MultipleSeqAlignment,
    reference_seq_str: str,
    gap_frac_cutoff: float = 0.2,
):
    gap_mask, score_mask = cons_tools.make_score_mask(
        alignment,
        reference_seq_str,
        gap_frac_cutoff=gap_frac_cutoff,
    )
    gap_mask = [bool(i) for i in gap_mask]
    score_mask = [bool(i) for i in score_mask]
    return gap_mask, score_mask


def score_alignment(
    reference_seqrecord: SeqRecord.SeqRecord,
    alignment: Align.MultipleSeqAlignment,
    matrix_df: pd.DataFrame,
):
    scores = cons_tools.asymmetric_valdar_score_df_unweighted(
        reference_seqrecord,
        alignment,
        matrix_df,
    )
    return scores


def main(
    input_alignment_file,
    output_file,
    reference_id,
    matrix_name,
    gap_frac_cutoff=0.2,
    overwrite=False,
    **kwargs,
):
    """
    calculate the asymmetric sum-of-pairs score for each column in an alignment

    Parameters
    ----------
    input_alignment_file : str|Path
        input alignment file (fasta format)
    output_file : str|Path
        output file (json). The output file will contain the following:
        "gap_mask" : list[bool]
            a list of bools indicating whether a column is masked by the gap mask
        "score_mask" : list[bool]
            a list of bools indicating whether a column is masked by the score mask. The score mask is the gap mask with additional columns masked if they are gaps in the reference sequence
        "scores" : list[float]
            a list of scores for each column in the alignment
    reference_id : str
        the id of the sequence to use as the reference for the score mask (gaps in this sequence will be masked).
    matrix_name : str
        the name of the matrix to use for scoring.
        Available matrices:
            BLOSUM62_row_norm
            BLOSUM62_max_off_diagonal_norm
            EDSSMat50_row_norm
            EDSSMat50_max_off_diagonal_norm
    gap_frac_cutoff : float
        A number between 0 and 1. he fraction of gaps allowed in a column. If column has >gap_frac_cutoff gaps, it is masked in the gap mask. by default 0.2
    overwrite : bool, optional
        if True, overwrites the `output_file` if it already exists, by default False
    """
    with open(input_alignment_file, "r") as f:
        alignment = AlignIO.read(f, "fasta")

    # check if output file exists
    if Path(output_file).exists() and not overwrite:
        print(f"{output_file} exists and overwrite is False")
        print("exiting...")
        return

    matrix_df = submats.load_precomputed_matrix_df(matrix_name)
    alignment_seqrecord_dict = {seqrecord.id: seqrecord for seqrecord in alignment}

    # if reference_id is None:
    #     reference_id = alignment[0].id
    reference_seqrecord = alignment_seqrecord_dict[reference_id]
    reference_seq_str = str(reference_seqrecord.seq)
    score_dict = {}
    score_dict["gap_mask"], score_dict["score_mask"] = mask_alignment(
        alignment,
        reference_seq_str,
        gap_frac_cutoff=gap_frac_cutoff,
    )
    score_dict["scores"] = score_alignment(reference_seqrecord, alignment, matrix_df)

    with open(output_file, "w") as f:
        json.dump(score_dict, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="calculate asymmetric sum of pairs conservation scores (valdar) for an input alignment file",
        formatter_class=RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        metavar="<file>",
        required=True,
        help="input alignment file (fasta format)",
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        metavar="<file>",
        required=True,
        help="output file (json)",
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="overwrite existing results"
    )
    # optional reference id to calculate conservation scores relative to. If not provided, the first sequence in the alignment will be used
    parser.add_argument(
        "-r",
        "--reference_id",
        type=str,
        metavar="<str>",
        default=None,
        help="""reference id to calculate conservation scores relative to. If not provided, the first sequence in the alignment will be used""",
    )
    parser.add_argument(
        "-m",
        "--matrix",
        type=str,
        metavar="<str>",
        default="EDSSMat50_max_off_diagonal_norm",
        help="""matrix name (e.g. blosum62)
    available matrices:
        BLOSUM62
        BLOSUM62_norm
        BLOSUM62_row_norm
        BLOSUM62_max_off_diagonal
        BLOSUM62_max_off_diagonal_norm
        EDSSMat50
        EDSSMat50_norm
        EDSSMat50_row_norm
        EDSSMat50_max_off_diagonal
        EDSSMat50_max_off_diagonal_norm
        grantham_similarity_norm
if you are going to be calculating a z-score, do not use a matrix that has different values in the diagonal. All of the diagonals must be the same value or else the z-score will be completely biased.""",
    )
    parser.add_argument(
        "-g",
        "--gap_frac_cutoff",
        type=float,
        metavar="<float>",
        default=0.2,
        help="""fraction of gaps allowed in a column. If column has >gap_frac_cutoff, it is masked in the gap mask, default=0.2""",
    )
    args = parser.parse_args()
    main(
        args.input,
        args.output_file,
        reference_id=args.reference_id,
        matrix_name = args.matrix,
        gap_frac_cutoff = args.gap_frac_cutoff,
        overwrite=args.overwrite,
    )
