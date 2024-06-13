import pandas as pd
from alfpy import word_distance, word_pattern, word_vector
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords as alf_seqrecords
from Bio import Align, AlignIO, Seq, SeqIO
from Bio.SeqRecord import SeqRecord


def score_alignment(
    seq1: str, seq2: str, subs_mat_df: pd.DataFrame, gap_open=-10, gap_extend=-0.5
) -> float:
    """returns the score of an alignment between two sequences

    treats the following situation as opening 2 gaps:
        AAAA---
        ----BBB

    you can use this function to compute a similarity score between two sequences
    if you normalize the score by the max possible score which might be something like:
    >>> max_score = max(score_alignment(seq1, seq1), score_alignment(seq2, seq2))

    so the similarity score would be:
    >>> max_score = max(score_alignment(seq1, seq1), score_alignment(seq2, seq2))
    >>> similarity_score = score_alignment(seq1, seq2) / max_score

    Parameters
    ----------
    seq1 : str
        first sequence
    seq2 : str
        second sequence
    subs_mat_df : pd.DataFrame
        the scoring matrix as a pandas dataframe
    gap_open : int, optional
        penalty for opening a gap (should be negative), by default -10
    gap_extend : int, optional
        penalty for extending a gap (should be negative), by default -0.5

    Returns
    -------
    float
        the score of the alignment
    """
    assert len(seq1) == len(seq2)
    score = 0
    gap1_open_flag = False
    gap2_open_flag = False
    for s1, s2 in zip(seq1, seq2):
        if s1 == "-":
            if gap1_open_flag:
                score += gap_extend
            else:
                score += gap_open
                gap1_open_flag = True
        elif s2 == "-":
            if gap2_open_flag:
                score += gap_extend
            else:
                score += gap_open
                gap2_open_flag = True
        else:
            score += float(subs_mat_df.loc[s1, s2])  # type: ignore
            gap1_open_flag = False
            gap2_open_flag = False
    return score


def score_alignment_from_alignment_obj(
    alignment_obj, subs_mat_df, gap_open, gap_extend
):
    seq1 = alignment_obj[0]
    seq2 = alignment_obj[1]
    return score_alignment(seq1, seq2, subs_mat_df, gap_open, gap_extend)


def pairwise_alignment(
    seq1,
    seq2,
    scoring_matrix_name="BLOSUM62",
    gap_opening_penalty=10,
    gap_extension_penalty=0.5,
):
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = Align.substitution_matrices.load(scoring_matrix_name)
    aligner.extend_gap_score = -gap_extension_penalty
    aligner.open_gap_score = -gap_opening_penalty
    alignment = aligner.align(seq1, seq2)[0]
    return alignment


def compute_pairwise_percent_id_from_msa(
    seqrecord1: SeqRecord, seqrecord2: SeqRecord
) -> float:
    """compute pairwise percent identity between two sequences
    - The sequences must be pre-aligned (i.e. they have the same length)
    - The returned percent identity is computed as the number of identical
    residues divided by the LENGTH OF compared alignment sites. Alignment sites
    where both sequences are gaps ('-' characters) are not compared, and thus
    do not contribute to the length.

    Parameters
    ----------
    seqrecord1 : seqrecord
        first sequence
    seqrecord2 : seqrecord
        second sequence

    Returns
    -------
    float
        percent identity between two sequences.
    """
    seq1 = str(seqrecord1.seq)
    seq2 = str(seqrecord2.seq)
    # print(seq1)
    # print(seq2)
    assert len(seq1) == len(seq2), "sequences are not the same length"
    num_same = 0
    length = len(seq1)
    for i in range(len(seq1)):
        if seq1[i] == "-" and seq2[i] == "-":
            length -= 1
            continue
        if seq1[i] == seq2[i]:
            num_same += 1
    return num_same / length


def percent_identity(seq1: str | SeqRecord, seq2: str | SeqRecord) -> float:
    """
    Returns a percent identity between 0 (no identical residues) and 1 (all residues are identical)
    - The sequences must be pre-aligned (i.e. they have the same length)
    - The returned percent identity is computed as the number of identical residues divided by the
    length of compared alignment sites. Alignment sites where both sequences are gaps ('-' characters)
    are not compared.

    Parameters
    ----------
    seq1 : str or seqrecord
        first sequence
    seq2 : str or seqrecord
        second sequence
    """
    assert len(seq1) == len(
        seq2
    ), "sequences are not the same length. Are they aligned?"
    num_same = 0
    length = len(seq1)
    for i in range(len(seq1)):
        if seq1[i] == "-" and seq2[i] == "-":
            length -= 1
            continue
        if seq1[i] == seq2[i]:
            num_same += 1
    pid = num_same / length
    return pid


def align_and_get_PID(seqrecord1, seqrecord2, **kwargs):
    """
    aligns the two sequences with Biopython's PairwiseAligner using the BLOSUM62 scoring matrix and gap_opening_penalty = 10 and gap_extension_penalty = 0.5 by default.
    The alignment parameters can be changed by passing them as keyword arguments:
    `scoring_matix_name`, `gap_opening_penalty`, and `gap_extension_penalty`
    They are passed to the `pairwise_alignment` function
    """
    aln = pairwise_alignment(seqrecord1.seq, seqrecord2.seq, **kwargs)
    pid = percent_identity(aln[0], aln[1])  # type: ignore
    return pid


def alfpy_distance_matrix(seqrecord_list, word_size=2):
    id_list = [i.id for i in seqrecord_list]
    seq_str_list = [str(i.seq) for i in seqrecord_list]
    alf_seq_records = alf_seqrecords.SeqRecords(id_list=id_list, seq_list=seq_str_list)
    p = word_pattern.create(alf_seq_records.seq_list, word_size=word_size)
    # counts = word_vector.Counts(alf_seq_records.length_list, p)
    freqs = word_vector.Freqs(alf_seq_records.length_list, p)
    dist = word_distance.Distance(freqs, "google")
    matrix = distmatrix.create(alf_seq_records.id_list, dist)
    return matrix
