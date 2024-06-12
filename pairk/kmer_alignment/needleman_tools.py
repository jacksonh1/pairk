from Bio import Align, Seq
import pairk.tools.sequence_utils as tools


def get_aligner(
    matrix: Align.substitution_matrices.Array,  # type: ignore
) -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -100000
    aligner.extend_gap_score = -100000
    aligner.query_end_gap_score = 0.0
    aligner.mode = "global"
    aligner.substitution_matrix = matrix
    return aligner


def kmer_align_aligner(
    kmer: str,
    sequence: str,
    aligner: Align.PairwiseAligner,
) -> Align.Alignment:
    alignment = aligner.align(Seq.Seq(sequence), Seq.Seq(kmer))[0]
    return alignment


def score_kmer_2_seq_no_gaps_needleman(
    kmer: str,
    sequence: str,
    aligner: Align.PairwiseAligner,
):
    if len(sequence) < len(kmer):
        raise ValueError("The sequence is shorter than the kmer")
    aln = kmer_align_aligner(kmer, sequence, aligner=aligner)
    s = tools.get_first_non_gap_index(aln[1])  # type: ignore
    best_subseq = aln[0][s : s + len(kmer)]
    best_pos = s
    return aln.score, best_subseq, best_pos  # type: ignore
