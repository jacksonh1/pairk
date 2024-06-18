from Bio import Align, Seq
from Bio.Align import Alignment
import pairk.backend.tools.sequence_utils as tools
import pairk.backend.tools.matrices as matrices
import pairk.backend.exceptions as _exceptions


def make_aligner(
    matrix_name: str,
) -> Align.PairwiseAligner:
    """generates a Bio.Align.PairwiseAligner object with the given matrix.

    The aligner is set to global mode, with open and extend gap scores set to
    extremely low values to prevent gaps from being introduced.

    Parameters
    ----------
    matrix_name : str
        the substitution matrix to be used in the alignment. The available matrices can be viewed with the function `print_available_matrices()` in `pairk.backend.tools.matrices`. Some of them may not be compatible with the Biopython aligner object, but they will raise an error if they are not.

    Returns
    -------
    Align.PairwiseAligner
        the aligner object with the given matrix
    """
    _exceptions.validate_matrix_name(matrix_name)
    matrix = matrices.load_matrix_for_aligner(matrix_name)
    aligner = Align.PairwiseAligner()
    aligner.open_gap_score = -1000000
    aligner.extend_gap_score = -1000000
    aligner.query_end_gap_score = 0.0
    aligner.mode = "global"
    aligner.substitution_matrix = matrix
    return aligner


def kmer_align_aligner(
    kmer: str,
    sequence: str,
    aligner: Align.PairwiseAligner,
) -> Alignment:
    alignment = aligner.align(Seq.Seq(sequence), Seq.Seq(kmer))[0]
    return alignment


def score_kmer_2_seq_needleman(
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
