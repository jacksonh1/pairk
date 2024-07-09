from pairk.backend.tools.sequence_utils import aln_2_idr_position_map
import pairk.backend.tools.sequence_utils as _tools  # hiding this module from the user
from pathlib import Path


def fasta_MSA_to_idr_dict(
    alignment_file: str | Path, idr_aln_start: int, idr_aln_end: int
) -> dict[str, str]:
    """import a multiple sequence alignment (MSA) in fasta format and return the IDR sequences as a dictionary. uses a single slice of the MSA to define the IDRs of each sequence.

    Parameters
    ----------
    alignment_file : str | Path
        a path to a fasta file containing a multiple sequence alignment
    idr_aln_start : int
        the start position of the IDR in the MSA
    idr_aln_end : int
        the end position of the IDR in the MSA

    Returns
    -------
    dict[str, str]
        a dictionary with sequence IDs as keys and the IDR region as values
    """
    faimporter = _tools.FastaImporter(alignment_file)
    aln = faimporter.import_as_alignment()
    idr_aln_list = list(aln[:, idr_aln_start : idr_aln_end + 1])
    idr_dict = {
        i.id: str(i.seq) for i in _tools.strip_dashes_from_sequences(idr_aln_list)  # type: ignore
    }
    return idr_dict  # type: ignore


def fasta_MSA_to_idr_map(
    alignment_file: str | Path, idr_aln_start: int, idr_aln_end: int
):
    """import a multiple sequence alignment (MSA) from a fasta file and return a
    dictionary of the start and end positions of the IDRs in the unaligned sequences.
    The IDRs are defined by a single slice of the alignment.

    Parameters
    ----------
    aln : str | Path
        a path to a fasta file containing a multiple sequence alignment
    idr_aln_start : int
        start position of the IDRs in the alignment
    idr_aln_end : int
        end position of the IDRs in the alignment

    Returns
    -------
    dict[str, list[int]]
        idr_position_map: dictionary of the start and end positions of the IDRs
        in the unaligned sequences. The keys are the sequence ids and the values
        are lists of the start and end positions of the IDRs in the unaligned
        sequences. The start and end positions are 0-indexed.
    """
    faimporter = _tools.FastaImporter(alignment_file)
    aln = faimporter.import_as_alignment()
    idr_map = aln_2_idr_position_map(aln, idr_aln_start, idr_aln_end + 1)
    return idr_map


def fasta_MSA_to_unaligned_sequences(alignment_file: str | Path) -> dict[str, str]:
    """import a multiple sequence alignment (MSA) from a fasta file and return the
    unaligned sequences as a dictionary.

    Parameters
    ----------
    alignment_file : str | Path
        a path to a fasta file containing a multiple sequence alignment

    Returns
    -------
    dict[str, str]
        a dictionary with sequence IDs as keys and the unaligned sequence as values
    """
    faimporter = _tools.FastaImporter(alignment_file)
    aln_dict = faimporter.import_as_str_dict()
    unaligned_dict = {k: _tools.strip_dashes_from_str(v) for k, v in aln_dict.items()}
    return unaligned_dict


# generate IDR map and IDR dict from
