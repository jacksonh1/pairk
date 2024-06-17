import pairk
from importlib.resources import files as _files
from pathlib import Path
import pairk.backend.tools.sequence_utils as tools


class ExampleAlignment:
    """
    Class to store an example alignment and the IDRs within it.

    Attributes:
    ----------
    alignment_file : pathlib.Path
        Path to the alignment file
    query_id : str
        ID of the query sequence in the alignment file
    aln_idr_start : int
        Start position of the IDR in the alignment
    aln_idr_end : int
        End position of the IDR in the alignment
    idr_start : int
        Start position of the IDR in the unaligned sequence
    idr_end : int
        End position of the IDR in the unaligned sequence
    alignment : Bio.Align.MultipleSeqAlignment
        Alignment object
    idr_dict : dict
        Dictionary of the IDRs (unaligned)
    """

    def __init__(
        self, alignment_file, query_id, aln_idr_start, aln_idr_end, idr_start, idr_end
    ):
        self.alignment_file = alignment_file
        self.query_id = query_id
        self.aln_idr_start = aln_idr_start
        self.aln_idr_end = aln_idr_end
        self.idr_start = idr_start
        self.idr_end = idr_end
        self.alignment = self._import_alignment()
        idr_aln_list = list(
            self.alignment[:, self.aln_idr_start : self.aln_idr_end + 1]
        )
        self.idr_dict = {
            i.id: str(i.seq) for i in tools.strip_dashes_from_sequences(idr_aln_list)  # type: ignore
        }
        self.full_length_dict = {
            i.id: str(i.seq)
            for i in tools.strip_dashes_from_sequences(list(self.alignment))
        }
        self.idr_position_map = self._find_unaligned_idr_position_map()

    def __str__(self):
        return (
            f"Alignment file: {self.alignment_file}\n"
            f"Query ID: {self.query_id}\n"
            f"IDR start in MSA: {self.aln_idr_start}\n"
            f"IDR end in MSA: {self.aln_idr_end}\n"
            f"IDR start: {self.idr_start}\n"
            f"IDR end: {self.idr_end}\n"
        )

    def __repr__(self):
        return f"ExampleAlignment({self.alignment_file}, {self.aln_idr_start}, {self.aln_idr_end}, {self.idr_start}, {self.idr_end})"

    def _import_alignment(self):
        faimporter = pairk.FastaImporter(self.alignment_file)
        return faimporter.import_as_alignment()

    def _find_unaligned_idr_position_map(self):
        idr_position_map = {}
        for seqrecord in list(self.alignment):
            p = tools.find_alnslice_positions_in_unaln(
                str(seqrecord.seq), self.aln_idr_start, self.aln_idr_end + 1
            )
            idr_position_map[seqrecord.id] = [p[0], p[-1]]
        return idr_position_map


example1_alignment_file = Path(
    _files("pairk.data")
    .joinpath("example_alignment_9606_0_00294e-idraln-555_to_971-idr-440_to_665.fasta")
    .__fspath__()  # type: ignore
)

example1 = ExampleAlignment(
    alignment_file=example1_alignment_file,
    query_id="9606_0:00294e",
    aln_idr_start=555,
    aln_idr_end=971,
    idr_start=440,
    idr_end=665,
)
