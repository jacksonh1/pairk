import io
import json

from Bio import AlignIO

import local_seqtools.general_utils as tools


class jch_alignment:
    """This class stores an alignment as a AlignIO.MultipleSeqAlignment object and columnwise annotations.
    It stores an index for converting between the unaligned and aligned sequence positions. The "query" or
    "reference sequence must be the first sequence in the alignment (for now).
    The main goal is to allow for easy slicing of the alignment and the column annotations.
    """

    def __init__(
        self,
        alignment: AlignIO.MultipleSeqAlignment,
        query_id: str,
        column_annotations: dict[str, list] | None = None,
    ):
        # Store the alignment
        self.alignment = alignment
        # query sequence is the first sequence in the alignment
        self.query_seqrecord = [
            self.alignment[i]
            for i in range(len(self.alignment))
            if self.alignment[i].id == query_id
        ][0]
        self.query_sequence_str = str(self.query_seqrecord.seq)
        self.query_sequence_id_str = self.query_seqrecord.id
        # Store the index and unaligned sequences
        self.query_unaligned_str, self.index = tools.reindex_alignment_str(
            self.query_sequence_str
        )
        if column_annotations is not None:
            self.alignment.column_annotations = column_annotations

    def __str__(self):
        s = f"jch_alignment object\n{self.alignment.get_alignment_length()} columns\n{len(self.alignment)} sequences\ncolumn annotations: {self.alignment.column_annotations.keys()}"
        return s

    def unaligned_annotation_to_column_annotation(self, unaligned_annotation: list):
        """Takes a list of annotations that correspond to the unaligned sequence and returns a list of annotations that correspond to the aligned sequence"""
        column_annotation = []
        unaligned_pos = 0
        for char in self.query_sequence_str:
            if char == "-":
                column_annotation.append("-")
            else:
                column_annotation.append(unaligned_annotation[unaligned_pos])
                unaligned_pos += 1
        return column_annotation

    def add_unaligned_annotation_to_column_annotations(
        self, key: str, unaligned_annotation: list
    ):
        """Takes a list of annotations that correspond to the unaligned sequence and adds it to the column annotations for the alignment"""
        column_annotation = self.unaligned_annotation_to_column_annotation(
            unaligned_annotation
        )
        self.alignment.column_annotations[key] = column_annotation

    def slice_by_unaligned_positions_inclusive(self, start, end):
        """INCLUSIVE. Returns a new alignment object that is a slice of the current
        alignment.  The slice is defined by the start and end positions of the
        query sequence."""
        # Get the start and end positions in the alignment
        start_aln, end_aln = self.index[start], self.index[end]
        # Slice the alignment
        aln_slice = self.alignment[:, start_aln : end_aln + 1]
        return jch_alignment(aln_slice)

    def slice_by_aln_positions(self, start, end):
        """Returns a new alignment object that is a slice of the current
        alignment.  The slice is defined by the start and end positions of the
        alignment."""
        # Slice the alignment
        aln_slice = self.alignment[:, start:end]
        return jch_alignment(aln_slice)

    def format_export_dict(self):
        export_dict = {}
        export_dict["alignment"] = format(self.alignment, "fasta")
        export_dict["column_annotations"] = self.alignment.column_annotations
        return export_dict

    def write_json(self, filename):
        export_dict = self.format_export_dict()
        with open(filename, "w") as f:
            json.dump(export_dict, f, indent=4)

    @classmethod
    def from_seq_file(cls, filename, format="fasta"):
        """Reads in an alignment from a file.  The file can be in any format
        that is supported by biopython."""
        aln = AlignIO.read(filename, format)
        return cls(aln)

    @classmethod
    def from_seqrecord_list(cls, seqrecord_list):
        """Creates an alignment from a list of SeqRecord objects."""
        aln = AlignIO.MultipleSeqAlignment(seqrecord_list)
        return cls(aln)

    @classmethod
    def from_json(cls, filename):
        with open(filename, "r") as f:
            data = json.load(f)
        text_fa = data["alignment"]
        aln = AlignIO.read(io.StringIO(text_fa), "fasta")
        column_annotations = data["column_annotations"]
        return cls(aln, column_annotations)

    @classmethod
    def from_dict(cls, aln_dict: dict):
        text_fa = aln_dict["alignment"]
        aln = AlignIO.read(io.StringIO(text_fa), "fasta")
        column_annotations = aln_dict["column_annotations"]
        return cls(aln, column_annotations)


# class _RestrictedDict(dict):
#     """Dict which only allows sequences of given length as values (PRIVATE).

#     This simple subclass of the Python dictionary is used in the SeqRecord
#     object for holding per-letter-annotations.  This class is intended to
#     prevent simple errors by only allowing python sequences (e.g. lists,
#     strings and tuples) to be stored, and only if their length matches that
#     expected (the length of the SeqRecord's seq object).  It cannot however
#     prevent the entries being edited in situ (for example appending entries
#     to a list).

#     >>> x = _RestrictedDict(5)
#     >>> x["test"] = "hello"
#     >>> x
#     {'test': 'hello'}

#     Adding entries which don't have the expected length are blocked:

#     >>> x["test"] = "hello world"
#     Traceback (most recent call last):
#     ...
#     TypeError: We only allow python sequences (lists, tuples or strings) of length 5.

#     The expected length is stored as a private attribute,

#     >>> x._length
#     5

#     In order that the SeqRecord (and other objects using this class) can be
#     pickled, for example for use in the multiprocessing library, we need to
#     be able to pickle the restricted dictionary objects.

#     Using the default protocol, which is 3 on Python 3,

#     >>> import pickle
#     >>> y = pickle.loads(pickle.dumps(x))
#     >>> y
#     {'test': 'hello'}
#     >>> y._length
#     5

#     Using the highest protocol, which is 4 on Python 3,

#     >>> import pickle
#     >>> z = pickle.loads(pickle.dumps(x, pickle.HIGHEST_PROTOCOL))
#     >>> z
#     {'test': 'hello'}
#     >>> z._length
#     5
#     """

#     def __init__(self, length):
#         """Create an EMPTY restricted dictionary."""
#         dict.__init__(self)
#         self._length = int(length)

#     def __setitem__(self, key, value):
#         # The check hasattr(self, "_length") is to cope with pickle protocol 2
#         # I couldn't seem to avoid this with __getstate__ and __setstate__
#         if (
#             not hasattr(value, "__len__")
#             or not hasattr(value, "__getitem__")
#             or (hasattr(self, "_length") and len(value) != self._length)
#         ):
#             raise TypeError(
#                 "We only allow python sequences (lists, tuples or strings) "
#                 f"of length {self._length}."
#             )
#         dict.__setitem__(self, key, value)

#     def update(self, new_dict):
#         # Force this to go via our strict __setitem__ method
#         for (key, value) in new_dict.items():
#             self[key] = value


# class MultipleSeqAlignment:

#     def __init__(
#         self, records, column_annotations=None
#     ):
#         self.records = records
#         # Annotations about each column of the alignment
#         if column_annotations is None:
#             column_annotations = {}
#         # Handle this via the property set function which will validate it
#         self.column_annotations = column_annotations

#     def _set_per_column_annotations(self, value):
#         if not isinstance(value, dict):
#             raise TypeError(
#                 "The per-column-annotations should be a (restricted) dictionary."
#             )
#         # Turn this into a restricted-dictionary (and check the entries)
#         if len(self):
#             # Use the standard method to get the length
#             expected_length = self.get_alignment_length()
#             self._per_col_annotations = _RestrictedDict(length=expected_length)
#             self._per_col_annotations.update(value)


#     def _get_per_column_annotations(self):
#         if self._per_col_annotations is None:
#             # This happens if empty at initialisation
#             if len(self):
#                 # Use the standard method to get the length
#                 expected_length = self.get_alignment_length()
#             else:
#                 # Should this raise an exception? Compare SeqRecord behaviour...
#                 expected_length = 0
#             self._per_col_annotations = _RestrictedDict(length=expected_length)
#         return self._per_col_annotations

#     column_annotations = property(
#         fget=_get_per_column_annotations,
#         fset=_set_per_column_annotations,
#         doc="""Dictionary of per-letter-annotation for the sequence.""",
#     )

#     def __len__(self):
#         """Return the number of sequences in the alignment.

#         Use len(alignment) to get the number of sequences (i.e. the number of
#         rows), and alignment.get_alignment_length() to get the length of the
#         longest sequence (i.e. the number of columns).

#         This is easy to remember if you think of the alignment as being like a
#         list of SeqRecord objects.
#         """
#         return len(self._records)


# class ColumnAnnotations(dict):
#     def __init__(self, alignment: AlignIO.MultipleSeqAlignment, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.alignment = alignment

#     def __setitem__(self, key, value):
#         if len(value) != len(self.alignment[0]):
#             raise ValueError(f"problem with entry `{key}`, The number of values in column_annotations must match the alignment length.")
#         super().__setitem__(key, value)


# class jch_alignment:
#     '''The main function of this class should be to store the alignment and
#     index between the alignment and the unaligned sequences.'''

#     def __init__(self, alignment: AlignIO.MultipleSeqAlignment, column_annotations: dict[str, list]|None=None):
#         # Store the alignment
#         self.alignment = alignment
#         # query sequence is the first sequence in the alignment
#         self.query_seqrecord = self.alignment[0]
#         # Store the index and unaligned sequences
#         self.query_unaligned_str, self.index = tools.reindex_alignment(self.query_seqrecord)
#         self.column_annotations = ColumnAnnotations(self.alignment)
#         if column_annotations is not None:
#             for key, value in column_annotations.items():
#                 self.column_annotations[key] = value

#     # @property
#     # def column_annotations(self):
#     #     # make sure that the column annotations are the same length as the alignment
#     #     if len(self._column_annotations) != self.alignment.get_alignment_length():
#     #         raise ValueError('Column annotations are not the same length as the alignment')
#     #     return self._column_annotations


#     def slice_by_unaligned_positions_inclusive(self, start, end):
#         '''Returns a new alignment object that is a slice of the current
#         alignment.  The slice is defined by the start and end positions of the
#         query sequence.'''
#         # Get the start and end positions in the alignment
#         start_aln, end_aln = self.index[start], self.index[end]
#         # Slice the alignment
#         aln_slice = self.alignment[:, start_aln:end_aln+1]
#         # slice the column annotations
#         col_annots = {}
#         for key, value in self.column_annotations.items():
#             col_annots[key] = value[start_aln:end_aln]
#         return jch_alignment(aln_slice, col_annots)

#     def slice_by_aln_positions(self, start, end):
#         '''Returns a new alignment object that is a slice of the current
#         alignment.  The slice is defined by the start and end positions of the
#         alignment.'''
#         # Slice the alignment
#         aln_slice = self.alignment[:, start:end]
#         # slice the column annotations
#         col_annots = {}
#         for key, value in self.column_annotations.items():
#             col_annots[key] = value[start:end]
#         return jch_alignment(aln_slice, col_annots)

#     def format_export_dict(self):
#         export_dict = {}
#         export_dict['alignment'] = format(self.alignment, 'fasta')
#         export_dict['column_annotations'] = self.column_annotations
#         return export_dict

#     def write_json(self, filename):
#         export_dict = self.format_export_dict()
#         with open(filename, 'w') as f:
#             json.dump(export_dict, f, indent=4)

#     @classmethod
#     def from_seq_file(cls, filename, format='fasta'):
#         '''Reads in an alignment from a file.  The file can be in any format
#         that is supported by biopython.'''
#         aln = AlignIO.read(filename, format)
#         return cls(aln)

#     @classmethod
#     def from_seqrecord_list(cls, seqrecord_list):
#         '''Creates an alignment from a list of SeqRecord objects.'''
#         aln = AlignIO.MultipleSeqAlignment(seqrecord_list)
#         return cls(aln)

#     @classmethod
#     def from_json(cls, filename):
#         with open(filename, 'r') as f:
#             data = json.load(f)
#         text_fa = data['alignment']
#         aln = AlignIO.read(io.StringIO(text_fa), 'fasta')
#         column_annotations = data['column_annotations']
#         return cls(aln, column_annotations)

#     @classmethod
#     def from_dict(cls, aln_dict: dict):
#         text_fa = data['alignment']
#         aln = AlignIO.read(io.StringIO(text_fa), 'fasta')
#         column_annotations = data['column_annotations']
#         return cls(aln, column_annotations)
