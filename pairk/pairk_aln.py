"""primary user-facing functions/objects for the pairk package"""

"""
exhaustive_pairk_align:
    - reference_id: str    
    - idr_dictionary: dict[str, str]
    - k: int
    - matrix: dict[str, dict[str, float]]


needleman_pairk_align:
    - reference_id: str    
    - idr_dictionary: dict[str, str]
    - k: int
    - aligner: Bio.Align.PairwiseAligner


input set 1:
    - reference_id: str
    - MSA: list[str]
    - idr_aln_st: int
    - idr_aln_end: int
    - k: int
    - alignment_type: str = ("exhaustive", "dynamic_programming")
    - matrix_name: str
        - or-
    - matrix_file: str = "grantham_similarity_normx100_aligner_compatible"
"""

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
import pairk.tools.sequence_utils as tools
import pairk.kmer_alignment.needleman_tools as needleman_tools
from pathlib import Path
from pairk.kmer_alignment.scoring_matrix_needleman import pairk_alignment_needleman
from pairk.kmer_alignment.scoring_matrix import pairk_alignment


def pairwise_kmer_alignment():
    pass
