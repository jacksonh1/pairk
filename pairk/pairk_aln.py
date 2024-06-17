"""primary user-facing functions/objects for the alignment step of the pairk package"""

"""
exhaustive_pairk_align:
    - query_id: str    
    - idr_dictionary: dict[str, str]
    - k: int
    - matrix: dict[str, dict[str, float]]


needleman_pairk_align:
    - query_id: str    
    - idr_dictionary: dict[str, str]
    - k: int
    - aligner: Bio.Align.PairwiseAligner


input set 1:
    - query_id: str
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
from pathlib import Path
from pairk.backend.kmer_alignment.needleman_tools import make_aligner
from pairk.backend.kmer_alignment.scoring_matrix_needleman import (
    pairk_alignment_needleman,
)
from pairk.backend.kmer_alignment.esm_embedding_distance import (
    pairk_alignment_embedding_distance,
)
from pairk.backend.kmer_alignment.scoring_matrix import pairk_alignment
from pairk.backend.tools.pairwise_tools import PairkAln
from pairk.backend.tools.matrices import print_available_matrices
from pairk.backend.tools.esm_tools import ESM_Model
