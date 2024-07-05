"""motif conservation in IDRs through pairwise k-mer alignment"""

# Add imports here
# from pairk.pairk_aln import (
#     pairk_alignment,
#     pairk_alignment_needleman,
#     pairk_alignment_embedding_distance,
#     make_aligner,
#     PairkAln,
#     print_available_matrices,
#     ESM_Model,
# )
# from pairk.pairk_conservation import (
#     PairkConservation,
#     calculate_conservation,
#     calculate_conservation_arrays,
#     capra_singh_functions,
# )
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

from pairk.backend.conservation.kmer_conservation import (
    PairkConservation,
    calculate_conservation,
    calculate_conservation_arrays,
)
import pairk.backend.conservation.capra_singh_functions as capra_singh_functions

from ._version import __version__
from pathlib import Path
from pairk.backend.tools.sequence_utils import (
    FastaImporter,
    strip_dashes_from_sequences,
)
from pairk.examples import example1
import pairk.utilities as utilities


# Define the __all__ variable
__all__ = [
    "PairkConservation",
    "capra_singh_functions",
    "calculate_conservation_arrays",
    "calculate_conservation",
    "PairkAln",
    "ESM_Model",
    "pairk_alignment_embedding_distance",
    "make_aligner",
    "pairk_alignment_needleman",
    "pairk_alignment",
    "print_available_matrices",
    "FastaImporter",
    "__version__",
]
