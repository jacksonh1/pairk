"""motif conservation in IDRs through pairwise k-mer alignment"""

# Add imports here
from pairk.pairk_aln import (
    pairk_alignment,
    pairk_alignment_needleman,
    pairk_alignment_embedding_distance,
    make_aligner,
    PairkAln,
    print_available_matrices,
    ESM_Model,
)
from pairk.pairk_conservation import (
    PairkConservation,
    calculate_conservation,
    calculate_conservation_arrays,
    capra_singh_functions,
)
from ._version import __version__
from pathlib import Path
from pairk.backend.tools.sequence_utils import (
    FastaImporter,
    strip_dashes_from_sequences,
)
from pairk.examples import example1


# Define the __all__ variable
__all__ = [
    "FastaImporter",
    "print_available_matrices",
    "pairk_alignment",
    "pairk_alignment_needleman",
    "pairk_alignment_embedding_distance",
    "make_aligner",
    "PairkAln",
    "PairkConservation",
    "calculate_conservation",
    "calculate_conservation_arrays",
    "capra_singh_functions",
    "ESM_Model",
    "__version__",
]
