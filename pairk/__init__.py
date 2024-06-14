"""motif conservation in IDRs through pairwise k-mer alignment"""

# Add imports here
from pairk.pairk_aln import (
    pairk_alignment,
    pairk_alignment_needleman,
    make_aligner,
    PairkAln,
    print_available_matrices,
)
from pairk.pairk_conservation import (
    PairkConservation,
    calculate_pairk_conservation,
    calculate_conservation,
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
    "make_aligner",
    "PairkAln",
    "PairkConservation",
    "calculate_pairk_conservation",
    "calculate_conservation",
    "__version__",
]
