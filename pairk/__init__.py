"""motif conservation in IDRs through pairwise k-mer alignment"""

# Add imports here
from .pairk_aln import pairk_alignment, pairk_alignment_needleman
from ._version import __version__


# Define the __all__ variable
__all__ = [
    "pairk_alignment",
    "pairk_alignment_needleman",
    "__version__",
]
