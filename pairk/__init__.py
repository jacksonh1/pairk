"""motif conservation in IDRs through pairwise k-mer alignment"""

# Add imports here
from .pairk_aln import (
    pairk_alignment,
    pairk_alignment_needleman,
    make_aligner,
    PairkAln,
)
from ._version import __version__
from importlib.resources import files as _files
from pathlib import Path

example_alignment_file = Path(
    _files("pairk.data")
    .joinpath("example_alignment_9606_0_00294e-idraln-555_to_971-idr-440_to_665.fasta")
    .__fspath__()  # type: ignore
)

# Define the __all__ variable
__all__ = [
    "pairk_alignment",
    "pairk_alignment_needleman",
    "make_aligner",
    "PairkAln",
    "__version__",
]
