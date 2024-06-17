"""primary user-facing functions/objects for scoring the conservation of k-mers from the results of the alignment step of the pairk package"""

from pairk.backend.conservation.kmer_conservation import (
    PairkConservation,
    calculate_conservation,
    calculate_conservation_arrays,
)

from typing import Callable

# from pairk.backend.conservation.capra_singh_functions.capra_singh_2007_scores import (
#     property_entropy as capra_singh_property_entropy,
# )
# from pairk.backend.conservation.capra_singh_functions.capra_singh_2007_scores import (
#     shannon_entropy as capra_singh_shannon_entropy,
# )
import pairk.backend.conservation.capra_singh_functions as capra_singh_functions
