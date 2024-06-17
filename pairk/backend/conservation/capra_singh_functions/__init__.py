"""Module to hold conservation methods from Capra and Singh 2007, DOI: 10.1093/bioinformatics/btm270

The methods in this module are used to calculate column-wise conservation scores and take a string of characters (an alignment column) as input.
Conservation scores are calculated and normalized to a range of 0 to 1, where 0 is the least conserved and 1 is the most conserved.
By default, gaps are penalized by multiplying the final conservation score by the fraction of non-gap characters in the column.

Methods:
--------
property_entropy : Callable
    Calculate the property entropy of a column of sequence characters
shannon_entropy : Callable
    Calculate the Shannon entropy of a column of sequence characters
"""

from typing import Callable
from pairk.backend.conservation.capra_singh_functions.capra_singh_2007_scores import (
    property_entropy,
)
from pairk.backend.conservation.capra_singh_functions.capra_singh_2007_scores import (
    shannon_entropy,
)


__all__ = ["property_entropy", "shannon_entropy"]


# from pairk.backend.conservation.capra_singh_functions.capra_singh_2007_scores import (
#     property_entropy,
#     shannon_entropy,
# )

# # __all__ = ["property_entropy","shannon_entropy"]


# class CapraSinghMethods:
#     """Class to hold conservation methods that operate on one column of sequence characters

#     Attributes:
#     ----------
#     property_entropy : Callable
#         Calculate the property entropy of a column of sequence characters. From Capra and Singh 2007, DOI: 10.1093/bioinformatics/btm270.
#     shannon_entropy : Callable
#         Calculate the Shannon entropy of a column of sequence characters. From Capra and Singh 2007, DOI: 10.1093/bioinformatics/btm270
#     """

#     def __init__(self):
#         self.property_entropy = property_entropy
#         self.shannon_entropy = shannon_entropy

#     def __getmethod__(self, key: str) -> Callable:
#         return getattr(self, key)


# capra_singh_scores = CapraSinghMethods()
