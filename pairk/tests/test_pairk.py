"""
Unit and regression test for the pairk package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import pairk


def test_pairk_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "pairk" in sys.modules


# import the example set or define one here
# run the alignment methods and check the outputs
# might want to test:
# - the length of the output dataframes
# - the data types of the output dataframes
# - the contents of the output dataframes - all k-mers length k.
# include:
# an IDR that is shorter than k
# a sequence with an AA that is not in the scoring matrix
# a sequence with no idr

# no idrs at all -> shouldn't give scores
