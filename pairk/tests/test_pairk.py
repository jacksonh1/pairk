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
