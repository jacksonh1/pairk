"""
Unit and regression test for the pairk package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import pairk

TEST_FULL_LENGTH_SEQ_DICT = {
    "seq1": "MSTVSDHEMATAA",
    "seq2": "TSTIETSVVVMVSAADDFEEEFDSA",
}
TEST_IDR_POSITION_MAP = {"seq1": [-1, -1], "seq2": [5, 10]}  # no idr
TEST_IDR_DICT = {
    k: v[TEST_IDR_POSITION_MAP[k][0] : TEST_IDR_POSITION_MAP[k][1] + 1]
    for k, v in TEST_FULL_LENGTH_SEQ_DICT.items()
}


def test_pairk_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "pairk" in sys.modules


# def test_encoding_functions():
#     mod = pairk.ESM_Model()
#     test_ids, test_seqs = zip(*TEST_FULL_LENGTH_SEQ_DICT.items())
#     test_embeddings = mod.encode_multiple_seqs(test_seqs)
#     print(test_embeddings.shape)
#     test_single_embeddings = [mod.encode(seq) for seq in test_seqs]
#     for i in test_single_embeddings:
#         print(i.shape)
#     assert len(test_embeddings) == len(test_ids)
#     assert len(test_single_embeddings) == len(test_ids)
#     for i in range(len(test_ids)):
#         print(test_embeddings[i].shape)
#         print(test_single_embeddings[i].shape)
#         t = test_embeddings[i] == test_single_embeddings[i]
#         assert t.all()


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
