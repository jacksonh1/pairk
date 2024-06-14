# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: pairk
#     language: python
#     name: python3
# ---

# %% [markdown]
# Pairk - Pairwise k-mer alignment 

# %%
import pairk
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # load an example dataset

# %% [markdown]
# load an example alignment

# %%
ex1 = pairk.example1

# %% [markdown]
# This example is imported as an object but just contains: an MSA, the query id, positions of query idrs, and the IDR sequences from the MSA themselves.

# %%
print(ex1)

# %% [markdown]
# For Pairk, you just need a set of IDR sequences in dictionary format.<br>
#
# Here I defined the IDRs via the query sequence IDR positions, extracted the relevant sequences from the MSA, and then unaligned them to get the homologous idrs. However, you can generate those IDRs in any way you like.

# %%
import pairk.backend.tools.sequence_utils as tools
import pairk.examples
faimporter = tools.FastaImporter(pairk.examples.example1_alignment_file)
alignment = faimporter.import_as_alignment()
idrs = tools.strip_dashes_from_sequences(list(alignment[:, 555:971]))
idr_dict = {i.id: str(i.seq) for i in idrs}
query_id = "9606_0:00294e"

# %% [markdown]
# # step 1: Pairk alignment

# %% [markdown]
# ## scoring matrix
#
# These methods use a scoring matrix to score the query k-mer to homolog k-mer matches.

# %% [markdown]
# ### method - exhaustive scoring - slow and not recommended for general use

# %% [markdown]
# The exhaustive scoring method actually scores all possible query k-mer to ortho k-mer matches. It's a brute force method that is generally slower than the needleman variant (see the next method). The results should actually be the same as the needleman method. I've included it in the `pairk` package in case people want to use the exhaustive method for custom analysis, such as including multiple high-scoring matches from each ortholog, etc.

# %% [markdown]
# all you need to run the method is:
# - the idr sequences in dictionary format (sequence id as key, sequence string as value)
# - the sequence id of the query sequence
# - the k-mer length to use for the alignment (k)
# - the scoring matrix to use for the alignment (default in EDSSMat50)

# %%
aln_results = pairk.pairk_alignment(
    idr_dict_in=ex1.idr_dict, 
    query_id=ex1.query_id,
    k=5, 
    matrix_name="EDSSMat50"
)

# %% [markdown]
# The results are returned as a `PairkAln` object.

# %%
type(aln_results)

# %%
# pairk.PairkAln?

# %% [markdown]
# The results of the pairwise alignments are stored in pandas DataFrames, which can be directly accessed via the `PairkAln` object.

# %%
aln_results.orthokmer_matrix.head()

# %%
aln_results.orthokmer_matrix.loc[1]

# %%
aln_results.score_matrix.head()

# %% [markdown]
# You can get the "pseudo-alignment" of any query k-mer via the `get_pseudo_alignment` method. <br>This method returns a list of the best-scoring ortholog k-mers for a query k-mer. The query k-mer is specified by its position in the query sequence (0-based). 
# <br>The returned list includes the query k-mer sequence

# %%
aln_results.get_pseudo_alignment(1)

# %% [markdown]
# you can search for a specific kmer to get its positions. You can then use the positions to query the matrices.

# %%
aln_results.find_query_kmer_positions('LPPPP')

# %%
aln_results.get_pseudo_alignment(75)

# %%
aln_results.orthokmer_matrix.loc[[75, 113, 127, 157]].T

# %% [markdown]
# Note - the k-mers are defined by position rather than sequence. You could easily make a variant of this method that uses the unique sequences instead. It would make the method slightly faster. <br>The reason that I didn't do this is because I wanted to mimic the LLM embedding version of Pairk, where identical k-mers have different embeddings and thus are treated as different k-mers.<br>Inclusion of duplicate k-mers does alter the final z-scores, so it's something to be aware of.

# %%
aligner = pairk.make_aligner('EDSSMat50')
r = pairk.pairk_alignment_needleman(idr_dict, query_id, 10, aligner=aligner)

# %%
import pairk.single_kmer as pairk_single
pairk_single.pairk_alignment_single_kmer("LPPPP", idr_dict)

# %%

# %%

import pairk
ex1 = pairk.example1
aln_results = pairk.pairk_alignment(
    idr_dict_in=ex1.idr_dict, 
    query_id=ex1.query_id,
    k=5, 
    matrix_name="EDSSMat50"
)

pairkcons = calculate_pairk_conservation(aln_results, cs.property_entropy)



# %%
import pairk
ex1 = pairk.example1
aln_results = pairk.pairk_alignment(
    idr_dict_in=ex1.idr_dict, 
    query_id=ex1.query_id,
    k=5, 
    matrix_name="EDSSMat50"
)

ok_arr, score_arr, z_score_arr = calculate_conservation(aln_results.orthokmer_matrix, cs.property_entropy)
pairk_conservation = PairkConservation(ok_arr, score_arr, z_score_arr)

# %%
pairk_conservation.plot_background_distribution()
plt.show()


# %%
pairk_conservation.print_array_hist(pairk_conservation.score_arr, np.linspace(0, 1, 30))
pairk_conservation.print_array_hist(pairk_conservation.z_score_arr, 30)

# %%
pairk_conservation.write_results_to_file("test.npz")
pkc = PairkConservation.read_results_from_file("test.npz")

# %%
print_array_hist(scores, np.linspace(0, 1, 30))
print_array_hist(z_scores, 30)


# %% [markdown]
# # Archive


# %%

import pairk.backend.kmer_matrix_scoring.conservation_tools.capra_singh_2007_scores as cs
import pandas as pd
from typing import Callable
import numpy as np
# %%

aln_results.orthokmer_matrix
# check that all elements in the orthokmer_matrix are strings

if not all([type(i) == str for i in aln_results.orthokmer_matrix.values.flatten()]):
    raise ValueError("Orthokmer matrix contains non-string elements")

# I think that numpy will be faster and it's pretty straightforward in this case.
orthokmer_arr = aln_results.orthokmer_matrix.values
orthokmer_arr[0, :]

# %%
def kmerlist_to_columns(seqlist: np.ndarray):
    col_strings = []
    for c in range(len(seqlist[0])):
        col = "".join([seqlist[i][c] for i in range(len(seqlist))])
        col_strings.append(col)
    return col_strings

kmerlist_to_columns(orthokmer_arr[0, :])

# %%
def score_pseudo_aln(
    seqlist: np.ndarray,
    score_func: Callable = cs.property_entropy,
):
    col_strings = kmerlist_to_columns(seqlist)
    return [score_func(c) for c in col_strings]

score_pseudo_aln(orthokmer_arr[0, :])

# %%

def orthokmer_df_to_score_arr(orthokmer_df: pd.DataFrame, score_func: Callable) -> np.ndarray:
    orthokmer_arr = orthokmer_df.copy().values
    k = len(orthokmer_arr[0,0])
    score_arr = np.zeros((orthokmer_arr.shape[0], k))
    for i in range(orthokmer_arr.shape[0]):
        score_arr[i, :] = score_pseudo_alns(orthokmer_arr[i, :])
    return score_arr

score_arr = orthokmer_df_to_score_arr(aln_results.orthokmer_matrix, cs.property_entropy)
bg_scores = score_arr.flatten()
n_bg_scores = len(bg_scores)
u = bg_scores.mean()
s = bg_scores.std()
z_score_arr = (score_arr - u) / s

# %%
# print a text representation of a histogram of the z-scores
def print_array_hist(a: np.ndarray, bins=30):
    hist, bin_edges = np.histogram(a.flatten(), bins=bins)
    for i, h in enumerate(hist):
        print(f"{bin_edges[i]:.2f} - {bin_edges[i+1]:.2f}: {'=' * h}")

print_array_hist(score_arr, np.linspace(0, 1, 30))
print_array_hist(score_arr, 10)




# %%

# convert the score matrix to a z-score matrix
def score_to_zscore(score_arr: np.ndarray, bg_scores: np.ndarray) -> np.ndarray:
    zscore_arr = np.zeros(score_arr.shape)
    for i in range(score_arr.shape[0]):
        zscore_arr[i, :] = (score_arr[i, :] - bg_scores.mean()) / bg_scores.std()
    return zscore_arr























