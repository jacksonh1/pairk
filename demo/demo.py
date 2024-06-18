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
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Pairk - tutorial notebook

# %% [markdown]
# Pairk - quantifying the conservation of motifs within intrinsically disordered regions.
#
# This notebook will go over the basics of how to use the pairk library using an example set of sequences (comes pre-installed with pairk).

# %% [markdown]
# ## brief overview of pairk
#
# Pairk can be divided into 2 parts:
# 1. `pairk_aln` - pairwise kmer alignment. Aligns each k-mer in a query sequence with each sequence in a set of homologs. The alignment is gapless and performed pairwise. For each k-mer in the query sequence, the result is a set of the best matching k-mers from each homolog. We call this set of k-mers a "pseudo-MSA". The k-mer alignment step can be performed via a scoring matrix (similar to a traditional alignment), or via ESM2 embeddings. These are explained in more detail below.
# 2. `pairk_conservation` - scores the conservation of the pseudo-MSAs via columnwise conservation scores.

# %% [markdown]
# # Step 1. pairk_aln

# %% [markdown]
# The basic input for step 1 is:
# - `idr_dict_in` - a dictionary of IDR sequences, where the keys are the sequence ids and the values are the sequences. Includes the query sequence (the sequence to split into k-mers and align with the homologs).
# - `query_id` - a query sequence id (the sequence to split into k-mers and align with the homologs). This id should be present in `idr_dict_in`.
# - `k` - the length of the k-mers

# %% [markdown]
# These inputs can be generated many ways, and there are helper functions in the pairk library to help with this. For this example, we will use the example sequences that come with pairk.

# %% [markdown]
# ## import pairk

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import pairk

# %% [markdown]
# ## load an example dataset

# %%
ex1 = pairk.example1

# %% [markdown]
# This example is imported as an object but it just holds the arguments that we need for `pairk_aln`:<br> 
# - the IDR sequences in a dictionary - `ex1.idr_dict`<br>
# - the query id - `ex1.query_id`<br>

# %%
for id,seq in ex1.idr_dict.items():
    print(id,seq)

# %%
print(ex1.query_id)

# %% [markdown]
# ## the scoring matrix method
#
# These methods use a scoring matrix to score the query k-mer to homolog k-mer matches and select the best scoring match from each homolog.
#
# There are 2 implementations of the scoring matrix method:
# 1. `pairk.pairk_alignment` - the original implementation. This is a bit slow because it does an exhaustive comparison of all k-mers in the query sequence with all k-mers in the homologs.
# 2. `pairk.pairk_alignment_needleman` - a faster implementation that uses the Needleman-Wunsch algorithm (as implemented in Biopython) to align the k-mers. This is faster and should yield the same results.
#
#
# To specify the scoring matrix used, you can pass the name of the matrix to the `matrix_name` argument. <br>
# To see the available matrices, use the `pairk.print_available_matrices()` function.
#
#
# <font size="2">*Note*: These methods currently find the best scoring match from each homolog for each query k-mer. However, if there are multiple top-scoring matches, only one is returned. The `pairk.pairk_alignment` method could be modified to return all top-scoring matches if needed, however that is not currently implemented. This may be more difficult to implement for the `pairk.pairk_alignment_needleman` method, though it should still be possible since the Biopython method (`Bio.Align.PairwiseAligner.align`) returns many alignments. The relevant code files are i `pairk/backend/kmer_alignment/` if you want to try to implement this yourself.</font>

# %%
pairk.print_available_matrices()

# %% [markdown]
# ---

# %% [markdown]
# ### `pairk.pairk_alignment` - slower of the two methods

# %%
aln_results = pairk.pairk_alignment(
    idr_dict_in=ex1.idr_dict, 
    query_id=ex1.query_id,
    k=5, 
    matrix_name="EDSSMat50"
)

# %% [markdown]
# #### `pairk_aln` results

# %% [markdown]
# The results are returned as a `PairkAln` object.<br>
#
# The actual "alignments" are stored as matrices in the `PairkAln` object. The main matrices are:<br>
# - orthokmer_matrix - the best matching k-mers from each homolog for each query k-mer<br>
# - position_matrix - the positions of the best matching k-mers in the homologs<br>
# - score_matrix - the scores of the best matching k-mers<br>
#
# Each matrix is a pandas DataFrame where the index is the start position of the k-mer in the query sequence. The columns are the query k-mers + the homolog sequence ids.

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

# %%
for kmer, score in zip(aln_results.orthokmer_matrix.loc[0].values[1:], aln_results.score_matrix.loc[0].values[1:]):
    print(kmer, score)

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

# %% [markdown]
# you can also plot a heatmap of the matrices

# %%
import matplotlib.pyplot as plt

# %%
fig, ax = plt.subplots(figsize=(3,3))
aln_results.plot_position_heatmap(ax)
ax.xaxis.set_visible(False)

# %% [markdown]
# You can also save the results to a file using `write_to_file` and load them back into python using `from_file`.<br>
# example:

# %%
aln_results.write_to_file('./aln_results.json')

# %%
aln_results = pairk.PairkAln.from_file('./aln_results.json')
print(aln_results)

# %% [markdown]
# ---

# %% [markdown]
# ### `pairk.pairk_alignment_needleman` - faster

# %% [markdown]
# This method returns the same results as `pairk.pairk_alignment`, but it is faster.
#
# The difference is that the `pairk.pairk_alignment_needleman` method uses the Needleman-Wunsch algorithm to align the k-mers, while the `pairk.pairk_alignment` method uses a scoring matrix to exhaustively score the k-mer matches. It ensures that the alignment is gapless by using an extremely high gap opening and extension penalty (-1000000). This will ensure that the alignment is gapless, unless you use a really unusual scoring matrix with very high scores.
#
# This methods take similar arguments as `pairk.pairk_alignment`, accept that the `pairk.pairk_alignment_needleman` method takes an optional `aligner` argument. This allows you to create the aligner before calling the method, which is useful if you want to do multiprocessing, so that you're not creating a new aligner for each process. I've found that if you create a new aligner for each process, the memory usage gets very high, as if the memory isn't being released until the entire script finishes
#
# The `aligner` object can be created via the `pairk.create_aligner` function. This function takes the name of the scoring matrix as an argument and returns the aligner object.
# If you don't pass the `aligner` argument to the `pairk.pairk_alignment_needleman` method, it will create a new aligner using the `matrix_name` argument. This is fine if you're not doing multiprocessing. If you are doing multiprocessing, i would suggest creating the aligner before calling the method. Or using 1 aligner for each process. If the `aligner` argument is passed, the `matrix_name` argument is ignored.

# %%
# Making the aligner ahead of time to demonstrate
aligner = pairk.make_aligner('EDSSMat50')

# %%
aln_results_needleman = pairk.pairk_alignment_needleman(
    idr_dict_in=ex1.idr_dict, 
    query_id=ex1.query_id,
    k=5, 
    aligner=aligner
)

# %% [markdown]
# results are the same as the previous method

# %%
(aln_results.position_matrix == aln_results_needleman.position_matrix).all().all()

# %% [markdown]
# ## ESM2 embedding distance method

# %% [markdown]
# ### `pairk.pairk_alignment_embedding_distance`

# %% [markdown]
# This method uses the Euclidean distance between the query k-mer residue embeddings and homolog k-mer residue embeddings and selects the lowest distanc match from each homolog.<br>
# For each homolog, it calculates the distance between the query k-mer and each k-mer in the homolog. It then selects the k-mer with the lowest distance as the best match.<br>
#
# Because residue embeddings are used, the inputs are slightly different than the previous methods. The inputs are:<br>
# - `full_length_dict_in` - a dictionary of full-length sequences, where the keys are the sequence ids and the values are the sequences. This is used to generate the embeddings.
# - `idr_position_map` - a dictionary where the keys are the full-length sequence ids and the values are the start and end positions of the IDR in the full-length sequence (using python indexing). This is used slice out the IDR embeddings/sequences from the full length embeddings/sequences.
# - `query_id` - a query sequence id (the sequence to split into k-mers and align with the homologs). This id should be present in `idr_position_map` and `full_length_dict_in`.
# - `k` - the length of the k-mers
# - `mod` - a `pairk.ESM_Model` object. This is the ESM2 model used to generate the embeddings. The code for the ESM2 embeddings is adapted from the kibby conservation tool [link](https://github.com/esbgkannan/kibby) DOI: 10.1093/bib/bbac599
# - `device` - the "device" for pytorch to use to generate the embeddings, either "cpu" or "cuda". (default is "cuda"). If "cuda" fails, it will default to "cpu".
#
# The `mod` input is required so that you can preload the ESM model before running the method. <br><br>
# Full length sequences (`full_length_dict_in`) are required to generate the embeddings because each embedding is dependent upon the neighboring residues. The embeddings for just an IDR are different than the embeddings for a full length sequences. Thus, the full length embeddings are gathered first, and then the IDR embeddings are sliced out for the k-mer alignment. <br><br>The `idr_position_map` is used to slice out the IDR embeddings, and there must be IDR positions for each sequence in the input sequence set.

# %% [markdown]
#

# %% [markdown]
# There is currently no way to use pre-generated embeddings for this method. I may add this in the future.

# %% [markdown]
# #### loading an ESM2 model using `pairk.ESM_Model`

# %%
# pairk.ESM_Model?

# %%
mod = pairk.ESM_Model(threads=4)

# %% [markdown]
# Luckily, I already have a full length sequence dictionary and an IDR position map for the example sequences, so we can use those for this example

# %% [markdown]
# I will just use the cpu for this example, but you would set the device to "cuda" if you have cuda set up on your machine and want to use gpu

# %% [markdown]
# installation and set up of cuda is outside the scope of this documentation/tutorial

# %%
aln_results_embedding = pairk.pairk_alignment_embedding_distance(
    full_length_dict_in=ex1.full_length_dict,
    idr_position_map=ex1.idr_position_map,
    query_id=ex1.query_id,
    k=5,
    mod=mod,
    device='cpu',
)

# %% [markdown]
# ---

# %% [markdown]
# # step 2. pairk_conservation

# %% [markdown]
# The last step is to calculate the conservation scores for the pseudo-MSAs. This is done via the `pairk.calculate_conservation` method. It simply takes the `PairkAln` object as input, along with a columnwise conservation scoring function and returns a `PairkConservation` object.

# %%
# pairk.calculate_conservation?

# %%
conservation_results = pairk.calculate_conservation(aln_results_needleman)
print(conservation_results)

# %% [markdown]
# Any function that takes a string of residues representing a column of a sequence alignment can be used to calculate the conservation scores. <br>By default the property_entropy function from Capra and Singh 2007 (DOI: 10.1093/bioinformatics/btm270) is used.<br>A shannon entropy function is also available in the `pairk.pairk_conservation.capra_singh_functions` module.

# %%
from pairk.pairk_conservation import capra_singh_functions
# capra_singh_functions?

# %% [markdown]
# Example conservation scoring function

# %%
column = 'NNNNNNNNNKNSNNNNNNNNSSN'
print(capra_singh_functions.shannon_entropy(column))

# %%
conservation_results = pairk.calculate_conservation(
    pairk_aln_results=aln_results_needleman,
    score_func=capra_singh_functions.shannon_entropy
)
print(conservation_results)

# %% [markdown]
# The returned `PairkConservation` object has the following attributes:<br>

# %%
# pairk.PairkConservation?

# %% [markdown]
# The matrix structures are essentially the same structure as the `PairkAln` objects, except that they are numpy arrays instead of dataframes. Each row in the arrays correspond to the starting position of the query k-mer in the query idr. The `score_arr` attribute is the conservation scores for each position in each pseudo-MSA. The `z_score_arr` attribute is the z-scores for each position in each pseudo-MSA.

# %% [markdown]
# You can use the k-mer starting position in the query IDR to access the k-mers, scores, and z-scores for each k-mer

# %% [markdown]
# For example, for the k-mer at position 0

# %%
k_mer_position = 0
print(f"k-mer at position {k_mer_position}: {conservation_results.orthokmer_arr[k_mer_position, 0]}")
print(f"scores for each position of the k-mer at position {k_mer_position}:")
print(conservation_results.score_arr[k_mer_position, :])
print(f"z scores for each position of the k-mer at position {k_mer_position}:")
print(conservation_results.z_score_arr[k_mer_position, :])

# %% [markdown]
# ## plotting functions

# %% [markdown]
# There are several plotting functions available from the `pairk.PairkConservation` object shown below

# %%
plt.rcParams.update({'font.size': 12})

# %% [markdown]
# Plotting background score distributions

# %%
conservation_results.plot_background_distribution()

# %% [markdown]
# plotting the conservation scores

# %%
fig, ax = plt.subplots(figsize=(7,1.5))
conservation_results.plot_score_barplot(k_mer_position, ax=ax)
ax.set_title('conservation scores')
fig, ax = plt.subplots(figsize=(7,1.5))
conservation_results.plot_score_barplot(k_mer_position, score_type='z_score', ax=ax)
ax.set_title('conservation z-scores')

# %% [markdown]
# Plotting sequence logos

# %%
fig, ax = plt.subplots(figsize=(7,1.5))
conservation_results.plot_sequence_logo(k_mer_position, ax=ax)

# %% [markdown]
# plotting a conservation summary plot

# %%
fig, axd=conservation_results.plot_conservation_mosaic(
    position = 0,
    score_type='z_score',
    figsize=(9, 3)
)
plt.tight_layout(h_pad=0.5, w_pad=0.5)

# %% [markdown]
# ## average conservation scores

# %% [markdown]
# You can use the `get_average_score` function to get the average conservation score for a k-mer position.

# %%
conservation_results.get_average_score(0, score_type='z_score')

# %% [markdown]
# The `get_average_score` function takes a `position_mask` as an optional argument that will only consider the conservation scores for the positions in the mask when calculating the average score. This is useful if you want to exclude certain positions from the average score calculation.

# %%
position_mask = [0, 1, 0, 1, 0]
conservation_results.get_average_score(0, score_type='z_score', position_mask=position_mask)

# %% [markdown]
# you could also do a weighted average from manually extracted conservation scores

# %%
# weighted average
np.average(conservation_results.z_score_arr[0, :], weights=position_mask)

# %%
np.average(conservation_results.z_score_arr[0, :], weights=[0.1, 1, 0.5, 1, 10])

# %% [markdown]
# ## writing and reading results from files

# %% [markdown]
# you can save the results to a file with `write_results_to_file` and load them back in with `read_results_from_file`

# %%
conservation_results.write_results_to_file('./conservation_results.npz')

# %%
conservation_results.read_results_from_file('./conservation_results.npz')
