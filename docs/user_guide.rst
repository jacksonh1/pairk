========================
pairk Usage Guide
========================

.. generate a table of contents for this document
.. contents:: Table of Contents
    :depth: 5
    :local:

see the `tutorial notebook <https://github.com/jacksonh1/pairk/blob/main/demo/demo.ipynb>`_ for an interactive version of the guide (most of the code examples are taken from the notebook), where you can see the outputs of the code cells, the plots, and run the code yourself. The docstrings for the functions/classes look a bit better here though.

********
Overview
********

.. image:: ./images/fragpair_cartoon_v2.png
    :align: center
    :width: 800


Pairk is a Python package for quantifying the conservation of motifs within intrinsically disordered regions (IDRs). It will run on any input protein sequences that you provide but it is designed for disordered regions, where multiple sequence alignments are particularly difficult to interpret and unreliable

the pairk method can be split into 2 steps, where each step has configurable options:

`Step 1: pairwise k-mer alignment`_
- aligns all of the k-mers in the query sequence to all of the k-mers in the homolog sequences. For each query k-mer, the best matching k-mer from each homolog is selected and recorded. The results are returned in a `PairkAln` object which provides methods for accessing the results and writing/reading the results to/from files.

`Step 2: k-mer conservation`_
- calculates conservation scores for each query k-mer using the results of Step 1. The conservation score results are stored in a `PairkConservation` object which also provides methods for plotting the results as well as writing/reading the results to/from files.


************
Installation
************

To install the package, use the following command:

.. code-block:: bash

    pip install pairk

or for an editable install that you can modify:

.. code-block:: bash

    git clone https://github.com/jacksonh1/pairk.git
    cd pairk
    pip install -e .


*****************
Detailed guide
*****************

Most of the functionality of pairk is accessed through the main pairk module, which can be imported directly.

**Import pairk**

.. code-block:: python

    import pairk

pairk comes with some example data that we can use to test the methods. The example data is stored in the ``pairk.example1`` object.

.. code-block:: python

    ex1 = pairk.example1

This example is imported as a python object and holds data that is compatible with the main pairk methods

for example, we can access the IDR sequences in a dictionary with ``ex1.idr_dict``
or the query id with ``ex1.query_id``

.. code-block:: python

    for id, seq in ex1.idr_dict.items():
        print(id, seq)

    print(ex1.query_id)


Step 1: pairwise k-mer alignment
=================================

there are 2 main ways to run the k-mer alignment step:

1. `Scoring matrix alignment`_

   * *These methods use a scoring matrix to score the query k-mer to homolog k-mer matches and select the best scoring match from each homolog.*

2. `Embedding distance alignment`_

   * *This method uses the Euclidean distance between the query k-mer residue embeddings from a protein large language model (such as ESM2) and homolog k-mer residue embeddings and selects the lowest distance match from each homolog.*


Scoring matrix alignment
-------------------------------------

There are 2 implementations of the scoring matrix method:

1. ``pairk.pairk_alignment`` - the original implementation. This is a bit slow because it does an exhaustive comparison of all k-mers in the query sequence with all k-mers in the homologs.
2. ``pairk.pairk_alignment_needleman`` - a faster implementation that uses the Needleman-Wunsch algorithm (as implemented in Biopython) to align the k-mers. This is faster and should yield the same results.


The basic inputs for these functions are:

* ``idr_dict_in``: a dictionary of IDR sequences, where the keys are the sequence ids and the values are the sequences. Includes the query sequence (the sequence to split into k-mers and align with the homologs).
* ``query_id``: a query sequence id (the sequence to split into k-mers and align with the homologs). This id should be present in ``idr_dict_in``.
* ``k``: the length of the k-mers

These inputs can be generated many ways, and there are helper functions in the pairk library to help with this (coming soon). 


pairk.pairk_alignment - slower
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. autofunction:: pairk.pairk_alignment
    :no-index:

example usage:

.. code-block:: python

    aln_results = pairk.pairk_alignment(
        idr_dict_in=ex1.idr_dict,
        query_id=ex1.query_id,
        k=5,
        matrix_name="EDSSMat50"
    )


To specify the scoring matrix used, you can pass the name of the matrix to the ``matrix_name`` argument.
To see the available matrices, use the ``pairk.print_available_matrices()`` function.

.. code-block:: python

    pairk.print_available_matrices()



k-mer alignment results
"""""""""""""""""""""""

The results are returned as a ``PairkAln`` object.

The actual "alignments" are stored as matrices in the ``PairkAln`` object. The main matrices are:

* orthokmer_matrix - the best matching k-mers from each homolog for each query k-mer
* position_matrix - the positions of the best matching k-mers in the homologs
* score_matrix - the scores of the best matching k-mers

Each matrix is a pandas DataFrame where the index is the start position of the k-mer in the query sequence. The columns are the query k-mers + the homolog sequence ids.

The ``PairkAln`` object has some useful methods for accessing the data. For example, you can get the best matching k-mers for a query k-mer by its position in the query sequence using the ``.get_pseudo_alignment`` method (or by directly accessing the dataframes). You can also plot the matrices as heatmaps, save the results to a json file, and load the results from that file. examples are shown below

.. autoclass:: pairk.PairkAln
   :members:
   :undoc-members:
   :no-index:

example: accessing the DataFrames from the ``PairkAln`` object directly

.. code-block:: python

    aln_results.score_matrix

example: access the best matching k-mers for the query k-mer at position 4:

.. code-block:: python

    aln_results.orthokmer_matrix.loc[4]

example: access the best matching k-mers for the query k-mer at position 4 using the ``get_pseudo_alignment`` method.
(the returned list includes the query k-mer sequence)

.. code-block:: python

    aln_results.get_pseudo_alignment(4)

example: plot a heatmap of the matrices:

.. code-block:: python

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3,3))
    aln_results.plot_position_heatmap(ax)
    ax.xaxis.set_visible(False)

example: save the results to a file using ``write_to_file`` and load them back into python using ``from_file``:

.. code-block:: python

    aln_results.write_to_file('./aln_results.json')
    aln_results = pairk.PairkAln.from_file('./aln_results.json')
    print(aln_results)


pairk.pairk_alignment_needleman - faster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This method returns the same results as ``pairk.pairk_alignment``, but it is faster.

The difference is that the ``pairk.pairk_alignment_needleman`` method uses the Needleman-Wunsch algorithm to align the k-mers, while the ``pairk.pairk_alignment`` method uses a scoring matrix to exhaustively score the k-mer matches. ``pairk.pairk_alignment_needleman`` ensures that the alignment is gapless by using an extremely high gap opening and extension penalty (-1000000). This will ensure that the alignment is gapless, unless you use a really unusual scoring matrix with very high scores.

This method takes similar arguments as ``pairk.pairk_alignment``, accept that the ``pairk.pairk_alignment_needleman`` method takes an optional ``aligner`` argument. This allows you to create the aligner before calling the method, which is useful if you want to do multiprocessing, so that you're not creating a new aligner for each process. I've found that if you create a new aligner for each process, the memory usage gets very high, as if the memory isn't being released until the entire script finishes

The ``aligner`` object can be created via the ``pairk.create_aligner`` function. This function takes the name of the scoring matrix as an argument and returns the aligner object. If you don't pass the ``aligner`` argument to the ``pairk.pairk_alignment_needleman`` method, it will create a new aligner using the ``matrix_name`` argument. This is fine if you're not doing multiprocessing. If you are doing multiprocessing, I would suggest creating the aligner before calling the method. If the ``aligner`` argument is passed, the ``matrix_name`` argument is ignored.

.. autofunction:: pairk.pairk_alignment_needleman
    :no-index: 

.. code-block:: python

    # Making the aligner ahead of time to demonstrate
    aligner = pairk.make_aligner('EDSSMat50')

.. code-block:: python

    aln_results_needleman = pairk.pairk_alignment_needleman(
        idr_dict_in=ex1.idr_dict,
        query_id=ex1.query_id,
        k=5,
        aligner=aligner
    )

Results are the same as the previous method:

.. code-block:: python

    (aln_results.position_matrix == aln_results_needleman.position_matrix).all().all()


Embedding distance alignment
------------------------------

This method uses the Euclidean distance between the query k-mer residue embeddings and homolog k-mer residue embeddings and selects the lowest distance match from each homolog. For each homolog, it calculates the distance between the query k-mer and each k-mer in the homolog. It then selects the k-mer with the lowest distance as the best match.


pairk.pairk_alignment_embedding_distance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This method generates residue embeddings using the ESM2 protein large language model. 

Because residue embeddings are used, the inputs are slightly different than the previous methods. 

The inputs are:

* ``full_length_dict_in``: a dictionary of full-length sequences, where the keys are the sequence ids and the values are the sequences. This is used to generate the embeddings.
* ``idr_position_map``: a dictionary where the keys are the full-length sequence ids and the values are the start and end positions of the IDR in the full-length sequence (using python indexing). This is used to slice out the IDR embeddings/sequences from the full-length embeddings/sequences.
* ``query_id``: a query sequence id (the sequence to split into k-mers and align with the homologs). This id should be present in ``idr_position_map`` and ``full_length_dict_in``.
* ``k`` - the length of the k-mers
* ``mod`` - a ``pairk.ESM_Model`` object. This is the ESM2 model used to generate the embeddings. The code for the ESM2 embeddings is adapted from the kibby conservation tool [link](https://github.com/esbgkannan/kibby) DOI: 10.1093/bib/bbac599
* ``device`` - whether to use cuda or your cpu for pytorch, should be either "cpu" or "cuda". (default is "cuda"). If "cuda" fails, it will default to "cpu". This argument is passed to the ``pairk.ESM_Model.encode`` method

The ``mod`` input is required so that you can preload the ESM model before running the method. <br><br>

Full length sequences (``full_length_dict_in``) are required to generate the embeddings because each embedding is dependent upon the neighboring residues. The embeddings for just an IDR are different than the embeddings for a full length sequences. Thus, the full length embeddings are gathered first, and then the IDR embeddings are sliced out for the k-mer alignment. 

The ``idr_position_map`` is used to slice out the IDR embeddings, and there must be IDR positions for each sequence in the input sequence set.

There is currently no way to use pre-generated embeddings for this method, but this functionality would be very easy to add.

.. autofunction:: pairk.pairk_alignment_embedding_distance
    :no-index: 

The ``pairk.pairk_alignment_embedding_distance`` method returns a ``PairkAln`` object, just like the previous methods

.. autoclass:: pairk.ESM_Model
   :members:
   :undoc-members:
   :no-index: 


example usage: loading the ESM2 model and running the method

.. code-block:: python

    mod = pairk.ESM_Model()
    aln_results_embedding = pairk.pairk_alignment_embedding_distance(
        full_length_dict_in=ex1.full_length_dict,
        idr_position_map=ex1.idr_position_map,
        query_id=ex1.query_id,
        k=5,
        mod=mod,
        device="cpu"
    )


Step 2: k-mer conservation
=================================

In this step, the query k-mer and the best matching homolog k-mers are treated as a gapless multiple sequence alignment with 'k' columns, which we call a "pseudo-MSA". Column-wise conservation scores are calculated for each position in each pseudo-MSA. All of the conservation scores are then converted to z-scores to give the relative conservation of each k-mer position compared to the rest of the query IDR. The conservation score results are stored in a ``PairkConservation`` object which also provides methods for plotting the results and reading/writing the results from/to files.


pairk.calculate_conservation
------------------------------

the main method for Step 2 is the ``pairk.calculate_conservation`` method. It simply takes the ``PairkAln`` object as input, along with a columnwise conservation scoring function and returns a ``PairkConservation`` object.

.. autofunction:: pairk.calculate_conservation
   :no-index: 
 
The columnwise conservation scoring function can be any function that takes a string of residues (a column of an alignment) as an input and returns a float (conservation score). You can use custom functions here, but pairk comes with a few built-in functions from Capra and Singh 2007 (DOI: 10.1093/bioinformatics/btm270) available in the ``pairk.pairk_conservation.capra_singh_functions`` module. The ``pairk.pairk_conservation.capra_singh_functions.property_entropy`` is the default function used by ``pairk.calculate_conservation``.

.. autoclass:: pairk.pairk_conservation.capra_singh_functions
   :members:
   :undoc-members:
   :no-index: 


example usage:

.. code-block:: python

    aln_results = pairk.pairk_alignment(
        idr_dict_in=ex1.idr_dict,
        query_id=ex1.query_id,
        k=5,
        matrix_name="EDSSMat50"
    )
    conservation_results = pairk.calculate_conservation(
        aln_results,
    )

example usage: using a different conservation scoring function:

.. code-block:: python

    from pairk.pairk_conservation import capra_singh_functions
    column = 'NNNNNNNNNKNSNNNNNNNNSSN'
    print(capra_singh_functions.shannon_entropy(column))

.. code-block:: python

    aln_results = pairk.pairk_alignment(
        idr_dict_in=ex1.idr_dict,
        query_id=ex1.query_id,
        k=5,
    )
    conservation_results = pairk.calculate_conservation(
        pairk_aln_results=aln_results,
        score_func=capra_singh_functions.shannon_entropy
    )


k-mer conservation results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``pairk.calculate_conservation`` method returns a ``PairkConservation`` object.
The returned ``PairkConservation`` object has matrices with similar structure as ``PairkAln`` object matrices, except that they are numpy arrays instead of pandas dataframes.

* orthokmer_arr - the best matching k-mers from each homolog for each query k-mer - analogous to the orthokmer_matrix in the ``PairkAln`` object
* score_arr - the conservation scores for each position in the pseudo-MSA of each query k-mer
* z_score_arr - the conservation score z-scores for each position in the pseudo-MSA of each query k-mer

If ``n`` is the number of k-mers in the query sequence and ``m`` is the number of homologs, the matrices will have the dimensions:

* orthokmer_arr: (n, m)
* score_arr: (n, k)
* z_score_arr: (n, k)

The row index of the arrays correspond to the starting position of the query k-mer in the query IDR. For example, to access the conservation scores for the k-mer at position 4 in the query IDR, you would access the 4th row of the arrays: ``.score_arr[4, :]``.

.. autoclass:: pairk.PairkConservation
   :members:
   :undoc-members:
   :no-index: 


accessing the results
"""""""""""""""""""""""

You can use the k-mer starting position in the query IDR to access the k-mers, scores, and z-scores for each k-mer.

For example, for the k-mer at position 4:

.. code-block:: python

    k_mer_position = 4
    print(f"query k-mer at position {k_mer_position}: {conservation_results.orthokmer_arr[k_mer_position, 0]}")
    print(f"pseudo-MSA for the query k-mer at position {k_mer_position} (including the query k-mer): {conservation_results.orthokmer_arr[k_mer_position, :]}")
    print(f"scores for each position of the k-mer at position {k_mer_position}:")
    print(conservation_results.score_arr[k_mer_position, :])
    print(f"z scores for each position of the k-mer at position {k_mer_position}:")
    print(conservation_results.z_score_arr[k_mer_position, :])


plotting the results
"""""""""""""""""""""""
There are several plotting functions available from the ``pairk.PairkConservation`` object shown below.

example usage: plotting background score distributions

.. code-block:: python

    conservation_results.plot_background_distribution()

example usage: plotting conservation scores

.. code-block:: python

    fig, ax = plt.subplots(figsize=(7,1.5))
    conservation_results.plot_score_barplot(k_mer_position, ax=ax)
    ax.set_title('conservation scores')
    fig, ax = plt.subplots(figsize=(7,1.5))
    conservation_results.plot_score_barplot(k_mer_position, score_type='z_score', ax=ax)
    ax.set_title('conservation z-scores')

example usage: display a pseudo-MSA as a sequence logo

.. code-block:: python

    fig, ax = plt.subplots(figsize=(7,1.5))
    conservation_results.plot_sequence_logo(k_mer_position, ax=ax)


example usage: plotting a conservation summary plot

.. code-block:: python

    fig, axd = conservation_results.plot_conservation_mosaic(
        position=0,
        score_type='z_score',
        figsize=(9, 3)
    )
    plt.tight_layout(h_pad=0.5, w_pad=0.5)


getting the average conservation score for a query k-mer
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

You can use the ``get_average_score`` function to get the average conservation score for a k-mer position.

example: get the average conservation score for the query k-mer at position 4

.. code-block:: python

    conservation_results.get_average_score(4, score_type='z_score')

The ``get_average_score`` function takes a ``position_mask`` as an optional argument that will only consider the conservation scores for the positions in the mask when calculating the average score. This is useful if you want to exclude certain positions from the average score calculation.

example: get the average conservation score for the query k-mer at position 0, but only consider the conservation scores for positions 1, and 3 within the k-mer

.. code-block:: python

    position_mask = [0, 1, 0, 1, 0]
    conservation_results.get_average_score(0, score_type='z_score', position_mask=position_mask)

You could also do a weighted average from manually extracted conservation scores.

example: get the weighted average conservation score for the query k-mer at position 0 (using some arbitrary weights)

.. code-block:: python

    np.average(conservation_results.z_score_arr[0, :], weights=[0.1, 1, 0.5, 1, 10])


Writing and Reading Results from Files
"""""""""""""""""""""""""""""""""""""""

You can save the results to a file with ``write_results_to_file`` and load them back in with ``read_results_from_file``.

example usage: save the results to a file and load them back in

.. code-block:: python

    conservation_results.write_results_to_file('./conservation_results.npz')

    conservation_results.read_results_from_file('./conservation_results.npz')


advanced customization
=================================
Step 1 and Step 2 can be modified, but this requires a bit of knowledge of the inner workings of the pairk package and the source code would have to be modified directly and probably installed in editable mode (see `Installation`_). If you want to modify the package, here's roughly how I organized the code in the ``pairk`` directory (the source code is available on `github <https://github.com/jacksonh1/pairk/>`_). Below file paths are relative to the github repo root directory:

* ``pairk/backend/`` - the main code for the package is located here. The main pairwise k-mer alignment and k-mer conservation functions are defined in files within this directory.
* ``pairk/__init__.py`` - User-facing functions are imported into the main ``pairk/__init__.py`` file so that they are accessible when the package is imported. I think it also simplifies the import statements for users. Use this __init__ file to find where pairk's main functions are defined within the directory structure if you want to modify any of the functions above. You could also modify the __init__ file to make any new functions you create easy to access.
* ``pairk/data/`` - data installed along with the package is stored here. This includes the scoring matrices and example data. The scoring matrices are stored in the ``pairk/data/matrices/`` folder.

The easiest customization to make would be to add a new scoring matrix. To do this, you would add a new matrix file to the ``pairk/data/matrices/`` folder. The tools should be able to automatically find the matrix file and add it to the available scoring matrices. It will be named after the name of the file. Use ``pairk.print_available_matrices()`` to confirm (make sure you've installed pairk as an editable install for changes to take affect). You could then use the new matrix in relevant methods by passing the matrix name as an argument. If this doesn't work, you may need to modify the code that reads the matrices in the ``pairk.backend.tools.matrices.py`` file.


