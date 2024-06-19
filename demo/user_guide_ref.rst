========================
pairk Usage Guide
========================

.. generate a table of contents for this document
.. contents:: Table of Contents
    :depth: 4
    :local:


********
Overview
********

.. image:: ./images/fragpair_cartoon_v2.png
    :align: center
    :width: 800


Pairk is a Python package for quantifying the conservation of motifs within intrinsically disordered regions (IDRs). It will run on any input protein sequences that you provide but it is designed for disordered regions, where multiple sequence alignments are particularly difficult to interpret and unreliable

the pairk method can be split into 2 steps, where each step has configurable options:

step 1: pairwise k-mer alignment
- aligns all of the k-mers in the query sequence to all of the k-mers in the homolog sequences. The best matching k-mer from each homolog is selected for each query k-mer. The results are returned in a `PairkAln` object.

step 2: k-mer conservation
- calculates conservation scores for each query k-mer. The query k-mer and the best matching homolog k-mers are treated as a gapless multiple sequence alignment with 'k' columns, which we call a "pseudo-MSA". Column-wise conservation scores are calculated for each position in each pseudo-MSA. All of the conservation scores are then converted to z-scores to give the relative conservation. The conservation score results are stored in a `PairkConservation` object which also provides methods for plotting the results.


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

This tutorial will guide you through the basics of using the pairk library using an example set of sequences.

**Import pairk**

.. code-block:: python

    import pairk

**Load an example dataset**

For this example, we will use the example sequences that come with pairk.

.. code-block:: python

    ex1 = pairk.example1

This example is imported as a python object and holds the data that we need for inputs to the k-mer alignment methods for convenience.

for example, we can access the IDR sequences in a dictionary with ``ex1.idr_dict``
or the query id with ``ex1.query_id``

.. code-block:: python

    for id, seq in ex1.idr_dict.items():
        print(id, seq)

    print(ex1.query_id)


Step 1: k-mer alignment
=======================

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
- orthokmer_matrix - the best matching k-mers from each homolog for each query k-mer
- position_matrix - the positions of the best matching k-mers in the homologs
- score_matrix - the scores of the best matching k-mers

Each matrix is a pandas DataFrame where the index is the start position of the k-mer in the query sequence. The columns are the query k-mers + the homolog sequence ids.


you can access the dataframes directly. For example, the score matrix:

.. code-block:: python

    aln_results.score_matrix

example - access the best matching k-mers for the query k-mer at position 1:

.. code-block:: python

    aln_results.orthokmer_matrix.loc[1]


You can get the "pseudo-alignment" of any query k-mer via the ``get_pseudo_alignment`` method.
This method returns a list of the best-scoring ortholog k-mers for a query k-mer. The query k-mer is specified by its position in the query sequence (0-based). The returned list includes the query k-mer sequence.

.. code-block:: python

    aln_results.get_pseudo_alignment(1)
    aln_results.find_query_kmer_positions('LPPPP')
    aln_results.get_pseudo_alignment(75)
    aln_results.orthokmer_matrix.loc[[75, 113, 127, 157]].T

You can also plot a heatmap of the matrices:

.. code-block:: python

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(3,3))
    aln_results.plot_position_heatmap(ax)
    ax.xaxis.set_visible(False)

You can also save the results to a file using ``write_to_file`` and load them back into python using ``from_file``:

.. code-block:: python

    aln_results.write_to_file('./aln_results.json')
    aln_results = pairk.PairkAln.from_file('./aln_results.json')
    print(aln_results)

pairk.pairk_alignment_needleman - faster
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This method returns the same results as ``pairk.pairk_alignment``, but it is faster.

The difference is that the ``pairk.pairk_alignment_needleman`` method uses the Needleman-Wunsch algorithm to align the k-mers, while the ``pairk.pairk_alignment`` method uses a scoring matrix to exhaustively score the k-mer matches. It ensures that the alignment is gapless by using an extremely high gap opening and extension penalty (-1000000). This will ensure that the alignment is gapless, unless you use a really unusual scoring matrix with very high scores.

This methods take similar arguments as ``pairk.pairk_alignment``, accept that the ``pairk.pairk_alignment_needleman`` method takes an optional ``aligner`` argument. This allows you to create the aligner before calling the method, which is useful if you want to do multiprocessing, so that you're not creating a new aligner for each process. I've found that if you create a new aligner for each process, the memory usage gets very high, as if the memory isn't being released until the entire script finishes

The ``aligner`` object can be created via the ``pairk.create_aligner`` function. This function takes the name of the scoring matrix as an argument and returns the aligner object. If you don't pass the ``aligner`` argument to the ``pairk.pairk_alignment_needleman`` method, it will create a new aligner using the ``matrix_name`` argument. This is fine if you're not doing multiprocessing. If you are doing multiprocessing, I would suggest creating the aligner before calling the method. Or using 1 aligner for each process. If the ``aligner`` argument is passed, the ``matrix_name`` argument is ignored.

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

Because residue embeddings are used, the inputs are slightly different than the previous methods. 

The inputs are:

* ``full_length_dict_in``: a dictionary of full-length sequences, where the keys are the sequence ids and the values are the sequences. This is used to generate the embeddings.
* ``idr_position_map``: a dictionary where the keys are the full-length sequence ids and the values are the start and end positions of the IDR in the full-length sequence (using python indexing). This is used to slice out the IDR embeddings/sequences from the full-length embeddings/sequences.
* ``query_id``: a query sequence id (the sequence to split into k-mers and align with the homologs). This id should be present in ``idr_position_map`` and ``full_length_dict_in``.
* ``k``
