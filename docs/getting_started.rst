.. _getting_started:

===============
Getting Started
===============

Pairk is a Python package designed to quantify the conservation of short motifs in disordered regions, where traditional multiple sequence alignments (MSAs) are difficult. It takes a sequence of interest (the query) and a set of homologous sequences as input and quantifies the conservation of the query sequence without using an MSA. Pairk is designed to be simple to use and easy to install.

See our manuscript for more details on the method: `PairK: Pairwise k-mer alignment for quantifying protein motif conservation in disordered regions <https://www.biorxiv.org/content/10.1101/2024.07.23.604860v1>`_

to install

.. code-block:: bash

    pip install pairk

or for an editable install that you can modify:

.. code-block:: bash

    git clone https://github.com/jacksonh1/pairk.git
    cd pairk
    pip install -e .


**************
Pairk Overview
**************

The pairk pipeline:

.. image:: ./images/fragpair_cartoon_v2.png
    :align: center
    :width: 800


SLiM conservation from an MSA vs. pairk:

.. image:: ./images/f1-example_MSA_problems.png
    :align: center
    :width: 400


.. raw:: html

   <br>
   <br>
   <br>


**The pairk pipeline is composed of two main steps:**


1: Pairwise k-mer alignment
=================================

For each k-mer in the query IDR, step 1 finds the best scoring length k fragment from each homolog IDR in a gapless and pairkwise manner. The position, sequence, and score of the best scoring fragments are stored in the results.

in informal pseudocode, the algorithm looks like this:

.. code-block:: none

    for each k-mer at each position in the query IDR:
        for each homologous IDR:
            for each k-mer at each position in the homologous IDR:
                score the alignment of the homolog k-mer - query k-mer match (with no gaps)
            store the score, position, and sequence of the best scoring homologous k-mer
        construct a "pseudo-MSA" composed of the query k-mer and the best scoring k-mers from each homologous IDR


2: k-mer conservation
=================================

For each k-mer "pseudo-MSA" from step 1, step 2 calculates the conservation of each position in the pseudo-MSA using a column-wise scoring function. The scores are then converted to z-scores to give the conservation relative to other residues in the query IDR. The z-score background distribution is all of the column-wise scores from all of the pseudo-MSAs.

.. code-block:: none

    for each k-mer in the query IDR:
        for each position in the pseudo-MSA:
            score the conservation of the position using a column-wise scoring function
    Convert all scores for all k-mer pseudo-MSAs to z-scores


**********
quickstart
**********

Here's a quick example to get you started:

.. code-block:: python

    import pairk

    # Load example dataset
    ex1 = pairk.example1

    # Perform k-mer alignment
    aln_results = pairk.pairk_alignment(idr_dict=ex1.idr_dict, query_id=ex1.query_id, k=5, matrix_name="EDSSMat50")


The resulting pseudo-MSAs are stored in the ``aln_results`` object (an instance of the :class:`pairk.PairkAln` class). You can access the results from this object's attributes directly (i.e. the :attr:`pairk.PairkAln.orthokmer_matrix`, :attr:`pairk.PairkAln.position_matrix`, and :attr:`pairk.PairkAln.score_matrix` DataFrames). 

Let's say that we are interested in the k-mer "LPPPP" which starts at position 75 in the query sequence. We can access the "pseudo-MSA" for this k-mer from the :attr:`pairk.PairkConservation.orthokmer_matrix` DataFrame or using the :attr:`pairk.PairkConservation.get_pseudo_alignment()` method:

.. code-block:: python

    In [4]: print(aln_results.get_pseudo_alignment(75))
    ['LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'PPMPP', 'LPPPP', 'LPDRP', 'APSPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'LPPPP', 'IPPPP']

    In [5]: print(aln_results.orthokmer_matrix.loc[75])
    query_kmer          LPPPP
    9793_0:005123       LPPPP
    1706337_0:000fc7    LPPPP
    51337_0:001b5a      LPPPP
    9568_0:004ae1       PPMPP
    43346_0:004190      LPPPP
    885580_0:00488c     LPDRP
    10181_0:00305d      APSPP
    1415580_0:000900    LPPPP
    61221_0:00105a      LPPPP
    7897_0:0033c5       LPPPP
    8407_0:002bff       LPPPP
    173247_0:004550     LPPPP
    30732_0:0046dd      LPPPP
    241271_0:0048e4     LPPPP
    8103_0:0045e4       LPPPP
    56723_0:00152f      LPPPP
    210632_0:004c0c     LPPPP
    31033_0:00264e      LPPPP
    63155_0:004c86      LPPPP
    7994_0:004d71       LPPPP
    109280_0:00369f     LPPPP
    150288_0:004e5a     IPPPP
    Name: 75, dtype: object

To calculate pairk conservation, an instance of the :class:`pairk.PairkAln` object is used as input to the :func:`pairk.calculate_conservation` function.

.. code-block:: python

    # Calculate conservation
    conservation_results = pairk.calculate_conservation(aln_results)

The conservation results are stored in ``conservation_results`` (an instance of the :class:`pairk.PairkConservation` class). You can access the results from this object's attributes (e.g. the :attr:`pairk.PairkConservation.orthokmer_arr`, :attr:`pairk.PairkConservation.score_arr`, and :attr:`pairk.PairkConservation.z_score_arr` arrays)

The :class:`pairk.PairkConservation` object also contains some plotting functions. for example, :func:`pairk.PairkConservation.plot_conservation_mosaic`:

.. code-block:: python

    # Plot conservation scores
    conservation_results.plot_conservation_mosaic(position=75) # start position of the query k-mer of interest


.. image:: ./images/mosaic_plot_annotated.png
    :align: center
    :width: 800


|

The above example output of the :func:`pairk.PairkConservation.plot_conservation_mosaic` is annotated with explanation of each element of the plot

The :class:`pairk.PairkConservation` also has methods to write the results to a file or read the results from a file:
 
.. code-block:: python
    
    # save the results
    conservation_results.write_to_file('example1_results.npz')

    # read the results
    conservation_results = pairk.PairkConservation.read_results_from_file('example1_results.npz')



*********************************************
PairK's main functions and classes
*********************************************

* step 1, pairwise k-mer alignment functions

    * :func:`pairk.pairk_alignment` - uses a scoring matrix to align the k-mers to each homolog
    * :func:`pairk.pairk_alignment_needleman` - uses a scoring matrix to align the k-mers to each homolog (use pairk.make_aligner to create the aligner object before using this function)
    * :func:`pairk.pairk_alignment_embedding_distance` - uses embeddings from ESM2 (or user provided embeddings) to align the k-mers to each homolog. To use ESM2, use :func:`pairk.ESM_Model` to load the model before using this function (provided as the ``mod`` argument)
    * :class:`pairk.PairkAln` - pairkwise k-mer alignment results are returned as an instance of this class. See the associated methods for more details on how to interact with the results.

* step 2, k-mer conservation functions

    * :func:`pairk.calculate_conservation` - calculates the conservation of the k-mers from a pairk.PairkAln object
    * :class:`pairk.PairkConservation` - conservation results are returned as an instance of this class. See the associated methods for more details on how to interact with the results.


* utility functions

    * :func:`pairk.print_available_matrices` - prints the available scoring matrices
    * :class:`pairk.FastaImporter` - class to import fasta files and return the sequences in different formats


.. see the `User Guide <https://pairk.readthedocs.io/en/latest/user_guide.html>`_ page for more details on how to use pairk.

See our `tutorial notebook <https://github.com/jacksonh1/pairk/blob/main/demo/pairk_tutorial.ipynb>`_ for a notebook-based tutorial on how to use pairk.

See the `API documentation <https://pairk.readthedocs.io/en/latest/api.html>`_ for more details on the functions and classes in pairk.


