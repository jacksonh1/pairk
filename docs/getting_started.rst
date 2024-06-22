Getting Started
===============

Pairk is a Python package designed to quantify the conservation of short motifs in disordered regions, where traditional multiple sequence alignments (MSAs) are difficult. It takes a sequence of interest (the query) and a set of homologous sequences as input and quantifies the conservation of the query sequence without using an MSA. Pairk is designed to be simple to use and easy to install.

See our manuscript for more details on the method: <coming soon>

to install

.. code-block:: bash

    pip install pairk

or for an editable install that you can modify:

.. code-block:: bash

    git clone https://github.com/jacksonh1/pairk.git
    cd pairk
    pip install -e .


Pairk Overview
""""""""""""""

.. image:: ./images/fragpair_cartoon_v2.png
    :align: center
    :width: 800



.. .. image:: ./images/f1-example_MSA_problems.png
..     :align: center
..     :width: 800


The pairk pipeline is composed of two main steps:


1: Pairwise k-mer alignment
"""""""""""""""""""""""""""""""""

For each k-mer in the query IDR, step 1 finds the best scoring length k fragment from each homolog IDR in a gapless and pairkwise manner. The position, sequence, and score of the best scoring fragments are stored in the results.

in informal pseudocode, the algorithm looks like this:

.. code-block:: none

    for each k-mer at each position in the query IDR:
       for each homologous IDR:
           for each k-mer in the homologous IDR:
               score the alignment of the two k-mers (with no gaps)
           store the score, position, and sequence of the best scoring homologous k-mer in DataFrames


2: k-mer conservation
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For each k-mer in the query IDR, step 2 constructs a pseudo-MSA composed of the query k-mer and the best scoring k-mers from each homologous IDR. The conservation of each position in the pseudo-MSA is scored using a column-wise scoring function. The scores are then converted to z-scores to give the conservation relative to other residues in the query IDR.

.. code-block:: none

    for each k-mer in the query IDR:
        from the step 1 results, construct a psuedo-MSA composed of the query k-mer and the best scoring k-mers from each homologous IDR
        for each position in the pseudo-MSA:
            score the conservation of the position using a column-wise scoring function
    Convert all scores for all k-mers pseudo-MSAs to a z-score


quickstart
""""""""""

Here's a quick example to get you started:

.. code-block:: python

    import pairk

    # Load example dataset
    ex1 = pairk.example1

    # Perform k-mer alignment
    aln_results = pairk.pairk_alignment(idr_dict_in=ex1.idr_dict, query_id=ex1.query_id, k=5, matrix_name="EDSSMat50")

    # Calculate conservation
    conservation_results = pairk.calculate_conservation(aln_results)

    # Plot conservation scores
    conservation_results.plot_conservation_mosaic(position=0)

    # save the results
    conservation_results.write_to_file('example1_results.csv')


Check out PairK's main functions and classes:

* step 1, pairwise k-mer alignment functions

    * :func:`pairk.pairk_alignment` - uses a scoring matrix to align the k-mers to each homolog
    * :func:`pairk.pairk_alignment_needleman` - uses a scoring matrix to align the k-mers to each homolog (use pairk.make_aligner to create the aligner object before using this function)
    * :func:`pairk_alignment_embedding_distance` - uses embeddings from ESM2 to align the k-mers to each homolog. (use pairk.ESM_Model to load the model before using this function)
    * :class:`pairk.PairkAln` - pairkwise k-mer alignment results are returned as an instance of this class. See the associated methods for more details on how to interact with the results.

* step 2, k-mer conservation functions

    * :func:`pairk.calculate_conservation` - calculates the conservation of the k-mers from a pairk.PairkAln object
    * :class:`pairk.PairkConservation` - conservation results are returned as an instance of this class. See the associated methods for more details on how to interact with the results.


* utility functions

    * :func:`pairk.print_available_matrices` - prints the available scoring matrices
    * :class:`pairk.FastaImporter` - class to import fasta files and return the sequences in different formats


see the `User Guide <https://pairk.readthedocs.io/en/latest/user_guide.html>`_ page for more details on how to use pairk.

See our `tutorial notebook <https://github.com/jacksonh1/pairk/blob/main/demo/demo.ipynb>`_ for a notebook-based tutorial on how to use pairk.

See the `API documentation <https://pairk.readthedocs.io/en/latest/api.html>`_ for more details on the functions and classes in pairk.

