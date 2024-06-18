Getting Started
===============

.. image:: ./images/fragpair_cartoon_v2.png
    :align: center
    :width: 800


.. code-block:: bash

    pip install pairk

or for an editable install that you can modify:

.. code-block:: bash

    git clone https://github.com/jacksonh1/pairk.git
    cd pairk
    pip install -e .



See our `tutorial notebook <https://github.com/jacksonh1/pairk/blob/main/demo/demo.ipynb>`_ for a notebook-based tutorial.

Pairk is a simple tool to align k-mers from one IDR to all homologous IDRs in a set of homologous sequences. The alignment is gapless and performed in pairwise manner. The most simple version of the alorithm works as follows (in informal pseudocode):


1: Pairwise k-mer alignment
"""""""""""""""""""""""""""""""""

.. code-block:: none

    for each k-mer at each position in the query IDR:
       for each homologous IDR:
           for each k-mer in the homologous IDR:
               score the alignment of the two k-mers (with no gaps)
           store the score, position, and sequence of the best scoring homologous k-mer in DataFrames


For each k-mer in the query IDR, step 1 finds the best scoring length k fragment from each homolog IDR. The position, sequence, and score of the best scoring fragments are accessible from the results.


2: score the relative conservation of each k-mer in the query IDR
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: none

    for each k-mer in the query IDR:
        from the step 1 results, construct a psuedo-MSA composed of the query k-mer and the best scoring k-mers from each homologous IDR
        for each position in the pseudo-MSA:
            score the conservation of the position using a column-wise scoring function
    Convert all scores for all k-mers pseudo-MSAs to a z-score





