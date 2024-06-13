Getting Started
===============

.. code-block:: bash

    pip install pairk

or for an editable install that you can modify:

.. code-block:: bash

    git clone https://github.com/jacksonh1/pairk.git
    cd pairk
    pip install -e .


Pairk is a simple tool to align k-mers from one IDR to all homologous IDRs in a set of homologous sequences. The alignment is gapless and performed in pairwise manner. The most simple version of the alorithm works as follows:

1: Pairwise k-mer alignment
"""""""""""""""""""""""""""""""""

.. code-block:: none

    for each k-mer in the query IDR:
       for each homologous IDR:
           for each k-mer in the homologous IDR:
               score the alignment of the two k-mers (with no gaps)
           store the score, position, and sequence of the best scoring homologous k-mer


2: score the relative conservation of each k-mer in the query IDR
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

.. code-block:: none

    for each k-mer in the query IDR:
        construct a psuedo-MSA composed of the query k-mer and the best scoring k-mers from each homologous IDR
        for each position in the psuedo-MSA:
            score the conservation of the position using a column-wise scoring function
    Convert all scores for all k-mers to a z-score





