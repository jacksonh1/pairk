=================
API Documentation
=================


.. contents:: Table of Contents
    :depth: 5
    :local:




.. autosummary::
   :toctree: autosummary

****************************************
PairK convencience functions  
****************************************

.. autoclass:: pairk.FastaImporter
   :members:
   :undoc-members:


.. autofunction:: pairk.print_available_matrices


****************************************
pairwise k-mer alignment
****************************************

.. autofunction:: pairk.pairk_alignment  

.. autofunction:: pairk.pairk_alignment_needleman  

.. autofunction:: pairk.make_aligner

.. autofunction:: pairk.pairk_alignment_embedding_distance

.. autoclass:: pairk.ESM_Model
   :members:


.. autoclass:: pairk.PairkAln
   :members:
..    :undoc-members:


..    :undoc-members:


****************************************
k-mer conservation
****************************************


.. autofunction:: pairk.calculate_conservation

.. autofunction:: pairk.calculate_conservation_arrays

.. autoclass:: pairk.PairkConservation
   :members:



.. .. automodule:: pairk
..    :members:


*****************
utility functions
*****************

.. automodule:: pairk.utilities
   :members:


   .. :show-inheritance:


**********************
single k-mer functions
**********************

.. automodule:: pairk.single_kmer
   :members:



   .. :show-inheritance:
