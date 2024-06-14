"""Pairk for single kmer analysis. This is available for custom analysis. 
Relative conservation scores (z-scores) not available for single kmers but it's 
considerably faster to align 1 kmer than all kmers in a sequence. We make this 
available for custom analysis."""

from pairk.backend.kmer_alignment.scoring_matrix import pairk_alignment_single_kmer

# Define the __all__ variable
__all__ = [
    'pairk_alignment_single_kmer'
]

