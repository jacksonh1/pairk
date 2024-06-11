from .aln_asym_sum_of_pairs import main as score_aln_asym_sum_of_pairs
from .aln_property_entropy import main as score_aln_property_entropy
from .aln_shannon_entropy import main as score_aln_shannon_entropy
from .fragment_pairwise_gapless import main_save_file as score_fragment_pairwise_gapless
from .frag_pairwise_gapless_embedding import main as score_frag_pairwise_gapless_embedding
from .pairwise_matrix_to_kmer_scores import matrix_json_2_pairwise_scores
from .pairwise_matrix_to_kmer_scores_rbm_penalty import matrix_json_2_pairwise_scores as matrix_json_2_pairwise_scores_rbm_penalty
from functools import partial
from local_conservation_scores.tools import capra_singh_2007_scores


class ConservationScoreMethods:
    def __init__(self):
        self.aln_asym_sum_of_pairs = score_aln_asym_sum_of_pairs
        self.aln_property_entropy = score_aln_property_entropy
        self.aln_shannon_entropy = score_aln_shannon_entropy
        # self.fragment_pairwise_gapless = score_fragment_pairwise_gapless
    
    def __getitem__(self, key):
        return getattr(self, key)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
    
    def __delitem__(self, key):
        delattr(self, key)

# CONSERVATION_SCORE_METHODS = conservation_score_methods()

class PairwiseMatrixMethods:
    def __init__(self):
        self.fragment_pairwise_gapless = score_fragment_pairwise_gapless
        self.frag_pairwise_gapless_embedding = score_frag_pairwise_gapless_embedding
    
    def __getitem__(self, key):
        return getattr(self, key)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
    
    def __delitem__(self, key):
        delattr(self, key)


class PairwiseMatrixKmerScoreMethods:
    def __init__(self):
        self.pairwise_matrix_to_kmer_scores = matrix_json_2_pairwise_scores
        self.pairwise_matrix_to_kmer_scores_rbm_penalty = matrix_json_2_pairwise_scores_rbm_penalty

    def __getitem__(self, key):
        return getattr(self, key)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
    
    def __delitem__(self, key):
        delattr(self, key)


class ColumnwiseScoreMethods:
    def __init__(self):
        self.shannon_entropy = capra_singh_2007_scores.shannon_entropy
        self.property_entropy = capra_singh_2007_scores.property_entropy
        self.vn_entropy = capra_singh_2007_scores.vn_entropy

    def __getitem__(self, key):
        return getattr(self, key)
    
    def __setitem__(self, key, value):
        setattr(self, key, value)
    
    def __delitem__(self, key):
        delattr(self, key)