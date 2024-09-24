===============================
References and acknowledgements
===============================

special thanks to Foster Birnbaum for help accelerating the embedding distance functions

Python modules utilized by PairK are documented in the requirements

**********************
Additional references:
**********************

* ESM2 (the model used to generate the embeddings): 

    * Z. Lin, H. Akin, R. Rao, B. Hie, Z. Zhu, W. Lu, N. Smetanin, R. Verkuil, O. Kabeli, Y. Shmueli, A. Dos Santos Costa, M. Fazel-Zarandi, T. Sercu, S. Candido, A. Rives, Evolutionary-scale prediction of atomic-level protein structure with a language model. Science 379, 1123–1130 (2023).

* Some of the ESM model sequence encoding functions are adapted from the kibby tool: https://github.com/esbgkannan/kibby

    * W. Yeung, Z. Zhou, S. Li, N. Kannan, Alignment-free estimation of sequence conservation for identifying functional sites using protein sequence embeddings. Brief Bioinform 24 (2023).

* Pairk's built-in conservation scoring functions are adapted from code released with this study: 

    * J. A. Capra, M. Singh, Predicting functionally important residues from sequence conservation. Bioinformatics 23, 1875–1882 (2007)

* Pairk's built-in scoring matrix "EDSSMat50" is from this study: 

    * R. Trivedi, H. A. Nagarajaram, Amino acid substitution scoring matrices specific to intrinsically disordered regions in proteins. Sci Rep 9, 16380 (2019)

* Pairk's built-in "grantham" matrices (including "grantham", "grantham_similarity_norm", and "grantham_similarity_normx100_aligner_compatible") are from or derived from the distance matrix in this study: 

    * R. Grantham, Amino acid difference formula to help explain protein evolution. Science 185, 862–864 (1974).

* blosum62 matrix is from biopython:

    * P. J. A. Cock, T. Antao, J. T. Chang, B. A. Chapman, C. J. Cox, A. Dalke, I. Friedberg, T. Hamelryck, F. Kauff, B. Wilczynski, M. J. L. de Hoon, Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25, 1422–1423 (2009).
    * S. Henikoff, J. G. Henikoff, Amino acid substitution matrices from protein blocks. Proc Natl Acad Sci U S A 89, 10915–10919 (1992).

