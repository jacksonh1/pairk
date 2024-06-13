from Bio import Align, Seq
import copy
import pairk.backend.tools.sequence_utils as tools
import pairk.backend.tools.pairwise_tools as pairwise_tools
import pairk.backend.tools.matrices as matrices
import pairk.backend.kmer_alignment.needleman_tools as needleman_tools
import pairk.backend.exceptions as _exceptions
import pandas as pd


def run_pairwise_kmer_alignment_needleman(
    query_idr: str,
    ortholog_idrs: dict[str, str],
    k: int,
    aligner: Align.PairwiseAligner,
):
    kmers = tools.gen_kmers(query_idr, k)
    positions = list(range(len(kmers)))

    score_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    orthokmer_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    pos_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    rbm_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    score_df.loc[positions, "query_kmer"] = kmers
    orthokmer_df.loc[positions, "query_kmer"] = kmers
    pos_df.loc[positions, "query_kmer"] = kmers
    rbm_df.loc[positions, "query_kmer"] = kmers
    for position, kmer in zip(positions, kmers):
        for ortholog_id, ortholog_idr in ortholog_idrs.items():
            if len(ortholog_idr) < k:
                orthokmer_df.loc[position, ortholog_id] = "-" * k
                continue
            try:
                best_score, best_subseq, best_pos = (
                    needleman_tools.score_kmer_2_seq_needleman(
                        kmer, ortholog_idr, aligner
                    )
                )
            except ValueError as e:
                orthokmer_df.loc[position, ortholog_id] = "-" * k
                print(e)
                continue
            except KeyError as e:
                orthokmer_df.loc[position, ortholog_id] = "-" * k
                print(e)
                continue
            score_df.loc[position, ortholog_id] = best_score
            orthokmer_df.loc[position, ortholog_id] = best_subseq  # type: ignore
            pos_df.loc[position, ortholog_id] = best_pos
            try:
                (
                    reci_best_score,
                    reci_best_subseq,
                    _,
                ) = needleman_tools.score_kmer_2_seq_needleman(
                    best_subseq, query_idr, aligner  # type: ignore
                )
            except ValueError as e:
                print(e)
                continue
            except KeyError as e:
                print(e)
                continue
            rbm_df.loc[position, ortholog_id] = reci_best_subseq == kmer
    return score_df, orthokmer_df, pos_df, rbm_df


def pairk_alignment_needleman(
    idr_dict_in: dict[str, str],
    query_id: str,
    k: int,
    aligner: Align.PairwiseAligner | None = None,
    matrix_name: str = "EDSSMat50",
) -> pairwise_tools.PairkAln:
    """
    run pairwise k-mer alignment method using the needleman-wunsch algorithm as implemented in Biopython.

    Parameters
    ----------
    idr_dict_in : dict[str, str]
        input sequences in dictionary format with the key being the sequence id and the value being the sequence as a string
    query_id : str
        the id of the query sequence within the `idr_dict_in` dictionary
    k : int
        the length of the k-mers to use for the alignment
    aligner : Align.PairwiseAligner | None, optional
        The Biopython pairwise aligner object to use in the pairwise gapless alignments, by default None. If None, then an aligner object will be created using the scoring matrix specified in `matrix_name`. If an aligner object is provided, it will take precedence over the `matrix_name` parameter, i.e. the `matrix_name` parameter will be ignored.
    matrix_name : str, optional
        The name of the scoring matrix to use in the algorithm, by default "EDSSMat50". The available matrices can be viewed with the function `print_available_matrices()` in `pairk.backend.tools.matrices`. If an aligner object is provided, this parameter **will be ignored**.
    rbm : bool, optional
        whether to calculate the reciprocal best match (RBM) for each k-mer, by default False.
        It true, the RBM matrix will be calculated by aligning the best scoring
        ortholog k-mer to the query sequence and checking if the best
        scoring query k-mer is the same as the original query k-mer.
        This is currently unused in the pairk package. The RBM matrix will be
        included in the output `PairkAln` object, it is a boolean dataframe
        indicating whether the query k-mer ortholog k-mer match is reciprocal.

    Returns
    -------
    pairwise_tools.PairkAln
        an object containing the alignment results. See the `pairk.PairkAln` class for more information.
    """
    _exceptions.check_queryid_in_idr_dict(idr_dict_in, query_id)
    if aligner is None:
        _exceptions.validate_matrix_name(matrix_name)
        aligner = needleman_tools.make_aligner(matrix_name)
    idr_dict = copy.deepcopy(idr_dict_in)
    query_idr = idr_dict.pop(query_id)
    score_df, orthokmer_df, pos_df, rbm_df = run_pairwise_kmer_alignment_needleman(
        query_idr,
        idr_dict,
        k,
        aligner=aligner,
    )
    return pairwise_tools.PairkAln(
        orthokmer_df=orthokmer_df,
        pos_df=pos_df,
        score_df=score_df,
        rbm_df=rbm_df,
    )
