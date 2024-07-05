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
    score_df.loc[positions, "query_kmer"] = kmers
    orthokmer_df.loc[positions, "query_kmer"] = kmers
    pos_df.loc[positions, "query_kmer"] = kmers
    for position, kmer in zip(positions, kmers):
        for ortholog_id, ortholog_idr in ortholog_idrs.items():
            if len(ortholog_idr) < k:
                orthokmer_df.loc[position, ortholog_id] = "-" * k
                continue
            best_score, best_subseq, best_pos = (
                needleman_tools.score_kmer_2_seq_needleman(kmer, ortholog_idr, aligner)
            )
            score_df.loc[position, ortholog_id] = best_score
            orthokmer_df.loc[position, ortholog_id] = best_subseq  # type: ignore
            pos_df.loc[position, ortholog_id] = best_pos
            (
                reci_best_score,
                reci_best_subseq,
                _,
            ) = needleman_tools.score_kmer_2_seq_needleman(
                best_subseq, query_idr, aligner  # type: ignore
            )
    return score_df, orthokmer_df, pos_df


def pairk_alignment_needleman(
    idr_dict: dict[str, str],
    query_id: str,
    k: int,
    aligner: Align.PairwiseAligner | None = None,
    matrix_name: str = "EDSSMat50",
) -> pairwise_tools.PairkAln:
    """
    run pairwise k-mer alignment method using the needleman-wunsch algorithm as
    implemented in Biopython. Each query k-mer is scored against each ortholog
    k-mer to find the best matching ortholog k-mer in each ortholog. If a ortholog
    IDR is shorter than the k-mer, a string of "-" characters  ("-"\\*k) is
    assigned as the best matching ortholog k-mer for that ortholog

    **Note**: if there are multiple top-scoring matches, only one is returned.

    Parameters
    ----------
    idr_dict : dict[str, str]
        input sequences in dictionary format with the key being the sequence id and
        the value being the sequence as a string
    query_id : str
        the id of the query sequence within the `idr_dict` dictionary
    k : int
        the length of the k-mers to use for the alignment
    aligner : Align.PairwiseAligner | None, optional
        The Biopython pairwise aligner object to use in the pairwise gapless alignments,
        by default None. If None, then an aligner object will be created using the scoring
        matrix specified in `matrix_name`. If an aligner object is provided, it will
        take precedence over the `matrix_name` parameter, i.e. the `matrix_name` parameter
        will be ignored.
    matrix_name : str, optional
        The name of the scoring matrix to use in the algorithm, by default "EDSSMat50".
        The available matrices can be viewed with the function `print_available_matrices()`
        in `pairk.backend.tools.matrices`. If an aligner object is provided,
        this parameter **will be ignored**.

    Returns
    -------
    pairwise_tools.PairkAln
        an object containing the alignment results. See the `pairk.PairkAln` class for more information.
    """
    _exceptions.check_queryid_in_idr_dict(idr_dict, query_id)
    idr_str_dict = copy.deepcopy(idr_dict)
    if aligner is None:
        _exceptions.validate_matrix_name(matrix_name)
        aligner = needleman_tools.make_aligner(matrix_name)
    query_idr = idr_str_dict.pop(query_id)
    score_df, orthokmer_df, pos_df = run_pairwise_kmer_alignment_needleman(
        query_idr,
        idr_str_dict,
        k,
        aligner=aligner,
    )
    return pairwise_tools.PairkAln(
        orthokmer_df=orthokmer_df,
        pos_df=pos_df,
        score_df=score_df,
    )
