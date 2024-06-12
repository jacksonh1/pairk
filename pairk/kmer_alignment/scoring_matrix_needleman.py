from Bio import Align, Seq
import copy
import pairk.tools.sequence_utils as tools
import pairk.tools.pairwise_tools as pairwise_tools
import pairk.tools.matrices as matrices
import pairk.kmer_alignment.needleman_tools as needleman_tools
import pandas as pd


def needleman_pairwise_kmer_alignment(
    ref_idr: str,
    ortholog_idrs: dict[str, str],
    k: int,
    aligner: Align.PairwiseAligner,
):
    kmers = tools.gen_kmers(ref_idr, k)
    positions = list(range(len(kmers)))

    score_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    subseq_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    pos_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    rbm_df = pairwise_tools.make_empty_kmer_ortho_df(
        positions, list(ortholog_idrs.keys())
    )
    score_df.loc[positions, "reference_kmer"] = kmers
    subseq_df.loc[positions, "reference_kmer"] = kmers
    pos_df.loc[positions, "reference_kmer"] = kmers
    rbm_df.loc[positions, "reference_kmer"] = kmers
    for position, kmer in zip(positions, kmers):
        for ortholog_id, ortholog_idr in ortholog_idrs.items():
            if len(ortholog_idr) < k:
                subseq_df.loc[position, ortholog_id] = "-" * k
                continue
            try:
                best_score, best_subseq, best_pos = (
                    needleman_tools.score_kmer_2_seq_no_gaps_needleman(
                        kmer, ortholog_idr, aligner
                    )
                )
            except ValueError as e:
                subseq_df.loc[position, ortholog_id] = "-" * k
                print(e)
                continue
            except KeyError as e:
                subseq_df.loc[position, ortholog_id] = "-" * k
                print(e)
                continue
            score_df.loc[position, ortholog_id] = best_score
            subseq_df.loc[position, ortholog_id] = best_subseq  # type: ignore
            pos_df.loc[position, ortholog_id] = best_pos
            try:
                (
                    reci_best_score,
                    reci_best_subseq,
                    _,
                ) = needleman_tools.score_kmer_2_seq_no_gaps_needleman(
                    best_subseq, ref_idr, aligner  # type: ignore
                )
            except ValueError as e:
                print(e)
                continue
            except KeyError as e:
                print(e)
                continue
            rbm_df.loc[position, ortholog_id] = reci_best_subseq == kmer
    return score_df, subseq_df, pos_df, rbm_df


def pairk_alignment_needleman(
    idr_dict_in: dict[str, str],
    reference_id: str,
    k: int,
    matrix_name: str = "grantham_similarity_norm",
) -> dict[str, pd.DataFrame]:
    """run pairwise k-mer alignment method using the needleman-wunsch algorithm as implemented in Biopython

    Parameters
    ----------
    idr_dict_in : dict[str, str]
        input sequences in dictionary format with the key being the sequence id and the value being the sequence as a string
    reference_id : str
        the id of the reference sequence within the `idr_dict_in` dictionary
    k : int
        the length of the k-mers to use for the alignment
    matrix_name : str, optional
        The name of the scoring matrix to use in the algorithm. The available matrices can be viewed with the function `print_available_matrices()` in pairk.tools.matrices, by default "grantham_similarity_norm"

    Returns
    -------
    dict[str, pd.DataFrame]
        results of the pairwise alignment in dictionary format, where the keys are the names of the dataframes and the values are the dataframes
    """
    idr_dict = copy.deepcopy(idr_dict_in)
    ref_idr = idr_dict.pop(reference_id)
    matrix = matrices.load_matrix_for_aligner(matrix_name)
    aligner = needleman_tools.get_aligner(matrix)
    score_df, subseq_df, pos_df, rbm_df = needleman_pairwise_kmer_alignment(
        ref_idr,
        idr_dict,
        k,
        aligner=aligner,
    )
    output_dict = {
        "score_dataframe": score_df,
        "subseq_dataframe": subseq_df,
        "position_dataframe": pos_df,
        "reciprocal_best_match_dataframe": rbm_df,
    }
    return output_dict
