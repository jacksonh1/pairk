from Bio import Align, Seq
import copy
import pairk.tools.sequence_utils as tools
import pairk.tools.pairwise_tools as pairwise_tools
import pairk.tools.matrices as matrices
import pairk.kmer_alignment.needleman_tools as needleman_tools
import pandas as pd
import pairk.exceptions as _exceptions


def run_pairwise_kmer_alignment_needleman(
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
                    needleman_tools.score_kmer_2_seq_needleman(
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
                ) = needleman_tools.score_kmer_2_seq_needleman(
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
    matrix_name: str = "EDSSMat50",
) -> dict[str, pd.DataFrame]:
    """
    run pairwise k-mer alignment method using the needleman-wunsch algorithm as implemented in Biopython.

    Parameters
    ----------
    idr_dict_in : dict[str, str]
        input sequences in dictionary format with the key being the sequence id and the value being the sequence as a string
    reference_id : str
        the id of the reference sequence within the `idr_dict_in` dictionary
    k : int
        the length of the k-mers to use for the alignment
    matrix_name : str, optional
        The name of the scoring matrix to use in the algorithm, by default "EDSSMat50". The available matrices can be viewed with the function `print_available_matrices()` in `pairk.tools.matrices`.

    Returns
    -------
    dict[str, pd.DataFrame]
        results of the pairwise alignment in dictionary format, where the keys are the names of the dataframes and the values are the dataframes. All dataframes have the same structure. One column is the reference k-mer sequence ('reference_kmer'). The other columns are named as the ortholog sequence ids. The dataframe indexes are the reference k-mer start position in the reference sequence. The returned dataframes are:\n
        - 'score_dataframe': the alignment scores for each k-mer in the reference sequence against the corresponding best matching ortholog k-mer.
        - 'subseq_dataframe': the best scoring k-mer from each ortholog for each reference k-mer.
        - 'position_dataframe': the start position of the best scoring k-mer from each ortholog for each reference k-mer.
        - 'reciprocal_best_match_dataframe': a boolean dataframe indicating whether the reference k-mer is the reciprocal best scoring k-mer to the ortholog k-mer.

    """
    _exceptions.validate_matrix_name(matrix_name)
    _exceptions.check_refid_in_idr_dict(idr_dict_in, reference_id)
    idr_dict = copy.deepcopy(idr_dict_in)
    ref_idr = idr_dict.pop(reference_id)
    matrix = matrices.load_matrix_for_aligner(matrix_name)
    aligner = needleman_tools.get_aligner(matrix)
    score_df, subseq_df, pos_df, rbm_df = run_pairwise_kmer_alignment_needleman(
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
