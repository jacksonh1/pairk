import json
import pandas as pd


def make_empty_kmer_ortho_df(positions, ortholog_ids: list[str]):
    cols = ["reference_kmer"] + ortholog_ids
    df = pd.DataFrame(
        index=positions,
        columns=cols,
    )
    return df


def matrix_dict_2_df(matrix_dict, mat_key):
    df = pd.DataFrame(
        matrix_dict[mat_key]["data"],
        columns=matrix_dict[mat_key]["columns"],
        index=matrix_dict[mat_key]["index"],
    )
    return df


def import_pairwise_matrices(
    filepath,
    matrix_keys: list[str] = [
        "score_dataframe",
        "subseq_dataframe",
        "reciprocal_best_match_dataframe",
        "position_dataframe",
    ],
):
    with open(filepath, "r") as json_file:
        data = json.load(json_file)
    matrices = {}
    for k in matrix_keys:
        if k in data:
            matrices[k] = matrix_dict_2_df(data, k)
    if len(matrices) == 0:
        raise ValueError("No matrices found in the json file")
    return matrices
