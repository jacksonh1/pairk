import os
from pathlib import Path

import pandas as pd
from Bio import Align
from importlib_resources import files

# MATRIX_DIR = Path(str(files('pairk.matrices')))
# AVAILABLE_MATRIX_FILES = {i.stem: i for i in MATRIX_DIR.glob("*") if not i.name.endswith(".py") and not i.name.startswith("__")}
AVAILABLE_MATRIX_FILES = {i.stem: i for i in files("pairk.data.matrices").iterdir()}  # type: ignore

BIOPYTHON_MATRICES = Align.substitution_matrices.load()  # type: ignore
AA_ORDER = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
    "B",
    "Z",
    "X",
    "*",
]


def convert_matrix_array_2_df(
    matrix_array: Align.substitution_matrices.Array,  # type: ignore
) -> pd.DataFrame:
    # convert to pandas dataframe
    aas = list(set([i[0] for i in matrix_array.keys()]))
    # sort aas to match the order in AA_ORDER
    aas = sorted(aas, key=lambda x: AA_ORDER.index(x))
    mat_df = pd.DataFrame(index=aas, columns=aas)
    for k in matrix_array.keys():
        mat_df.loc[k[0], k[1]] = matrix_array[k]
    mat_df = mat_df.astype(float)
    return mat_df


def load_matrix_as_df(matrix_name: str) -> pd.DataFrame:
    # convert to pandas dataframe
    if matrix_name in BIOPYTHON_MATRICES:
        mat = Align.substitution_matrices.load(matrix_name)  # type: ignore
        return convert_matrix_array_2_df(mat)  # type: ignore
    try:
        mat = Align.substitution_matrices.read(AVAILABLE_MATRIX_FILES[matrix_name])  # type: ignore
        return convert_matrix_array_2_df(mat)
    except AssertionError:
        return pd.read_csv(AVAILABLE_MATRIX_FILES[matrix_name], index_col=0)  # type: ignore


def load_matrix_for_aligner(matrix_name: str) -> Align.substitution_matrices.Array:  # type: ignore
    # convert to pandas dataframe
    if matrix_name in BIOPYTHON_MATRICES:
        return Align.substitution_matrices.load(matrix_name)  # type: ignore
    try:
        return Align.substitution_matrices.read(AVAILABLE_MATRIX_FILES[matrix_name])  # type: ignore
    except AssertionError as e:
        raise ValueError(
            f"matrix {matrix_name}: ({AVAILABLE_MATRIX_FILES[matrix_name]}) is probably not compatible with Bio.Align Aligner"
        ) from e


def matrix_df_to_dict(matrix_df: pd.DataFrame) -> dict[str, dict[str, float]]:
    """
    It's much faster to access a dictionary than a pandas dataframe
    """
    matrix_dict = {}
    for i, row in matrix_df.iterrows():
        if i in matrix_dict:
            raise ValueError(f"Duplicate key {i}")
        matrix_dict[i] = {}
        for j, val in row.items():
            matrix_dict[i][j] = val
    return matrix_dict


def print_available_matrices():
    """print available matrices that can be used in the pairk aligner"""
    print("biopython-builtin matrices (aligner compatible):")
    for k in BIOPYTHON_MATRICES:
        print(k)
    print("\nother matrices:")
    for k, v in AVAILABLE_MATRIX_FILES.items():
        # print(f"{k}\n - {v}")
        print(f"{k}")
