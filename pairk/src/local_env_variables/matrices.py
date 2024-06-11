import os
from pathlib import Path

import pandas as pd
from Bio import Align

src_dir = Path(os.path.dirname(__file__)).parent
MATRIX_DIR = src_dir / 'local_data/substitution_matrices/'

# these are compatible with biopython's Align aligner
MATRIX_FILE_DICT = {
    'grantham_similarity': MATRIX_DIR / "grantham_similarity_normx100_aligner_compatible",
    'grantham_similarity_norm': MATRIX_DIR / "grantham_similarity_norm.csv",
    'EDSSMat50': MATRIX_DIR / "EDSSMat50",
    'BLOSUM62': MATRIX_DIR / "BLOSUM62",
}

ALIGNER_COMPATIBLE_MATRICES = ['grantham_similarity', 'EDSSMat50']
BIOPYTHON_MATRICES = Align.substitution_matrices.load()

def convert_matrix_array_2_df(matrix_array: Align.substitution_matrices.Array) -> pd.DataFrame:
    # convert to pandas dataframe
    aas = list(set([i[0] for i in matrix_array.keys()]))
    mat_df = pd.DataFrame(index=aas, columns=aas)
    for k in matrix_array.keys():
        mat_df.loc[k[0], k[1]] = matrix_array[k]
    mat_df = mat_df.astype(float)
    return mat_df


def load_matrix_as_df(matrix_name: str) -> pd.DataFrame:
    # convert to pandas dataframe
    if matrix_name in ALIGNER_COMPATIBLE_MATRICES:
        mat = Align.substitution_matrices.read(MATRIX_FILE_DICT[matrix_name])
        mat_df = convert_matrix_array_2_df(mat)
    elif matrix_name in BIOPYTHON_MATRICES:
        mat = Align.substitution_matrices.load(matrix_name)
        mat_df = convert_matrix_array_2_df(mat)
    else:
        mat_df = pd.read_csv(MATRIX_FILE_DICT[matrix_name], index_col=0)
    return mat_df


def load_matrix_for_aligner(matrix_name: str) -> Align.substitution_matrices.Array:
    if matrix_name in ALIGNER_COMPATIBLE_MATRICES:
        mat = Align.substitution_matrices.read(MATRIX_FILE_DICT[matrix_name])
    elif matrix_name in BIOPYTHON_MATRICES:
        mat = Align.substitution_matrices.load(matrix_name)
    else:
        raise ValueError(f"Matrix {matrix_name} is not compatible with Aligner")
    return mat


def matrix_df_to_dict(matrix_df: pd.DataFrame) -> dict[str, dict[str, float]]:
    '''
    It's much faster to access a dictionary than a pandas dataframe
    '''
    matrix_dict = {}
    for i, row in matrix_df.iterrows():
        if i in matrix_dict:
            raise ValueError(f"Duplicate key {i}")
        matrix_dict[i] = {}
        for j, val in row.items():
            matrix_dict[i][j] = val
    return matrix_dict

    
def print_available_matrices():
    print("aligner-compatible matrices:")
    for k in ALIGNER_COMPATIBLE_MATRICES:
        print(k)
    print("biopython-builtin matrices (aligner compatible):")
    for k in BIOPYTHON_MATRICES:
        print(k)
    print("other matrices:")
    for k in MATRIX_FILE_DICT:
        if k not in ALIGNER_COMPATIBLE_MATRICES and k not in BIOPYTHON_MATRICES:
            print(k)

