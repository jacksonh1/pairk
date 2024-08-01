import pairk.backend.tools.matrices as matrices


class InvalidKeywordValue(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def validate_matrix_name(matrix_name):
    available_matrices = list(matrices.AVAILABLE_MATRIX_FILES.keys())
    biopython_matrices = list(matrices.BIOPYTHON_MATRICES)
    if matrix_name not in biopython_matrices:
        if matrix_name not in available_matrices:
            raise InvalidKeywordValue(
                f"matrix '{matrix_name}' not found in allowed matrices. use one of the following:\n    {available_matrices}\nor one of the built-in biopython matrices:\n    {biopython_matrices}"
            )


def check_queryid_in_idr_dict(idr_dict: dict, query_id: str):
    if query_id not in idr_dict:
        raise ValueError(f"query_id '{query_id}' not found in the input dictionary")


def check_sequence_characters_dict(sequence_dict: dict[str, str], matrix_name: str):
    matrix_df = matrices.load_matrix_as_df(matrix_name)
    allowed_characters = list(matrix_df.index)
    for i, seq in sequence_dict.items():
        for aa in seq:
            if aa not in allowed_characters:
                raise ValueError(
                    f"sequence '{i}' contains character '{aa}' not found in the matrix '{matrix_name}'"
                )


def check_sequence_characters_list(sequence_list: list[str], matrix_name: str):
    matrix_df = matrices.load_matrix_as_df(matrix_name)
    allowed_characters = list(matrix_df.index)
    for seq in sequence_list:
        for aa in seq:
            if aa not in allowed_characters:
                raise ValueError(
                    f"sequence contains character '{aa}' not found in the matrix '{matrix_name}'"
                )
