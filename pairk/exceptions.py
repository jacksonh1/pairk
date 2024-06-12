import pairk.tools.matrices as matrices


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


def check_refid_in_idr_dict(idr_dict_in, reference_id):
    if reference_id not in idr_dict_in:
        raise ValueError(
            f"reference_id '{reference_id}' not found in the idr_dict_in dictionary"
        )
