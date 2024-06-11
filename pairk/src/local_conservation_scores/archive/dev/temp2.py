# %%
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
from local_conservation_scores.tools import capra_singh_2007_scores, general
from local_env_variables import matrices
from local_seqtools import alignment_tools as aln_tools
from local_seqtools import general_utils as tools
from local_seqtools import jch_alignment as jch_aln

# ==============================================================================
# //
# ==============================================================================

# %%


# Assuming df1 and df2 are your DataFrames representing matrices
data1 = {"col1": [1, 2, 3], "col2": [4, 5, 6], "col3": [7, 8, 9]}
data2 = {"colA": [10, 11, 12], "colB": [13, 14, 15], "colC": [16, 17, 18]}

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)

# Create a dictionary to store multiple matrices
matrices_dict = {"matrix1": df1, "matrix2": df2}

# Save the dictionary to a JSON file
with open("matrices.json", "w") as json_file:
    json_file.write(pd.Series(matrices_dict).to_json(orient="index", indent=4))


# %%

import json

import pandas as pd

data1 = {"col1": [1, 2, 3], "col2": [4, 5, 6], "col3": [7, 8, 9]}
data2 = {"colA": [10, 11, 12], "colB": [13, 14, 15], "colC": [16, 17, 18]}

df1 = pd.DataFrame(data1)
df2 = pd.DataFrame(data2)

matrices_dict = {
    "matrix1": df1.to_dict(orient="split"),
    "matrix2": df2.to_dict(orient="split"),
    "additional_info": {"key1": "value1", "key2": "value2"},
}
with open("matrices_with_info.json", "w") as json_file:
    json.dump(matrices_dict, json_file, indent=4)

# ==============================================================================
# ==============================================================================

with open("matrices_with_info.json", "r") as json_file:
    data = json.load(json_file)

df1_reconstructed = pd.DataFrame(
    data["matrix1"]["data"],
    columns=data["matrix1"]["columns"],
    index=data["matrix1"]["index"],
)
df2_reconstructed = pd.DataFrame(
    data["matrix2"]["data"],
    columns=data["matrix2"]["columns"],
    index=data["matrix2"]["index"],
)
additional_info = data["additional_info"]

# %%