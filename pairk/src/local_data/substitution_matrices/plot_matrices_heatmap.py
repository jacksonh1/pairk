import json
import os
import re
import sys
from pathlib import Path
from Bio import AlignIO, Seq, SeqIO, Align
import numpy as np
import pandas as pd

# %load_ext autoreload
# %autoreload 2
import matplotlib.pyplot as plt
plt.style.use('custom_standard')
plt.style.use('custom_small')
import seaborn as sns
# pd.options.plotting.backend = "plotly"
# %%

def plot_matrix(mat_df, title=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 7))
    if title is not None:
        ax.set_title(title, fontsize = 20)
    sns.heatmap(
        mat_df,
        annot=True,
        # cmap="coolwarm",
        cmap="vlag",
        fmt=".1g",
        # fmt=".2f",
        linewidths=0.5,
        ax=ax, 
        annot_kws={"size": 8, 'color': 'black', 'weight': 'bold'},
        cbar=False,
        linecolor='black',
        clip_on=False,
        # vmin=0,
        # vmax=1,
    )
    return ax

def plot_and_save_matrix(filename, output_dir, title=None):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    filename = Path(filename)
    mat_df = pd.read_csv(filename, index_col=0)
    if title is None:
        title = filename.stem
    plot_matrix(mat_df, title=title)
    plt.savefig(output_dir / filename.with_suffix('.png'), dpi=300, bbox_inches='tight')
    plt.close()

output_dir = Path('./plots')
output_dir.mkdir(exist_ok=True)
for filename in Path('./').glob('*.csv'):
    plot_and_save_matrix(filename, output_dir)





