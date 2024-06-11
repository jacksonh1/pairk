# %%
import json
from pathlib import Path

import matplotlib.pyplot as plt
from Bio import AlignIO, SeqIO

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_conservation_scores.tools.score_plots as score_plots

# %%


# %%


def multi_level_plots(json_file, score_key, score_type='zscore', num_bg_scores_cutoff=20):
    og = group_tools.ConserGene(json_file)
    og.load_aln_scores(score_key, num_bg_scores_cutoff=num_bg_scores_cutoff)
    output_folder = Path(og.info_dict["analysis_folder"])

    # multi_plot_filename_1 = score_plots.ogscreen_plot_from_og_folder_class(
    #     og,
    #     big_plot_output_folder=output_folder,
    #     score_name=f"{score_key}-stripped",
    #     bar_ylim=[-0.5, 2.5],
    #     strip_gaps=True,
    # )
    # og.add_item_to_json(
    #     f"multilevel_plot_file-gapless-{score_key}",
    #     str(multi_plot_filename_1),
    #     save_json=True,
    # )

    multi_plot_filename_2 = score_plots.ogscreen_plot_from_og_folder_class(
        og,
        big_plot_output_folder=output_folder,
        score_name=f"{score_key}",
        score_type=score_type,
        bar_ylim=[-0.5, 2.5],
        strip_gaps=False,
    )
    og.add_item_to_json(
        f"multilevel_plot_file-{score_key}", str(multi_plot_filename_2), save_json=True
    )
    plt.close("all")


# lvlo.load_scores('property_entropy')

# dir(lvlo)

# group_tools.conservation_level_score(og.info_dict['orthogroups'][level])

# %%
