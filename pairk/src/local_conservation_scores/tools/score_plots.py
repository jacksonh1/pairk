import json
import os
import re
import sys
from pathlib import Path

import logomaker as lm
import matplotlib.pyplot as plt
# plt.style.use("custom_standard")
# plt.style.use("custom_small")
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import Align, AlignIO, Seq, SeqIO

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_env_variables.env_variables as env
import local_seqtools.pssms as pssms


def format_bar_plot(ax, xlabel_sequence: str, bar_ylim=[0, 1], labelsize=16):
    """format bar plot"""
    _ = ax.set_xticks(
        list(range(len(xlabel_sequence))),
        labels=list(xlabel_sequence),
    )
    ax.set_xlim(-0.5, len(xlabel_sequence) - 0.5)
    ax.tick_params(axis="x", which="major", labelsize=labelsize)
    if bar_ylim is not None:
        ax.set_ylim(bar_ylim)
    return ax


def get_non_gap_indexes(hit_seq_in_alignment_str):
    """get list of nongap positions in `hit_seq_in_alignment_str`"""
    return [c for c, i in enumerate(hit_seq_in_alignment_str) if i != "-"]


def strip_gaps_from_slice(alignment_slice, query_seq_in_alignment_str):
    """strip gaps from alignment slice"""
    new_sequences = [""] * (len(alignment_slice))
    non_gap_indices = get_non_gap_indexes(query_seq_in_alignment_str)
    new_query_slice = "".join(
        [query_seq_in_alignment_str[i] for i in non_gap_indices]
    )
    for c in non_gap_indices:
        for j in range(len(new_sequences)):
            new_sequences[j] += alignment_slice[j, c]
    return new_sequences, new_query_slice, non_gap_indices


def plot_alignment_slice_conservation2query(
    query_seq_in_alignment_str: str,
    alignment: Align.MultipleSeqAlignment,
    score_list: list|np.ndarray,
    score_mask: list|None = None,
    slice_coordinates: list|None = None,
    strip_gaps: bool = True,
    bar_ylim=[0, 1],
    axes = None,
    figsize=(20, 5),
    **tick_params,
):
    """plot alignment slice conservation to query"""
    assert len(query_seq_in_alignment_str) == len(score_list)
    if axes is None:
        fig, (ax1, ax2) = plt.subplots(figsize=figsize, nrows=2)
    else:
        ax1, ax2 = axes
    if score_mask is not None:
        assert len(score_list) == len(score_mask)
        score_list = np.array(score_list)
        score_list[[not i for i in score_mask]] = 0
        score_list = list(score_list)
    if slice_coordinates is not None:
        alignment = alignment[:, slice_coordinates[0] : slice_coordinates[1]+1] # type: ignore
        query_seq_in_alignment_str = query_seq_in_alignment_str[
            slice_coordinates[0] : slice_coordinates[1]+1
        ]
        score_list = score_list[slice_coordinates[0] : slice_coordinates[1]+1]
        if score_mask is not None:
            score_mask = score_mask[slice_coordinates[0] : slice_coordinates[1]+1]
    if strip_gaps:
        (
            alignment_str_list,
            query_seq_in_alignment_str,
            non_gap_indices,
        ) = strip_gaps_from_slice(
            alignment, query_seq_in_alignment_str
        )
        score_list = [score_list[i] for i in non_gap_indices]
    else:
        alignment_str_list = [str(i.seq) for i in alignment]
    if bar_ylim == "auto":
        limit = max(abs(min(score_list)), abs(max(score_list)))
        bar_ylim = [-limit, limit]

    ax1.bar(
        list(range(len(query_seq_in_alignment_str))),
        score_list,
    )
    ax1 = format_bar_plot(
        ax1, query_seq_in_alignment_str, bar_ylim=bar_ylim, **tick_params
    )

    counts = pssms.alignment_2_counts(
        alignment_str_list, show_plot=False, heatmap=False
    )
    lm.Logo(counts, color_scheme="chemistry", ax=ax2)
    ax2.set_ylim(0, len(alignment_str_list))
    return ax1, ax2, counts


def bg_score_distro_plot(score_list: list, ax):
    sns.histplot(score_list, bins=20, ax=ax)
    ax.set_xlim([0, 1])
    # add a vertical line at the mean
    ax.axvline(np.mean(score_list), color="k", linewidth=1)
    # add a vertical line at one standard deviation above and below the mean
    ax.axvline(
        np.mean(score_list) + np.std(score_list), # type: ignore
        color="k",
        linestyle="dashed",
        linewidth=1,
    )
    ax.axvline(
        np.mean(score_list) - np.std(score_list), # type: ignore
        color="k",
        linestyle="dashed",
        linewidth=1,
    )

# ==============================================================================
# // z-score plots
# ==============================================================================
def build_og_level_screen_mosaic_z_score(levels):
    mos_vector = []
    for level in levels:
        mos_vector.append([f"bg_dist-{level}"] + [f"scores-{level}"] * 3)
        mos_vector.append([f"bg_dist-{level}"] + [f"logo-{level}"] * 3)
    fig, axd = plt.subplot_mosaic(mos_vector, figsize=(10, 10), layout="constrained")
    plt.tight_layout()
    return fig, axd


def add_z_score_plots_to_mosaic(
    axd, level, lvl_o: group_tools.LevelAlnScore, bar_ylim=[-2.5, 2.5], highlight_positions=None, **plot_kwargs
):
    ax1, ax2, counts = plot_alignment_slice_conservation2query(
        query_seq_in_alignment_str=lvl_o.query_aln_sequence,
        alignment=lvl_o.aln,
        score_list=lvl_o.z_scores, # type: ignore
        score_mask=lvl_o.score_mask,
        slice_coordinates=[
            lvl_o.hit_aln_start,
            lvl_o.hit_aln_end,
        ],
        bar_ylim=bar_ylim,
        axes=[axd[f"scores-{level}"], axd[f"logo-{level}"]],
        **plot_kwargs,
    )
    if highlight_positions is not None:
        for position in highlight_positions:
            for ax in [ax1, ax2]:
                ax.axvspan(
                    position+0.5,
                    position-0.5,
                    color="red",
                    alpha=0.3,
                )


def big_z_score_plot(og: group_tools.ConserGene, **kwargs):
    levels = list(og.info_dict["orthogroups"].keys())
    if all([x in env.PHYLOGENY_LVL_ORDERING for x in levels]):
        levels = sorted(
            levels,
            key=lambda x: env.PHYLOGENY_LVL_ORDERING.index(x),
        )
    fig, axd = build_og_level_screen_mosaic_z_score(levels)
    for level in levels:
        if level not in og.levels_passing_filters:
            message = f"{og.query_gene_id}: not enough sequences"
            axd[f"scores-{level}"].text(
                0.5,
                0.5,
                message,
                horizontalalignment="center",
                verticalalignment="center",
                transform=axd[f"scores-{level}"].transAxes,
                fontsize=11,
            )
            # axd[f'scores-{level}'].set_title(, fontsize=11)
            continue
        if og.aln_score_objects[level].z_score_failure is not None:
            message = f"{og.query_gene_id} - {og.aln_score_objects[level].z_score_failure}: not enough background scores for z-score"
            axd[f"scores-{level}"].text(
                0.5,
                0.5,
                message,
                horizontalalignment="center",
                verticalalignment="center",
                transform=axd[f"scores-{level}"].transAxes,
                fontsize=11,
            )
            continue            
        lvl_o = og.aln_score_objects[level]
        add_z_score_plots_to_mosaic(axd, level, lvl_o, **kwargs)
        bg_score_distro_plot(
            lvl_o.bg_scores, axd[f"bg_dist-{level}"] # type: ignore
        )
        plot_title = f"{og.query_gene_id} - {level} - {len(lvl_o.aln)} sequences - z-score (IDR bg using {len(lvl_o.bg_scores)} residues)" # type: ignore
        axd[f"scores-{level}"].set_title(plot_title, fontsize=11)
        axd[f"bg_dist-{level}"].set_title(f"{level} - {og.query_gene_id}")
    return fig, axd


# ==============================================================================
# // just score plots (not z-scores)
# ==============================================================================
def build_og_level_screen_mosaic(levels):
    mos_vector = []
    for level in levels:
        mos_vector.append([f"scores-{level}"] * 3)
        mos_vector.append([f"logo-{level}"] * 3)
    fig, axd = plt.subplot_mosaic(mos_vector, figsize=(10, 10), layout="constrained")
    plt.tight_layout()
    return fig, axd


def add_score_plots_to_mosaic(
    axd, level, lvl_o, bar_ylim=[0, 1], highlight_positions=None, **plot_kwargs
):
    ax1, ax2, counts = plot_alignment_slice_conservation2query(
        query_seq_in_alignment_str=lvl_o.query_aln_sequence,
        alignment=lvl_o.aln,
        score_list=lvl_o.scores,
        score_mask=lvl_o.score_mask,
        slice_coordinates=[
            lvl_o.hit_aln_start,
            lvl_o.hit_aln_end,
        ],
        bar_ylim=bar_ylim,
        axes=[axd[f"scores-{level}"], axd[f"logo-{level}"]],
        **plot_kwargs,
    )
    if highlight_positions is not None:
        for position in highlight_positions:
            for ax in [ax1, ax2]:
                ax.axvspan(
                    position+0.5,
                    position-0.5,
                    color="red",
                    alpha=0.3,
                )


def big_score_plot(og: group_tools.ConserGene, **kwargs):
    levels = list(og.info_dict["orthogroups"].keys())
    if all([x in env.PHYLOGENY_LVL_ORDERING for x in levels]):
        levels = sorted(
            levels,
            key=lambda x: env.PHYLOGENY_LVL_ORDERING.index(x),
        )
    fig, axd = build_og_level_screen_mosaic(levels)
    for level in levels:
        if level not in og.levels_passing_filters:
            message = f"{og.query_gene_id}: not enough sequences"
            axd[f"scores-{level}"].text(
                0.5,
                0.5,
                message,
                horizontalalignment="center",
                verticalalignment="center",
                transform=axd[f"scores-{level}"].transAxes,
                fontsize=11,
            )
            # axd[f'scores-{level}'].set_title(, fontsize=11)
            continue
        lvl_o = og.aln_score_objects[level]
        add_score_plots_to_mosaic(axd, level, lvl_o, **kwargs)
        plot_title = f"{og.query_gene_id} - {level} - {len(lvl_o.aln)} sequences - NOT Z-SCORES"
        axd[f"scores-{level}"].set_title(plot_title, fontsize=11)
    return fig, axd


# %%
# ==============================================================================
# // big plot driver
# ==============================================================================
def ogscreen_plot_from_og_folder_class(
    og: group_tools.ConserGene, big_plot_output_folder, save_big_plot=True, score_name=None, score_type='zscore', **kwargs
):
    if score_type == 'zscore':
        fig, axd=big_z_score_plot(og, **kwargs)
    else:
        fig, axd=big_score_plot(og, **kwargs)
    big_plot_filename = (
        f"{og.reference_index}-{og.query_gene_id}_og_level_score_screen.png"
    )
    if save_big_plot:
        big_plot_output_folder = Path(big_plot_output_folder)
        big_plot_output_folder.mkdir(parents=True, exist_ok=True)
        if score_name is None:
            big_plot_filename = (
                big_plot_output_folder
                / f"{og.reference_index}-{og.query_gene_id.replace(':','')}_og_level_score_screen.png"
            )
        else:
            big_plot_filename = (
                big_plot_output_folder
                / f"{og.reference_index}-{og.query_gene_id.replace(':','')}-{score_name}_og_level_score_screen.png"
            )
        fig.savefig(big_plot_filename, bbox_inches="tight", dpi=300)
    return big_plot_filename
    # return fig, axd, big_plot_filename


# %%
# ==============================================================================
# // TITLE
# ==============================================================================
def build_mosaic_z_score_plot():
    mos_vector = []
    mos_vector.append(["bg_dist"] + ["scores"]*2)
    mos_vector.append(["bg_dist"] + ["logo"]*2)
    fig, axd = plt.subplot_mosaic(mos_vector, figsize=(15, 5), layout="constrained")
    plt.tight_layout()
    return fig, axd

def plot_score_bar_plot(ax, score_list, query_seq, mask=None):
    if mask is not None:
        score_list = np.array(score_list)
        score_list[[not i for i in mask]] = 0
        score_list = list(score_list)
    ax.bar(
        list(range(len(score_list))),
        score_list,
    )
    ax = _format_bar_plot(ax, query_seq)
    return ax

def _format_bar_plot(ax, xlabel_sequence: str, labelsize=16):
    """format bar plot"""
    _ = ax.set_xticks(
        list(range(len(xlabel_sequence))),
        labels=list(xlabel_sequence),
    )
    ax.set_xlim(-0.5, len(xlabel_sequence) - 0.5)
    ax.tick_params(axis="x", which="major", labelsize=labelsize)
    return ax

def plot_logo(ax, str_list, tick_label_str, labelsize=16):
    counts = pssms.alignment_2_counts(
        str_list, show_plot=False, heatmap=False
    )
    lm.Logo(counts, color_scheme="chemistry", ax=ax)
    ax.set_ylim(0, len(str_list))
    _=ax.set_xticks(
        list(range(len(str_list[0]))),
        labels=list(tick_label_str),
    )
    ax.tick_params(axis="x", which="major", labelsize=labelsize)
    return ax

def format_logo_xticks_with_str(ax, tick_label_str):
    _=ax.set_xticks(
        list(range(len(tick_label_str))),
        labels=list(tick_label_str),
    )
    return ax
