import pairk.backend.tools.pssms as pssms
import matplotlib.pyplot as plt

plt.style.use("pairk.data.pairk_plotstyle")
import matplotlib.pyplot as plt
import matplotlib.axes
import seaborn as sns
import logomaker as lm


def build_mosaic_z_score_plot(figsize=(15, 5)):
    mos_vector = []
    mos_vector.append(["background"] + ["scores"] * 2)
    mos_vector.append(["background"] + ["logo"] * 2)
    fig, axd = plt.subplot_mosaic(mos_vector, figsize=figsize, layout="constrained")
    # plt.tight_layout()
    return fig, axd


def plot_score_bar_plot(ax: matplotlib.axes.Axes, score_list: list, query_seq: str):
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
    # ax.tick_params(axis="x", which="major", labelsize=labelsize)
    return ax


def plot_logo(ax, str_list, tick_label_str, labelsize=16):
    counts = pssms.alignment_2_counts(str_list)
    lm.Logo(counts, color_scheme="chemistry", ax=ax)
    ax.set_ylim(0, len(str_list))
    _ = ax.set_xticks(
        list(range(len(str_list[0]))),
        labels=list(tick_label_str),
    )
    # ax.tick_params(axis="x", which="major", labelsize=labelsize)
    return ax


def format_logo_xticks_with_str(ax, tick_label_str):
    _ = ax.set_xticks(
        list(range(len(tick_label_str))),
        labels=list(tick_label_str),
    )
    return ax
