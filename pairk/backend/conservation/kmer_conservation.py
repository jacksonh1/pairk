from pairk.backend.conservation.capra_singh_functions.capra_singh_2007_scores import (
    property_entropy,
)
import pandas as pd
import numpy as np
from typing import Callable
from pairk.backend.tools.pairwise_tools import PairkAln

# import pairk.backend.tools.pssms as pssms
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.axes
import matplotlib.figure
import pairk.backend.tools.plotting_tools as plotting_tools

plt.style.use("pairk.backend.pairk_plotstyle")


# %%
def kmerlist_to_columns(seqlist: np.ndarray):
    col_strings = []
    for c in range(len(seqlist[0])):
        col = "".join([seqlist[i][c] for i in range(len(seqlist))])
        col_strings.append(col)
    return col_strings


def score_pseudo_aln(
    seqlist: np.ndarray,
    score_func: Callable = property_entropy,
):
    col_strings = kmerlist_to_columns(seqlist)
    return [score_func(c) for c in col_strings]


def orthokmer_arr_to_score_arr(
    orthokmer_arr: np.ndarray, score_func: Callable = property_entropy
) -> np.ndarray:
    k = len(orthokmer_arr[0, 0])
    score_arr = np.zeros((orthokmer_arr.shape[0], k))
    for i in range(orthokmer_arr.shape[0]):
        score_arr[i, :] = score_pseudo_aln(orthokmer_arr[i, :], score_func)
    return score_arr


def calculate_z_scores(score_arr: np.ndarray):
    bg_scores = score_arr.flatten()
    bg_mean = bg_scores.mean()
    bg_std = bg_scores.std()
    z_score_arr = (score_arr - bg_mean) / bg_std
    return z_score_arr


def calculate_conservation_arrays(
    orthokmer_df: pd.DataFrame,
    score_func: Callable = property_entropy,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not all([type(i) == str for i in orthokmer_df.values.flatten()]):
        raise ValueError("Orthokmer matrix contains non-string elements")
    # sort orthokmer_df by the index
    orthokmer_df = orthokmer_df.copy().sort_index()
    orthokmer_arr = orthokmer_df.copy().values
    score_arr = orthokmer_arr_to_score_arr(orthokmer_arr, score_func)
    z_score_arr = calculate_z_scores(score_arr)
    return orthokmer_arr, score_arr, z_score_arr


class PairkConservation:
    """a class to store the results of the conservation scoring

    The methods can be used to create plots of the conservation scores and
    sequence logos.

    Attributes
    ----------
    orthokmer_arr : np.ndarray
        the best scoring k-mer from each ortholog for each query k-mer.
    score_arr : np.ndarray
        the conservation scores for each k-mer position.
    z_score_arr : np.ndarray
        the z-scores for each k-mer position.
    query_kmers : list[str]
        the query k-mers.
    query_sequence : str
        the query sequence.
    k : int
        the length of the query k-mers.
    bg_scores : np.ndarray
        the background conservation scores used to calculate the z-scores.
    n_bg_scores : int
        the number of background scores used to calculate the z-scores.
    n_bg_kmers : int
        the number of k-mers used to calculate the z-scores.
    bg_mean : float
        the mean of the background scores.
    bg_std : float
        the standard deviation of the background scores.
    """

    def __init__(
        self,
        orthokmer_arr: np.ndarray,
        score_arr: np.ndarray,
        z_score_arr: np.ndarray,
    ):
        self.orthokmer_arr = orthokmer_arr
        self.score_arr = score_arr
        self.z_score_arr = z_score_arr
        self.query_kmers = list(self.orthokmer_arr[:, 0])
        self.query_sequence = (
            "".join([i[0] for i in self.query_kmers]) + self.query_kmers[-1][1:]
        )
        self.k = len(self.query_kmers[0])
        self.bg_scores = self.score_arr.flatten()
        self.n_bg_scores = len(self.bg_scores)
        self.n_bg_kmers = self.orthokmer_arr.shape[0]
        self.bg_mean = self.bg_scores.mean()
        self.bg_std = self.bg_scores.std()

    def __repr__(self):
        return (
            f"PairkConservation object\n"
            f"{self.orthokmer_arr.shape[0]} query kmers x {self.orthokmer_arr.shape[1]} orthologs\n"
            f"kmer length - {self.k}\n"
            f"{self.n_bg_scores} background scores for z-score calculation\n"
            f"{self.n_bg_kmers} kmers used for z-score calculation"
        )

    def get_average_score(
        self,
        position: int,
        score_type: str = "z_score",
        position_mask: np.ndarray | None = None,
    ):
        if score_type == "score":
            scores = np.copy(self.score_arr[position, :])
        elif score_type == "z_score":
            scores = np.copy(self.z_score_arr[position, :])
        else:
            raise ValueError("score_type must be 'score' or 'z_score'")
        if position_mask is None:
            position_mask = np.ones(scores.shape)
        maskarr = np.array(position_mask)
        mask = maskarr.nonzero()
        masked_scores = scores[mask]
        return np.mean(masked_scores)

    @staticmethod
    def print_array_hist(a: np.ndarray, bins: np.ndarray | int = 30):
        hist, bin_edges = np.histogram(a.flatten(), bins=bins)
        for i, h in enumerate(hist):
            print(f"{bin_edges[i]:.2f} - {bin_edges[i+1]:.2f}: {'=' * h}")

    def _create_axes_if_none(
        self, ax: matplotlib.axes.Axes | None = None
    ) -> matplotlib.axes.Axes:
        if ax is None:
            fig, ax = plt.subplots(figsize=(4, 4))
        return ax  # type: ignore

    def plot_background_distribution(
        self,
        ax: matplotlib.axes.Axes | None = None,
        bins: int = 20,
    ):
        ax = self._create_axes_if_none(ax)
        ax.hist(self.bg_scores, bins=bins)
        # ax.set_title("Background")
        ax.set_xlabel("Conservation score")
        ax.set_ylabel("Count")
        ax.axvline(self.bg_mean, color="red", linestyle="--", label="Mean", linewidth=2)
        ax.axvline(
            self.bg_mean + self.bg_std, color="black", label="Mean + 1 std", linewidth=2
        )
        ax.axvline(
            self.bg_mean - self.bg_std, color="black", label="Mean - 1 std", linewidth=2
        )
        # ax.legend()
        ax.set_xlim(0, 1)
        return ax

    def plot_score_barplot(
        self,
        position: int,
        score_type: str = "score",
        ax: matplotlib.axes.Axes | None = None,
    ):
        ax = self._create_axes_if_none(ax)
        if score_type == "score":
            scores = np.copy(self.score_arr[position, :])
        elif score_type == "z_score":
            scores = np.copy(self.z_score_arr[position, :])
        else:
            raise ValueError("score_type must be 'score' or 'z_score'")
        query_kmer = self.query_kmers[position]
        ax = plotting_tools.plot_score_bar_plot(ax, list(scores), query_kmer)
        return ax

    def plot_sequence_logo(
        self,
        position: int,
        ax: matplotlib.axes.Axes | None = None,
    ):
        ax = self._create_axes_if_none(ax)
        pseudo_aln = list(self.orthokmer_arr[position, :])
        query_kmer = self.query_kmers[position]
        ax = plotting_tools.plot_logo(ax, pseudo_aln, query_kmer)
        return ax

    def plot_conservation_mosaic(
        self,
        position: int,
        score_type: str = "z_score",
        figsize: tuple[int, int] = (15, 5),
    ) -> tuple[matplotlib.figure.Figure, dict[str, matplotlib.axes.Axes]]:
        """makes a mosaic plot of the conservation scores, sequence logos, and background scores for the pairk conservation results.

        Parameters
        ----------
        position : int
            starting position of the k-mer in the query sequence.
        score_type : str, optional
            either 'score' or 'z_score'. The type of score to plot on the bar
            plot, by default "score"

        Returns
        -------
        tuple[matplotlib.figure.Figure, dict[str, matplotlib.axes.Axes]]
            the figure and axes dictionary for the mosaic plot.
        """
        fig, axd = plotting_tools.build_mosaic_z_score_plot(figsize=figsize)
        axd["background"] = self.plot_background_distribution(axd["background"])
        axd["scores"] = self.plot_score_barplot(position, score_type, axd["scores"])
        axd["logo"] = self.plot_sequence_logo(position, axd["logo"])
        return fig, axd

    def write_results_to_file(self, filename: str | Path):
        # save np arrays to a file
        # to avoid having to pickle the np arrays, we convert the orthokmer_arr
        # to numpy strings (np.str_) before saving
        np.savez(
            filename,
            orthokmer_arr=self.orthokmer_arr.astype(np.str_),
            score_arr=self.score_arr,
            z_score_arr=self.z_score_arr,
        )
        return filename

    @classmethod
    def read_results_from_file(cls, filename: str | Path):
        # load np arrays from a file
        npzfile = np.load(filename)
        orthokmer_arr = npzfile["orthokmer_arr"]
        # convert the orthokmer_arr back to python strings.
        # it probably doesn't matter, but I want to keep everything consistent
        orthokmer_arr = orthokmer_arr.astype(str)
        score_arr = npzfile["score_arr"]
        z_score_arr = npzfile["z_score_arr"]
        return cls(orthokmer_arr, score_arr, z_score_arr)


# %%


# write a function that does the same thing as the above code
def calculate_conservation(
    pairk_aln_results: PairkAln,
    score_func: Callable = property_entropy,
) -> PairkConservation:
    """calculate the conservation scores for the k-mers in the PairkAln object. calculates the conservation scores and z-scores for each k-mer position.

    Parameters
    ----------
    pairk_aln_results : PairkAln
        the results of the pairk alignment step as a pairk.PairkAln object.
    score_func : Callable, optional
        A function to calculate conservation scores in a columnwise manner, by
        default it is the property_entropy function from Capra and Singh 2007,
        DOI: 10.1093/bioinformatics/btm270 located in the
        `pairk.pairk_conservation.capra_singh_functions` module.

    Returns
    -------
    PairkConservation
        PairkConservation object containing the conservation scores and z-scores for each k-mer position.
    """
    ok_arr, score_arr, z_score_arr = calculate_conservation_arrays(
        pairk_aln_results.orthokmer_matrix, score_func
    )
    return PairkConservation(ok_arr, score_arr, z_score_arr)


# %%
