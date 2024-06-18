import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import matplotlib.axes


def make_empty_kmer_ortho_df(positions, ortholog_ids: list[str]):
    cols = ["query_kmer"] + ortholog_ids
    return pd.DataFrame(
        index=positions,
        columns=cols,
    )


def matrix_dict_2_df(matrix_dict, mat_key):
    return pd.DataFrame(
        matrix_dict[mat_key]["data"],
        columns=matrix_dict[mat_key]["columns"],
        index=matrix_dict[mat_key]["index"],
    )


def import_pairwise_matrices(filepath, matrix_keys: list[str] | None = None):
    if matrix_keys is None:
        matrix_keys = [
            "score_matrix",
            "orthokmer_matrix",
            "position_matrix",
        ]
    with open(filepath, "r") as json_file:
        data = json.load(json_file)
    matrices = {}
    for k in matrix_keys:
        if k in data:
            matrices[k] = matrix_dict_2_df(data, k)
    if len(matrices) == 0:
        raise ValueError("No matrices found in the json file")
    return matrices


class PairkAln:
    """A class to store the results of the pairwise alignment.

    The primary data is stored in pandas dataframes. All dataframes have the
    same structure. One column is the query k-mer sequence
    ('query_kmer'). The other columns are named as the ortholog sequence
    ids. The dataframe indexes are the query k-mer start position in the
    query sequence.

    Attributes
    ----------
    orthokmer_matrix : pd.DataFrame
        the best scoring k-mer from each ortholog for each query k-mer.
    position_matrix : pd.DataFrame
        the start position of the best scoring k-mer from each ortholog for
        each query k-mer.
    score_matrix : pd.DataFrame | None
        the alignment scores for each k-mer in the query sequence against
        the corresponding best matching ortholog k-mer.
    query_kmers : list[str]
        the list of query k-mers that were aligned.
    query_sequence : str
        the full query sequence that was originally split into k-mers and aligned.
    k : int
        the k-mer size used for the alignment.
    """

    def __init__(
        self,
        orthokmer_df: pd.DataFrame,
        pos_df: pd.DataFrame,
        score_df: pd.DataFrame | None = None,
    ):
        self.orthokmer_matrix = orthokmer_df
        self.position_matrix = pos_df
        self.score_matrix = score_df
        # a little trick to get the query sequence from the orthokmer_matrix
        query_kmers = self.orthokmer_matrix["query_kmer"].to_list()
        self.query_kmers = query_kmers
        self.query_sequence = "".join([i[0] for i in query_kmers]) + query_kmers[-1][1:]
        self.k = len(query_kmers[0])

    def __str__(self):
        return (
            f"PairkAln object for {len(self.query_kmers)} query k-mers\n"
            f"query sequence: {self.query_sequence}\n"
            f"k-mer length: {self.k}\n"
        )

    @classmethod
    def from_file(cls, filepath: str | Path):
        """import the pairwise alignment matrices from a json file.

        Parameters
        ----------
        filepath : str|Path
            the path to the json file containing the pairwise alignment matrices.

        Returns
        -------
        Pairk.PairkAln
            PairkAln object containing the pairwise alignment matrices.
        """
        matrices = import_pairwise_matrices(filepath)
        return cls(
            orthokmer_df=matrices["orthokmer_matrix"],
            pos_df=matrices["position_matrix"],
            score_df=matrices.get("score_matrix"),
        )

    def get_pseudo_alignment(self, position: int):
        """get a list of the best scoring k-mers from each ortholog at a given query position.

        Parameters
        ----------
        position : int
            the position of the query k-mer in the query sequence (0-based index).

        Returns
        -------
        list[str]
            list of the best scoring k-mers from each ortholog for the query k-mer.
        """
        return self.orthokmer_matrix.loc[position, :].to_list()  # type: ignore

    def find_query_kmer_positions(self, kmer: str):
        """convenience function to search for the positions of a k-mer string.

        Parameters
        ----------
        kmer : str
            the k-mer string to search for.

        Returns
        -------
        list[int]
            the positions in the query sequence that match the input kmer.
        """
        return [i for i, x in enumerate(self.query_kmers) if x == kmer]

    def write_to_file(self, filepath: str | Path) -> str | Path:
        """save the PairkAln object matrices to a json file.

        Parameters
        ----------
        filepath : str | Path
            file path to save the json file.

        Returns
        -------
        str | Path
            file path of the saved json file.
        """
        output_dict = {
            "orthokmer_matrix": self.orthokmer_matrix.to_dict(orient="split"),
            "position_matrix": self.position_matrix.to_dict(orient="split"),
        }
        if self.score_matrix is not None:
            output_dict["score_matrix"] = self.score_matrix.to_dict(orient="split")
        with open(filepath, "w") as json_file:
            json.dump(
                output_dict,
                json_file,
            )
        return filepath

    def _create_axes_if_none(
        self, ax: matplotlib.axes.Axes | None = None
    ) -> matplotlib.axes.Axes:
        if ax is None:
            fig, ax = plt.subplots()
        return ax  # type: ignore

    def plot_score_heatmap(
        self, ax: matplotlib.axes.Axes | None = None, **kwargs
    ) -> matplotlib.axes.Axes:
        """plots a heatmap of the alignment scores.

        Parameters
        ----------
        ax : matplotlib.axes.Axes | None, optional
            The axis on which to plot the heatmap, by default None. if None, a new figure is created.
        kwargs : optional
            additional keyword arguments to pass to the seaborn heatmap function.

        Returns
        -------
        matplotlib.axes.Axes
            The axis on which the heatmap was plotted.

        Raises
        ------
        ValueError
            raised if no score matrix is found.
        """
        ax = self._create_axes_if_none(ax)
        if self.score_matrix is None:
            raise ValueError("No score matrix found.")
        df = self.score_matrix.copy()
        df = df.drop(columns="query_kmer")
        df = df.astype(float)
        if "cmap" not in kwargs:
            kwargs["cmap"] = "viridis"
        sns.heatmap(df, ax=ax, **kwargs)
        return ax  # type: ignore

    def plot_position_heatmap(
        self, ax: matplotlib.axes.Axes | None = None, **kwargs
    ) -> matplotlib.axes.Axes:
        """plot a heatmap of the start positions of the best scoring k-mers in each ortholog.

        Parameters
        ----------
        ax : matplotlib.axes.Axes | None, optional
            The axis to plot the heatmap on, by default None. If None, a new figure is created.
        kwargs : optional
            additional keyword arguments to pass to the seaborn heatmap function.

        Returns
        -------
        matplotlib.axes.Axes
            The axis on which the heatmap was plotted.
        """
        ax = self._create_axes_if_none(ax)
        df = self.position_matrix.copy()
        df = df.drop(columns="query_kmer")
        df = df.astype(float)
        if "cmap" not in kwargs:
            kwargs["cmap"] = "viridis"
        sns.heatmap(df, ax=ax, **kwargs)
        return ax  # type: ignore
