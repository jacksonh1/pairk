import logomaker as lm
import matplotlib.pyplot as plt
import seaborn as sns

amino_acids = sorted(["V","E","K","I","H","L","G","T","M","N","S","P","A","F","W","Y","Q","R","C","D"])

def seqfile2list(seqs_file):
    with open(seqs_file) as handle:
        f = handle.readlines()
        sequence_list = [i.strip() for i in f]
    return sequence_list


def write_seqlist(seqlist, filename):
    with open(filename, 'w') as handle:
        for seq in seqlist:
            handle.write(seq + '\n')


def union_2_lists(l1, l2):
    '''converts each list into a set. finds the union of the 2 sets and returns it as a list'''
    sl1 = set(l1)
    sl2 = set(l2)
    return list(sl1.union(sl2))


def count_sequences(seq_list):
    seq_counts = {i : seq_list.count(i) for i in set(seq_list)}
    return seq_counts


def alignment_2_counts(sequence_list, file_prefix=False, show_plot=False, logosize=[10,2.5], heatmap=False):
    """
    creates a counts matrix PSSM (logomaker format) from a list of sequences
    if an amino acid is absent at a position in the input sequence list, 
     it adds a 0 count for that entry in the matrix

    uses PSSM_plot to generate sequence logo and heatmap plots if specified

    Parameters
    ----------
    sequence_list : list
        list of sequences (strings) from which to generate PSSM
    file_prefix : :obj:`string`, optional
        name of output file prefix to save PSSM csv file
        False to not save
    show_plot : bool, optional
        if True, plots the sequence logo using logomaker
    heatmap : bool, optional
        if True, plots a heatmap of the counts matrix using seaborn

    Returns
    ------
    PSSM counts matrix (logomaker format)
    """
    AAs = [
        "V",
        "E",
        "K",
        "I",
        "H",
        "L",
        "G",
        "T",
        "M",
        "N",
        "S",
        "P",
        "A",
        "F",
        "W",
        "Y",
        "Q",
        "R",
        "C",
        "D",
    ]
    counts = lm.alignment_to_matrix(sequences=sequence_list, to_type="counts")
    for AA in AAs:
        if AA not in counts.columns:
            counts[AA] = 0.0
    counts = counts.reindex(sorted(counts.columns), axis=1)
    if show_plot:
        PSSM_plot(counts, file_prefix=file_prefix, heatmap=heatmap, logosize=logosize)
    if file_prefix:
        counts.to_csv(file_prefix + ".csv")
    return counts


def PSSM_plot(mat, heatmap=False, file_prefix=False, logosize=[10,2.5]):
    """
    Plots sequence logo and heatmap from logomaker style PSSM

    Parameters
    ----------
    mat
        PSSM matrix (logomaker format)
    heatmap : :obj:`bool`, optional
        whether or not to plot a heatmap of the pssm
    file_prefix : :obj:`string`, optional
        name of output file prefix to save images
        False to not save images
    """
    lm.Logo(mat, color_scheme="chemistry", figsize=logosize)
    if file_prefix:
        plt.savefig(file_prefix + "_logo.png", format="png", dpi=200)
    if heatmap:
        fig, ax = plt.subplots(figsize=[9, 8])
        sns.heatmap(
            mat.T, annot=True, fmt=".2g", linewidth=0.5, ax=ax, annot_kws={"size": 14}
        )
        ax.tick_params(axis=u"both", which=u"both", length=0, labelsize=16)
        ax.set_xlabel("position")
        if file_prefix:
            plt.savefig(file_prefix + "_heatmap.png", format="png", dpi=200)


def PSSM_score_sequence(sequence, PSSM):
    """
    PSSM is a PSSM dataframe of sequence position (rows) vs. residue (columns)
    """
    score = 0
    for pos, AA in enumerate(sequence):
        score = score + PSSM.loc[pos, AA]
    return score


