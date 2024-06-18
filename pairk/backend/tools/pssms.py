import logomaker as lm
import matplotlib.pyplot as plt
import seaborn as sns

amino_acids = sorted(
    [
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
)


def seqfile2list(seqs_file):
    with open(seqs_file) as handle:
        f = handle.readlines()
        sequence_list = [i.strip() for i in f]
    return sequence_list


def write_seqlist(seqlist, filename):
    with open(filename, "w") as handle:
        for seq in seqlist:
            handle.write(seq + "\n")


def union_2_lists(l1, l2):
    """converts each list into a set. finds the union of the 2 sets and returns it as a list"""
    sl1 = set(l1)
    sl2 = set(l2)
    return list(sl1.union(sl2))


def count_sequences(seq_list):
    seq_counts = {i: seq_list.count(i) for i in set(seq_list)}
    return seq_counts


def alignment_2_counts(sequence_list):
    """
    creates a counts matrix PSSM (logomaker format) from a list of sequences
    if an amino acid is absent at a position in the input sequence list,
     it adds a 0 count for that entry in the matrix


    Parameters
    ----------
    sequence_list : list
        list of sequences (strings) from which to generate PSSM

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
    return counts


def PSSM_score_sequence(sequence, PSSM):
    """
    PSSM is a PSSM dataframe of sequence position (rows) vs. residue (columns)
    """
    score = 0
    for pos, AA in enumerate(sequence):
        score = score + PSSM.loc[pos, AA]
    return score
