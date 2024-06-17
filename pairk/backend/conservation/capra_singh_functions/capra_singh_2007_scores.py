"""
From Capra and Singh 2007, DOI: 10.1093/bioinformatics/btm270
"""

import math
import numpy as np

amino_acids = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
    "-",
]
iupac_alphabet = [
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "U",
    "V",
    "W",
    "Y",
    "Z",
    "X",
    "*",
    "-",
]

# dictionary to map from amino acid to its row/column in a similarity matrix
AA_TO_INDEX = {}
for i, aa in enumerate(amino_acids):
    AA_TO_INDEX[aa] = i

PSEUDOCOUNT = 0.0000001


def read_fasta_alignment(filename):
    """Read in the alignment stored in the FASTA file, filename. Return two
    lists: the identifiers and sequences."""
    f = open(filename)

    names = []
    alignment = []
    cur_seq = ""

    for line in f:
        line = line[:-1]
        if len(line) == 0:
            continue
        if line[0] == ";":
            continue
        if line[0] == ">":
            names.append(line[1:].replace("\r", ""))

            if cur_seq != "":
                cur_seq = cur_seq.upper()
                for i, aa in enumerate(cur_seq):
                    if aa not in iupac_alphabet:
                        cur_seq = cur_seq.replace(aa, "-")
                alignment.append(
                    cur_seq.replace("B", "D").replace("Z", "Q").replace("X", "-")
                )
                cur_seq = ""
        elif line[0] in iupac_alphabet:
            cur_seq += line.replace("\r", "")

    # add the last sequence
    cur_seq = cur_seq.upper()
    for i, aa in enumerate(cur_seq):
        if aa not in iupac_alphabet:
            cur_seq = cur_seq.replace(aa, "-")
    alignment.append(cur_seq.replace("B", "D").replace("Z", "Q").replace("X", "-"))
    f.close()
    return names, alignment


def calculate_sequence_weights(msa):
    """Calculate the sequence weights using the Henikoff '94 method
    for the given msa."""

    # one element for each sequence in the msa (one weight per sequence)
    seq_weights = [0.0] * len(msa)

    # for each column
    for i in range(len(msa[0])):
        freq_counts = [0] * len(amino_acids)
        col = []
        # for each sequence in the alignment
        for j in range(len(msa)):
            if msa[j][i] != "-":  # ignore gaps
                freq_counts[AA_TO_INDEX[msa[j][i]]] += 1

        num_observed_types = 0
        for j in range(len(freq_counts)):
            if freq_counts[j] > 0:
                num_observed_types += 1

        for j in range(len(msa)):
            d = freq_counts[AA_TO_INDEX[msa[j][i]]] * num_observed_types
            if d > 0:
                seq_weights[j] += 1.0 / d

    for w in range(len(seq_weights)):
        seq_weights[w] /= len(msa[0])

    return seq_weights


def weighted_gap_penalty(col, seq_weights=None):
    """Calculate the simple gap penalty multiplier for the column. If the
    sequences are weighted, the gaps, when penalized, are weighted
    accordingly."""
    if seq_weights == None:
        seq_weights = [1.0] * len(col)
    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.0] * len(col)
    gap_sum = 0.0
    for i in range(len(col)):
        if col[i] == "-":
            gap_sum += seq_weights[i]
    return 1 - (gap_sum / sum(seq_weights))


def weighted_freq_count_pseudocount(col, seq_weights=None, pc_amount=0.0000001):
    """Return the weighted frequency count for a column--with pseudocount."""
    if seq_weights == None:
        seq_weights = [1.0] * len(col)
    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.0] * len(col)
    aa_num = 0
    freq_counts = len(amino_acids) * [pc_amount]  # in order defined by amino_acids
    for aa in amino_acids:
        for j in range(len(col)):
            if col[j] == aa:
                freq_counts[aa_num] += 1 * seq_weights[j]
        aa_num += 1
    for j in range(len(freq_counts)):
        freq_counts[j] = freq_counts[j] / (
            sum(seq_weights) + len(amino_acids) * pc_amount
        )
    return freq_counts


def property_entropy(col, seq_weights=None, gap_penalty=1):
    """Calculate the entropy of a column col relative to a partition of the
    amino acids.
    ---- New custom docstring ----
    Williamson '95 partitioning is used. If gap_penalty == 1, then gaps are
    penalized. The entropy will be between inverted so that it spans 0 to 1,
    with 1 being the most conserved.
    partitioning:
    ["V", "L", "I", "M"],
    ["F", "W", "Y"],
    ["S", "T"],
    ["N", "Q"],
    ["H", "K", "R"],
    ["D", "E"],
    ["A", "G"],
    ["P"],
    ["C"],
    ["-"],
    """

    # Mirny and Shakn. '99
    # property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]
    # Williamson '95
    property_partition = [
        ["V", "L", "I", "M"],
        ["F", "W", "Y"],
        ["S", "T"],
        ["N", "Q"],
        ["H", "K", "R"],
        ["D", "E"],
        ["A", "G"],
        ["P"],
        ["C"],
        ["-"],
    ]
    if seq_weights is None:
        seq_weights = [1.0] * len(col)
    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)
    # sum the aa frequencies to get the property frequencies
    prop_fc = [0.0] * len(property_partition)
    for p in range(len(property_partition)):
        for aa in property_partition[p]:
            prop_fc[p] += fc[AA_TO_INDEX[aa]]
    h = 0.0
    for i in range(len(prop_fc)):
        if prop_fc[i] != 0:
            h += prop_fc[i] * math.log(prop_fc[i])
    h /= math.log(min(len(property_partition), len(col)))
    inf_score = 1 - (-1 * h)
    if gap_penalty == 1:
        return inf_score * weighted_gap_penalty(col, seq_weights)
    else:
        return inf_score


def shannon_entropy(col, seq_weights=None, gap_penalty=1):
    """Calculates the Shannon entropy of the column col. sim_matrix  and
    bg_distr are ignored. If gap_penalty == 1, then gaps are penalized. The
    entropy will be between zero and one because of its base. See p.13 of
    Valdar 02 for details. The information score 1 - h is returned for the sake
    of consistency with other scores."""

    if seq_weights is None:
        seq_weights = [1.0] * len(col)
    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

    h = 0.0
    for i in range(len(fc)):
        if fc[i] != 0:
            h += fc[i] * math.log(fc[i])

    #    h /= math.log(len(fc))
    h /= math.log(min(len(fc), len(col)))

    inf_score = 1 - (-1 * h)

    if gap_penalty == 1:
        return inf_score * weighted_gap_penalty(col, seq_weights)
    else:
        return inf_score
