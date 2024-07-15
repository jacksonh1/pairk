import copy
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

# from Bio.Align import MultipleSeqAlignment


def percent_identity(seq1: SeqRecord, seq2: SeqRecord) -> float:
    """
    Returns a percent identity between 0 (no identical residues) and 1 (all residues are identical)
    - The sequences must be pre-aligned (i.e. they have the same length)
    - The returned percent identity is computed as the number of identical residues divided by the
    length of compared alignment sites. Alignment sites where both sequences are gaps ('-' characters)
    are not compared.

    Parameters
    ----------
    seq1 : str or seqrecord
        first sequence
    seq2 : str or seqrecord
        second sequence
    """
    assert len(seq1) == len(
        seq2
    ), "sequences are not the same length. Are they aligned?"
    num_same = 0
    length = len(seq1)
    for i in range(len(seq1)):
        if seq1[i] == "-" and seq2[i] == "-":
            length -= 1
            continue
        if seq1[i] == seq2[i]:
            num_same += 1
    pid = num_same / length
    return pid


def get_first_non_gap_index(s: str) -> int:
    """get the index of the first non-gap character (first character that is
    not `-`) in a string

    Parameters
    ----------
    s : str
        input string

    Returns
    -------
    int
        the index of the first non-gap character in the string (0-indexed)
    """
    index = 0
    while index < len(s) and s[index] == "-":
        index += 1
    return index


def sort_aln_by_pid2ref(aln: Align.MultipleSeqAlignment, refseq: SeqRecord):
    aln.sort(
        key=lambda record: percent_identity(record, refseq),
        reverse=True,
    )
    return aln


def subsample_seqrec_list(seqrec_list, fraction=0.1):
    """return a subsampled list of seqrecords"""
    n = int(len(seqrec_list) * fraction)
    interval = int(len(seqrec_list) / n)
    new_list = []
    count = 0
    for seq in seqrec_list:
        if count % interval == 0:
            new_list.append(seq)
        count += 1
    return new_list


def subsample_seqrec_list_target_num(
    seqrec_list, max_len=9999999999, target_num=10, skip_list=None
):
    """return a subsampled list of seqrecords"""
    interval = int(len(seqrec_list) / target_num)
    if skip_list is None:
        skip_list = []
    new_list = []
    count = 0
    for seq in seqrec_list:
        if count % interval == 0:
            if seq.id in skip_list:
                continue
            if len(seq.seq) > max_len:
                continue
            new_list.append(seq)
        count += 1
    return new_list


def pad_hit(seq: str, st_pos: int, end_pos: int, l_flank: int = 0, r_flank: int = 0):
    st = max(0, st_pos - l_flank)
    end = min(len(seq) - 1, end_pos + r_flank)
    return st, end, seq[st : end + 1]


def z_score(scores):
    scores = np.array(scores)
    return (scores - np.mean(scores)) / np.std(scores)


def z_score_comparison(scores, scores_ref):
    u = np.mean(scores_ref)
    s = np.std(scores_ref)
    if s == 0:
        print("standard deviation is 0")
    return [(i - u) / s for i in scores]


def calculate_z_score_bg_region(
    score_list,
    score_list_mask,
    bg_slicing_region=None,
    num_bg_scores_cutoff=20,
):
    """calculates the z-score of a score list compared to the background across a region of interest

    Uses the score_list_mask choose which positions to use for the background scores. I do not want to include extremely gappy regions in the background calculation, so I use the score_list_mask to choose which positions to use for the background calculation. The background is calculated as the mean and standard deviation of the background scores. The z-score is calculated as (score - mean) / std

    Parameters
    ----------
    score_list : list
        list of scores. Usually for the entire sequence
    score_list_mask : list
        boolean mask. Must be the same length as the `score_list`. This is primarally used to mask the positions used to calculate the background of the z-score
    bg_slicing_region : list
        region of interest to use as the background ([start_position, end_position]). Numbering must coorespond to the `score_list` and `score_list_mask`
    num_bg_scores_cutoff : int
        minimum number of background scores to calculate the z-score. If there are not enough background scores, then the z-score is not calculated

    Returns
    -------
    dictionary
        dictionary of the z-score, and the background scores (keys: "z_scores", "bg_scores")
    """
    score_list_mask = np.array(score_list_mask)
    if bg_slicing_region is not None:
        score_list_mask = score_list_mask.copy()
        # set mask to False for positions outside of the slicing region
        score_list_mask[: bg_slicing_region[0]] = False
        score_list_mask[bg_slicing_region[1] + 1 :] = False
    bg_scores = np.array(score_list)[score_list_mask]
    if len(bg_scores) < num_bg_scores_cutoff:
        raise ValueError(
            f"not enough background scores to calculate z-score. Require at least {num_bg_scores_cutoff} background scores. Only have {len(bg_scores)} background scores"
        )
    z_scores = z_score_comparison(score_list, bg_scores)
    return {"z_scores": list(z_scores), "bg_scores": list(bg_scores)}


def get_non_gap_indexes(aln_seq: str) -> list[int]:
    """get list of nongap positions in `aln_seq`"""
    return [c for c, i in enumerate(aln_seq) if i != "-"]


# ==============================================================================
# // sequence tools
# ==============================================================================
def str2seqrecord(seq_str, id=None, description=None):
    """convert a string to a SeqRecord object

    Parameters
    ----------
    seq_str : str
        sequence string
    id : str, optional
        id of the SeqRecord object, by default None
    description : str, optional
        description of the SeqRecord object, by default None

    Returns
    -------
    SeqRecord
        SeqRecord object with the sequence string
    """
    if description is None:
        description = ""
    return SeqIO.SeqRecord(Seq.Seq(seq_str), id=id, description=description)  # type: ignore


def pad_with_aas_or_gaps(seq: str, st: int, end: int, flank: int = 0) -> str:
    """
    slice a sequence at the specified positions
    flank the slice with the specified number of residues to the left and right
    if the flank would go out of bounds, slice to the end/beginning of the sequence and pad with dashes

    Parameters
    ----------
    seq : str
        sequence to be sliced
    st : int
        start position of the slice
    end : int
        end position of the slice
    flank : int, optional
        number of residues to flank the slice with, by default 0

    Returns
    -------
    str
        sliced sequence

    Examples
    --------
    >>> slice_seq('ABCDEF', 2, 4, flank=1)
    >>> s = '012FPpPp890'
    >>> st=10
    >>> end=7
    >>> flank=0
    >>> print(s[st:end+1])
    >>> for c, i in enumerate(slice_seq(s, st, end+1, flank)):
            print(c, i)
    >>> print(slice_seq(s, st, end+1, flank))

    """
    assert st <= end, "start must be less than or equal to end"
    assert flank >= 0, "flank must be positive"

    if st - flank < 0:
        left = "-" * abs(flank - st) + seq[: max(0, st)]
    else:
        left = seq[st - flank : st]
    if end + flank > len(seq):
        right = seq[max(0, st) : len(seq)] + "-" * (end + flank - len(seq))
    else:
        right = seq[max(0, st) : end + flank]
    # print(left)
    # print(right)
    return left + right


def gen_kmers(seq: str, k: int) -> list[str]:
    """
    generates list of length k "kmers" comprising `seq`

    Parameters
    ----------
    seq : str
        input sequence to split into k-mers. Should be entry in fasta file
    k : int
        length of fragments (k-mers)

    Returns
    -------
    kmers : list
        list of strings - k-mers generated from seq
    """
    k2 = k - 1
    kmers = []
    for i in range(len(seq) - k2):
        kmers.append(seq[i : i + k])
    return kmers


class FastaImporter:
    """import fasta file and return seqrecord objects in various formats

    Parameters
    ----------
    fasta_path : str
        file path to fasta file
    """

    def __init__(self, fasta_path: str | Path):
        self.fasta_path = fasta_path

    def import_as_list(self) -> list[SeqRecord]:
        """return list of SeqRecord objects for each sequence in the fasta file

        Returns
        -------
        List[SeqRecord]
            list of SeqRecord objects
        """
        with open(self.fasta_path) as handle:
            return list(SeqIO.parse(handle, "fasta"))

    def import_as_dict(self) -> dict[str, SeqRecord]:
        """return dictionary of SeqRecord objects for each sequence in the fasta file

        Returns
        -------
        dict[str, SeqRecord]
            dictionary of SeqRecord objects, keys are the sequence ids and values are the SeqRecord objects
        """
        with open(self.fasta_path) as handle:
            return SeqIO.to_dict(SeqIO.parse(handle, "fasta"))

    def import_as_str_dict(self) -> dict[str, str]:
        """return dictionary of strings for each sequence in the fasta file

        Returns
        -------
        dict[str, str]
            dictionary of sequence strings, keys are the sequence ids and values are the sequences as strings
        """
        with open(self.fasta_path) as handle:
            d = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        return {k: str(v.seq) for k, v in d.items()}

    def import_as_alignment(self) -> Align.MultipleSeqAlignment:
        """return multiple sequence alignment object

        Returns
        -------
        Align.MultipleSeqAlignment
            multiple sequence alignment object
        """
        with open(self.fasta_path) as handle:
            return AlignIO.read(handle, "fasta")


def strip_dashes_from_str(seq_str):
    """remove `-` characters from a sequence string

    Parameters
    ----------
    seq_str : str
        sequence string

    Returns
    -------
    str
        sequence string with `-` characters removed
    """
    return seq_str.replace("-", "")


def strip_dashes_from_sequences(sequences: list[SeqRecord]) -> list[SeqRecord]:
    """remove `-` characters from sequences in a list of SeqRecord objects

    Parameters
    ----------
    sequences : list
        list of SeqRecord objects

    Returns
    -------
    list
        list of SeqRecord objects with `-` characters removed from sequences
    """
    sequences = copy.deepcopy(sequences)
    for s in sequences:
        s.seq = Seq.Seq(strip_dashes_from_str(str(s.seq)))
    return sequences


def write_fasta(sequences, output_filename):
    """save list of SeqRecord objects to a fasta file

    Parameters
    ----------
    sequences : list
        list of SeqRecord objects
    output_filename : str
        file path to new output fasta file
    """
    with open(output_filename, "w") as handle:
        SeqIO.write(sequences, handle, "fasta")


def longest_common_substring(s1, s2):
    """find the largest substring shared by two strings"""
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest : x_longest]


# ==============================================================================
# // alignment tools
# ==============================================================================


def reindex_alignment(record):
    """convert indexes of an ungapped sequence to indexes of it's sequence with gaps

    Parameters
    ----------
    record : SeqRecord object
        biopython sequence record object with gaps present

    Returns
    -------
    str
        ungapped sequence
    list
        list of indexes of nongap characters in the gapped sequence. index the return list with
        non-gap indexes to get the index of the gapped sequence. If the gapped sequence is `A--A--A`, the return list will be `[0, 3, 6]`
    """
    s = str(record.seq)
    unal_seq = ""
    index_map = []
    for al_pos, i in enumerate(s):
        # print(al_pos, i)
        if i != "-":
            unal_seq = unal_seq + i
            index_map.append(al_pos)
    return unal_seq, index_map


def reindex_alignment_str(seq_str):
    """convert indexes of an ungapped sequence to indexes of it's sequence with gaps

    Parameters
    ----------
    seq_str : str
        string of a sequence with gaps present

    Returns
    -------
    str
        ungapped sequence
    list
        list of indexes of nongap characters in the gapped sequence. index the return list with
        non-gap indexes to get the index of the gapped sequence. If the gapped sequence is `A--A--A`, the return list will be `[0, 3, 6]`

    Examples
    --------
    >>> aligned = 'A--A--A'
    >>> unaligned, ind = reindex_alignment_str(aligned)
    >>> print(unaligned)
    'AAA'
    >>> print(ind)
    [0, 3, 6]
    >>> print(aligned[ind[0]:ind[-1]+1])
    'A--A--A'
    """
    unal_seq = ""
    index_map = []
    for al_pos, i in enumerate(seq_str):
        # print(al_pos, i)
        if i != "-":
            unal_seq = unal_seq + i
            index_map.append(al_pos)
    return unal_seq, index_map


def find_alnslice_positions_in_unaln(aln: str, aln_start: int, aln_end: int):
    """Go from a slice of an aligned sequence, to the positions in the unaligned
    sequence that correspond to the slice.

    Parameters
    ----------
    aln : str
        aligned sequence with gaps
    aln_start : int
        start position of the slice of the aligned sequence
    aln_end : int
        end position of the slice of the aligned sequence

    Returns
    -------
    list
        list of positions of the slice in the unaligned sequence
    """
    positions = []
    unaln_index = 0
    for i in range(len(aln)):
        if aln[i] != "-":
            if aln_start <= i < aln_end:
                positions.append(unaln_index)
            unaln_index += 1
    return positions


def aln_2_idr_position_map(
    aln: MultipleSeqAlignment, idr_aln_start: int, idr_aln_end: int
) -> dict[str, list[int]]:
    """from an input multiple sequence alignment, generate a dictionary of the
    start and end positions of the IDRs in the unaligned sequences. The IDRs
    are defined by single start/end positions in the alignment. In other words,
    the IDRs are defined by a single slice of the alignment.

    Parameters
    ----------
    aln : MultipleSeqAlignment
        BioPython MultipleSeqAlignment object
    idr_aln_start : int
        start position of the IDRs in the alignment
    idr_aln_end : int
        end position of the IDRs in the alignment

    Returns
    -------
    dict[str, list[int]]
        idr_position_map: dictionary of the start and end positions of the IDRs
        in the unaligned sequences. The keys are the sequence ids and the values
        are lists of the start and end positions of the IDRs in the unaligned
        sequences. The start and end positions are 0-indexed.
    """
    idr_position_map = {}
    for seqrecord in list(aln):
        p = find_alnslice_positions_in_unaln(
            str(seqrecord.seq), idr_aln_start, idr_aln_end
        )
        if len(p) == 0:
            idr_position_map[seqrecord.id] = [-1, -1]
            continue
        idr_position_map[seqrecord.id] = [p[0], p[-1]]
    return idr_position_map


# t = '--LPPPP---PP------TEST--'
# tun = tools.strip_dashes_from_str(t)
# start = 6
# end = 13
# print(t[start:end])
# unal_pos =find_alnslice_positions_in_unaln(t, start, end)
# print(unal_pos)
# if len(unal_pos) == 0:
#     print('No positions')
# else:
#     print(tun[unal_pos[0]:unal_pos[-1]+1])
# for i in unal_pos:
#     print(i, tun[i])


def find_subseq_in_aligned_seq_str(unaligned_subseq_str, aligned_seq_str):
    """
    Find the position of unaligned_subsequence in aligned_seq_str
    return the positions of the subsequence in the aligned sequence indeces
    returns:
        al_start: start position of the subsequence in the aligned sequence
        al_end: end position of the subsequence in the aligned sequence
        aligned_seq_str[al_start:al_end]: the subsequence in the aligned sequence
    """
    unaligned_seq_str, alignment_index_key = reindex_alignment_str(aligned_seq_str)
    if unaligned_subseq_str not in unaligned_seq_str:
        raise ValueError(f"{unaligned_subseq_str} not in `aligned_seq_str`")
    unal_start = unaligned_seq_str.find(unaligned_subseq_str)
    unal_end = unal_start + len(unaligned_subseq_str) - 1
    al_start = alignment_index_key[unal_start]
    al_end = alignment_index_key[unal_end]
    return al_start, al_end, aligned_seq_str[al_start : al_end + 1]


# ==============================================================================
# // parse accession ids
# ==============================================================================


# def split_uniprot(prot_id):
#     j = re.compile(r"^[st][pr]\|(.+)\|(.+)")
#     prot = j.findall(prot_id)[0]
#     accession = prot[0]
#     name = prot[1]
#     return name, accession


# def get_description(desc):
#     r = re.compile(r"^sp\|.+\|\w+ (.+) OS=.+")
#     return r.findall(desc)[0]


# def convert_id_to_name(uniprot_id):
#     import requests

#     url = f"https://www.uniprot.org/uniprot/{uniprot_id}.xml"
#     response = requests.get(url)

#     if response.status_code == 200:
#         content = response.text
#         start = content.find("<name>") + len("<name>")
#         end = content.find("</name>")
#         uniprot_name = content[start:end].strip()
#         return uniprot_name
#     else:
#         return "not_found"


# ==============================================================================
# // regex utilities
# ==============================================================================


def regex2overlapping(regex):
    new_regex = r"(?=(" + regex + "))"
    return new_regex


def get_regex_matches(regex_pattern: str, seq_str: str):
    """searches for all matches of a regex pattern in a sequence string
    returns a generator object that yields the match sequence, start index, and end index

    Parameters
    ----------
    regex_pattern : str
        regular expression pattern
    seq_str : str
        string to search for matches

    Yields
    ------
    tuple
        (match sequence, start index, end index)
    """
    p = re.compile(regex_pattern)
    for m in p.finditer(seq_str):
        if m.start() == m.end():
            # even if there are groups in the lookahead, the first group should be the full match b/c that group surrounds the entire regex
            # so this will work whether or not there are groups in the lookahead
            match_seq = m.groups()[0]
        else:
            match_seq = seq_str[m.start() : m.end()]
        yield match_seq, m.start(), m.start() + len(match_seq) - 1
