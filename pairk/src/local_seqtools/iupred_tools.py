import sys

import local_env_variables.env_variables as env
import local_seqtools.general_utils as tools

sys.path.append(str(env.IUPRED_DIR))
import iupred2a_lib as iup


def iupred_scores_from_seqstring(seqstring):
    return iup.iupred(seqstring)[0]


def get_idr_indexes(score_list, cutoff):
    idr_indexes = []
    in_idr = False
    st, end = 0, 0
    for c, i in enumerate(score_list):
        # print(c,i)
        if i > cutoff:
            if in_idr:
                end = c
            else:
                st = c
                end = c
                in_idr = True
        else:
            if in_idr:
                if st != end:
                    idr_indexes.append([st, end])
                in_idr = False
            else:
                # print(f'{c}-pass')
                pass
    if in_idr:
        idr_indexes.append([st, end])
    return idr_indexes


def merge_idr_regions(idr_regions, gap_merge_threshold=10):
    """
    if number of residues between two consecutive IDRs is less than or equal to this, then merge them
    """
    small_gap_merged_idr_regions = []
    if len(idr_regions) <= 1:
        return idr_regions
    region = idr_regions[0]
    for i in range(len(idr_regions) - 1):
        gapi = idr_regions[i + 1][0] - region[1]
        remaining = len(idr_regions) - i - 2
        if gapi <= gap_merge_threshold:
            region = [region[0], idr_regions[i + 1][1]]
            if remaining == 0:
                small_gap_merged_idr_regions.append(region)
        else:
            small_gap_merged_idr_regions.append(region)
            region = idr_regions[i + 1]
            if remaining == 0:
                small_gap_merged_idr_regions.append(region)
    return small_gap_merged_idr_regions


def len_filter_idr_regions(idr_regions, idr_min_length=10):
    new_idr_regions = []
    for indexes in idr_regions:
        # print(indexes[1]-indexes[0])
        if (indexes[1] + 1) - indexes[0] >= idr_min_length:
            new_idr_regions.append(indexes)
    return new_idr_regions


def main_find_idr_regions(
    sequence_str: str,
    iupred_cutoff: float = 0.4,
    gap_merge_threshold: int = 10,
    idr_min_length: int = 8,
) -> list[list[int]]:
    """finds the IDRs in a sequence string. IDRs are regions with IUPred scores above a cutoff. IDRs are merged if they are separated by a gap of less than or equal to a threshold. IDRs are filtered if they are shorter than a minimum length.

    Parameters
    ----------
    sequence_str : str
        amino acid sequence in string format
    iupred_cutoff : float, optional
        iupred score cutoff. Scores above this value are considered disordered, by default 0.4
    gap_merge_threshold : int, optional
        If there are `gap_merge_threshold` or fewer residues in between 2 idrs
        (regions with iupred score > `iupred_cutoff`), then the regions are
        merged into 1 region. For example [0, 10], [12, 58] -> [0, 58].
        by default 10
    idr_min_length : int, optional
        disordered regions that are shorter than this value are removed, by default 8

    Returns
    -------
    list[list[int]]]
        list of lists, each sublist is a pair of indexes that represent the start and end of an IDR
        The end index is not inclusive
    """
    idr_regions = get_idr_indexes(
        iupred_scores_from_seqstring(sequence_str), cutoff=iupred_cutoff
    )
    idr_regions = merge_idr_regions(
        idr_regions, gap_merge_threshold=gap_merge_threshold
    )
    idr_regions = len_filter_idr_regions(idr_regions, idr_min_length=idr_min_length)
    return idr_regions
