import json
from pathlib import Path

import local_seqtools.general_utils as tools


def pad_hit_sequence(st: int, end: int, query_sequence: str, target_hit_length=36):
    hit_sequence = query_sequence[st : end + 1]
    if len(hit_sequence) >= target_hit_length:
        return st, end, hit_sequence
    flank = round((target_hit_length - len(hit_sequence)) / 2)
    s_i_new = max(0, st - flank)
    s_e_new = min(len(query_sequence), end + flank)
    return s_i_new, s_e_new, query_sequence[s_i_new : s_e_new + 1]


# pad_hit_sequence(3, 6, 'abcdefghij', target_hit_length=6)
# print('abcdefghij'[3:6+1])


def hit_in_query_bool(hit_sequence: str, query_sequence: str):
    return hit_sequence in query_sequence


def find_all_substring_occurences(substring: str, string: str):
    substring_length = len(substring)
    occurences = []
    for i in range(len(string)):
        if string[i : i + substring_length] == substring:
            occurences.append((i, i + substring_length - 1))
    return occurences


def find_query_hit_sequence(subsequence: str, sequence: str):
    occurences = find_all_substring_occurences(subsequence, sequence)
    if len(occurences) == 0:
        raise ValueError(f"subsequence {subsequence} not found in query sequence")
    elif len(occurences) == 1:
        return occurences[0]
    else:
        raise ValueError(
            f"subsequence {subsequence} found in query sequence more than once"
        )


def save_out_json(info_dict, json_file):
    with open(json_file, "w") as f:
        json.dump(info_dict, f, indent=4)


def hit_in_idr(hit_start_position, hit_end_position, idr_regions):
    d = {}
    d["hit_in_idr"] = False
    for region in idr_regions:
        if hit_start_position >= region[0] and hit_end_position <= region[1]:
            d["hit_in_idr"] = True
            d["idr_start"] = region[0]
            d["idr_end"] = region[1]
    return d


def main(
    json_file,
    search_method="search",
    longest_common_subsequence=False,
    lcs_min_length=20,
    target_hit_length=0,
):
    with open(json_file, "r") as f:
        info_dict = json.load(f)
    query_sequence = info_dict["query_sequence"]
    hit_sequence = info_dict["hit_sequence"]
    if search_method == "search":
        if longest_common_subsequence:
            hit_sequence = tools.longest_common_substring(
                hit_sequence,
                query_sequence,
            )
            if len(hit_sequence) < lcs_min_length:
                info_dict[
                    "critical_error"
                ] = f"could not find a common substring of at least {lcs_min_length} characters between the hit sequence and the query sequence"
                save_out_json(info_dict, json_file)
                return "fail"

        if hit_sequence not in query_sequence:
            info_dict[
                "critical_error"
            ] = f"hit sequence {hit_sequence} not found in query sequence"
            save_out_json(info_dict, json_file)
            return "fail"
        try:
            hit_start_position, hit_end_position = find_query_hit_sequence(
                hit_sequence, query_sequence
            )
        except ValueError as e:
            info_dict["critical_error"] = str(e)
            save_out_json(info_dict, json_file)
            return "fail"

        if target_hit_length > 0:
            hit_start_position, hit_end_position, hit_sequence = pad_hit_sequence(
                hit_start_position, hit_end_position, query_sequence, target_hit_length
            )
        info_dict["hit_start_position"] = hit_start_position
        info_dict["hit_end_position"] = hit_end_position
        info_dict["hit_sequence"] = hit_sequence
    else:
        hit_start_position = info_dict["hit_start_position"]
        hit_end_position = info_dict["hit_end_position"]
    d = hit_in_idr(hit_start_position, hit_end_position, info_dict["idr_regions"])
    info_dict.update(d)
    if not info_dict["hit_in_idr"]:
        info_dict["critical_error"] = "hit sequence not in an idr"
    save_out_json(info_dict, json_file)
    if "critical_error" in info_dict:
        return "fail"
    else:
        return "pass"


# def driver(
#     query_sequence,
#     hit_sequence,
#     longest_common_subsequence=False,
#     lcs_min_length=20,
#     target_hit_length=0,
# ):
#     if longest_common_subsequence:
#         hit_sequence = tools.longest_common_substring(
#             hit_sequence,
#             query_sequence,
#         )
#         if len(hit_sequence) < lcs_min_length:
#             info_dict[
#                 "critical_error"
#             ] = f"could not find a common substring of at least {lcs_min_length} characters between the hit sequence and the query sequence"
#             return info_dict

#     if hit_sequence not in query_sequence:
#         info_dict[
#             "critical_error"
#         ] = f"hit sequence {hit_sequence} not found in query sequence {query_sequence}"
#         return info_dict
#     try:
#         hit_start_position, hit_end_position = find_query_hit_sequence(
#             hit_sequence, query_sequence
#         )
#     except ValueError as e:
#         info_dict["critical_error"] = str(e)
#         return info_dict

#     if target_hit_length > 0:
#         hit_start_position, hit_end_position, hit_sequence = pad_hit_sequence(
#             hit_start_position, hit_end_position, query_sequence, target_hit_length
#         )
#     info_dict["hit_start_position"] = hit_start_position
#     info_dict["hit_end_position"] = hit_end_position
#     info_dict["hit_sequence"] = hit_sequence
#     return info_dict


# def main(
#     json_file,
#     longest_common_subsequence=False,
#     lcs_min_length=20,
#     target_hit_length=0,
# ):
#     with open(json_file, "r") as f:
#         info_dict = json.load(f)
#     query_sequence = info_dict["query_sequence"]
#     hit_sequence = info_dict["hit_sequence"]
#     hit_dict = driver(
#         query_sequence,
#         hit_sequence,
#         longest_common_subsequence=longest_common_subsequence,
#         lcs_min_length=lcs_min_length,
#         target_hit_length=target_hit_length,
#     )
#     info_dict.update(hit_dict)
#     save_out_json(info_dict, json_file)
