import json
from pathlib import Path

import local_seqtools.iupred_tools as iup_tools

FIND_IDRS = True
IDR_MAP_FILE = None
IUPRED_CUTOFF = 0.4
GAP_MERGE_THRESHOLD = 10
IDR_MIN_LENGTH = 8


def main(
    info_json_file, find_idrs=True, idr_map_file: str | Path | None = None, **kwargs
):
    # either find_idrs or idr_map_file must be provided, but not both
    if find_idrs and idr_map_file is not None:
        raise ValueError(
            "find_idrs=True and idr_map_file was provided, only one of these is allowed"
        )
    if not find_idrs and idr_map_file is None:
        raise ValueError("either find_idrs or idr_map_file must be provided")

    with open(info_json_file, "r") as f:
        info_dict = json.load(f)
    query_sequence = info_dict["query_sequence"]
    if find_idrs:
        idr_map_file = None
        idrs = iup_tools.main_find_idr_regions(
            sequence_str=query_sequence,
            **kwargs,
        )
        info_dict["idr_regions"] = idrs
    else:
        with open(idr_map_file, "r") as f:  # type: ignore
            idr_map = json.load(f)
        info_dict["idr_regions"] = idr_map[info_dict["query_gene_id"]]
    with open(info_json_file, "w") as f:
        json.dump(info_dict, f, indent=4)
