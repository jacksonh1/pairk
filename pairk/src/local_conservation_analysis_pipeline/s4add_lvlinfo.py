import json
from pathlib import Path

from Bio import AlignIO

import local_env_variables.env_variables as env
import local_seqtools.general_utils as tools


def get_info_from_aln_str(hit_start, hit_end, idr_start, idr_end, aln_str):
    _, aln_index = tools.reindex_alignment_str(aln_str)
    hit_aln_start = aln_index[hit_start]
    hit_aln_end = aln_index[hit_end]
    idr_aln_start = aln_index[idr_start]
    idr_aln_end = aln_index[idr_end]
    return hit_aln_start, hit_aln_end, idr_aln_start, idr_aln_end


def main(
    json_file,
    min_num_orthologs,
):
    with open(json_file, "r") as f:
        info_dict = json.load(f)
    hit_start = info_dict["hit_start_position"]
    hit_end = info_dict["hit_end_position"]
    idr_start = info_dict["idr_start"]
    idr_end = info_dict["idr_end"]
    passing_levels = []
    for lvl in info_dict["orthogroups"].keys():
        lvl_dict = info_dict["orthogroups"][lvl]
        fa_importer = tools.FastaImporter(lvl_dict["alignment_file"])
        aln = fa_importer.import_as_dict()
        query_aln_str = str(aln[info_dict["query_gene_id"]].seq)
        hit_aln_start, hit_aln_end, idr_aln_start, idr_aln_end = get_info_from_aln_str(
            hit_start, hit_end, idr_start, idr_end, query_aln_str
        )
        lvl_dict["hit_aln_start"] = hit_aln_start
        lvl_dict["hit_aln_end"] = hit_aln_end
        lvl_dict["idr_aln_start"] = idr_aln_start
        lvl_dict["idr_aln_end"] = idr_aln_end
        lvl_dict["query_aln_sequence"] = query_aln_str
        lvl_dict["hit_aln_sequence"] = query_aln_str[hit_aln_start : hit_aln_end + 1]
        lvl_dict["num_clustered_ldos"] = len(aln)
        if lvl_dict["num_clustered_ldos"] >= min_num_orthologs:
            passing_levels.append(lvl)
    # sort passing levels by the ordering in another list
    if all([x in env.PHYLOGENY_LVL_ORDERING for x in passing_levels]):
        passing_levels = sorted(
            passing_levels,
            key=lambda x: env.PHYLOGENY_LVL_ORDERING.index(x),
        )
    else:
        passing_levels = passing_levels
    
    info_dict["levels_passing_filters"] = passing_levels
    with open(json_file, "w") as f:
        json.dump(info_dict, f, indent=4)

