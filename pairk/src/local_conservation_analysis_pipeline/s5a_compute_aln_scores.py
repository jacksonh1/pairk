from pathlib import Path

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
from local_conservation_scores import ConservationScoreMethods, PairwiseMatrixMethods
import local_seqtools.general_utils as tools
import numpy as np
# import gc


SCORES = ConservationScoreMethods()
PAIRWISESCORES = PairwiseMatrixMethods()


def get_hit_aln_scores(lvlo: group_tools.LevelAlnScore):
    hit_slice = slice(lvlo.hit_aln_start, lvlo.hit_aln_end + 1)
    hit_scores = lvlo.scores[hit_slice]
    hit_z_scores = lvlo.z_scores[hit_slice]
    hit_aln_seq = lvlo.query_aln_sequence[hit_slice]
    nongap_inds = tools.get_non_gap_indexes(hit_aln_seq)
    return list(np.array(hit_scores)[nongap_inds]), list(np.array(hit_z_scores)[nongap_inds])

# def get_hit_aln_scores(lvlo: group_tools.LevelAlnScore):
#     """
#     returns a list of the scores and a list of the z scores for the hit (non-gap) positions in query sequence
#     """
#     hit_slice = slice(lvlo.hit_aln_start, lvlo.hit_aln_end + 1)
#     hit_scores = lvlo.scores[hit_slice]
#     hit_aln_seq = lvlo.query_aln_sequence[hit_slice]
#     nongap_inds = tools.get_non_gap_indexes(hit_aln_seq)
#     return list(np.array(hit_scores)[nongap_inds])


def alignment_scores(
    og: group_tools.ConserGene,
    levels: list[str],
    score_key: str,
    score_function_name: str,
    score_params: dict,
    output_folder: str|Path|None = None,
):
    score_function = SCORES.__getitem__(score_function_name)
    if output_folder is None:
        output_folder = Path(og.info_dict["analysis_folder"])
    output_folder = Path(output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)
    for level in levels:
        if level not in og.level_objects:
            continue
        lvlo = og.level_objects[level]
        aln_file = lvlo.alignment_file
        # output_file = output_folder / f"{aln_file.stem}-{score_key}.json"
        output_file = output_folder / f"{og.reference_index}-{level}-{score_key}.json"
        if score_key in og.info_dict["orthogroups"][level]["conservation_scores"]:
            print(f"WARNING: {score_key} already exists in {og.reference_index}-{level}")
            print(f"skipping {score_key} for {og.reference_index}-{level}")
            continue
        score_function(
            input_alignment_file=aln_file,
            output_file=output_file,
            reference_id=og.query_gene_id,
            **score_params,
        )
        # if "conservation_scores" not in og.info_dict["orthogroups"][level]:
        #     og.info_dict["orthogroups"][level]["conservation_scores"] = {}
        og.info_dict["orthogroups"][level]["conservation_scores"][f"{score_key}"] = {}
        score_dict = og.info_dict["orthogroups"][level]["conservation_scores"][f"{score_key}"]
        score_dict["file"] = str(output_file)
        score_dict["score_function_name"] = score_function_name
        score_dict["score_params"] = score_params
        lvl_aln_o = group_tools.LevelAlnScore.from_conser_level(lvlo, score_key)
        if lvl_aln_o.z_score_failure is None:
            hit_scores, hit_z_scores = get_hit_aln_scores(lvl_aln_o)
            score_dict["hit_scores"] = hit_scores
            score_dict["hit_z_scores"] = hit_z_scores
        og._overwrite_json()
    # gc.collect()


def compute_aln_cons_scores(
    json_file: str | Path,
    score_key: str,
    score_function_name: str,
    score_params: dict,
    level: str|None = None,
    **kwargs,
    # device = None,
    # threads = None,
    # EsmMod = None,
):
    og = group_tools.ConserGene(json_file)
    if hasattr(og, "critical_error"):
        return
    og.load_levels()
    if level is not None:
        levels = [level]
    else:
        levels = list(og.level_objects.keys())
    if hasattr(SCORES, score_function_name):
        alignment_scores(
            og=og, 
            levels=levels,
            score_key=score_key, 
            score_function_name=score_function_name,
            score_params=score_params, 
        )
    # gc.collect()
