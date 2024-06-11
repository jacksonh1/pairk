from pathlib import Path

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
from local_conservation_scores import PairwiseMatrixKmerScoreMethods
from local_conservation_scores import ColumnwiseScoreMethods
from local_conservation_scores import PairwiseMatrixMethods
# import gc
from typing import Callable
from local_conservation_scores.tools import capra_singh_2007_scores as cs
from attrs import asdict, define, field, validators
from local_config import conservation_pipeline_parameters as conf
import traceback


PAIRWISEMETHODS = PairwiseMatrixKmerScoreMethods()
COLSCOREMETHODS = ColumnwiseScoreMethods()
PAIRWISEMATFUNCS = PairwiseMatrixMethods()


@define
class PairwiseScoreResults:
    hit_sequence: str
    hit_scores: list[float]
    hit_z_scores: list[float]
    flank_hit_sequence: str | None = None
    flank_hit_scores: list[float] | None = None
    flank_hit_z_scores: list[float] | None = None
    background_scores: list[float] | None = None


def lvlo_2_pairwise_scores(
    lvlo: group_tools.ConserLevel,
    score_key: str,
    # mat2score_func: str = "matrix_json_2_pairwise_scores",
    params: conf.PairMatrixToScoreConf,
):
    mat_function = PAIRWISEMETHODS.__getitem__(params.matrix_to_score_function_name)
    col_function = COLSCOREMETHODS.__getitem__(params.columnwise_score_function_name)
    
    pairdict = lvlo.conservation_scores[score_key]
    flanked_hit_scores = mat_function(
        pairdict["matrix_file"],
        pairdict["flanked_hit_start_position_in_idr"],
        columnwise_score_func=col_function,
        reciprocal_best_match=params.reciprocal_best_match,
        bg_cutoff=params.bg_cutoff,
        similarity_threshold=params.similarity_threshold,
        bg_kmer_cutoff=params.bg_kmer_cutoff,
    )
    hit_slice = slice(
        pairdict["original_hit_st_in_flanked_hit"],
        pairdict["original_hit_end_in_flanked_hit"] + 1,
    )
    scores = PairwiseScoreResults(
        hit_sequence=flanked_hit_scores.hit_sequence[hit_slice],
        hit_scores=flanked_hit_scores.hit_scores[hit_slice],
        hit_z_scores=flanked_hit_scores.hit_z_scores[hit_slice],
        flank_hit_sequence=flanked_hit_scores.hit_sequence,
        flank_hit_scores=flanked_hit_scores.hit_scores,
        flank_hit_z_scores=flanked_hit_scores.hit_z_scores,
        # background_scores=flanked_hit_scores.background_scores,
    )
    return scores


def pairwise_scores(
    og: group_tools.ConserGene,
    levels: list[str],
    score_key: str,
    params: conf.PairMatrixToScoreConf
):
    for level in levels:
        if level not in og.level_objects:
            continue
        lvlo = og.level_objects[level]
        if score_key not in lvlo.conservation_scores:
            raise ValueError(f"score_key {score_key} not in lvlo.conservation_scores")
        score_dict = og.info_dict["orthogroups"][level]["conservation_scores"][
            f"{score_key}"
        ]
        # scores = lvlo_2_pairwise_scores(
        #     lvlo=lvlo,
        #     score_key=score_key,
        #     params=params,
        # )
        try:
            scores = lvlo_2_pairwise_scores(
                lvlo=lvlo,
                score_key=score_key,
                params=params,
            )
        except ValueError as e:
            score_dict["score_error"] = str(e)
            continue
        score_dict["flanked_hit_sequence"] = scores.flank_hit_sequence
        score_dict["flanked_hit_scores"] = scores.flank_hit_scores
        score_dict["flanked_hit_z_scores"] = scores.flank_hit_z_scores
        score_dict["hit_sequence"] = scores.hit_sequence
        score_dict["hit_scores"] = scores.hit_scores
        score_dict["hit_z_scores"] = scores.hit_z_scores
        score_dict["mat2score_params"] = asdict(params)
    og._overwrite_json()
    # gc.collect()


def compute_hit_conservation_scores(
    json_file: str | Path,
    scoremethod: conf.ScoreMethod|conf.ScoreMethodEmbedding,
    params: conf.PairMatrixToScoreConf,
):
    og = group_tools.ConserGene(json_file)
    if hasattr(og, "critical_error"):
        return
    og.load_levels()
    if scoremethod.level is not None:
        levels = [scoremethod.level]
    else:
        levels = list(og.level_objects.keys())
    if hasattr(PAIRWISEMATFUNCS, scoremethod.score_function_name):
        pairwise_scores(
            og=og, 
            levels=levels,
            score_key=scoremethod.score_key, 
            params=params,
        )
