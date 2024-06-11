from pathlib import Path

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
from local_conservation_scores import PairwiseMatrixMethods
from local_seqtools import general_utils as tools
# import gc


PAIRWISEMATFUNCS = PairwiseMatrixMethods()


def pairwise_matrices(
    og: group_tools.ConserGene,
    levels: list[str],
    score_key: str,
    score_function_name: str,
    score_params: dict,
    lflank: int,
    rflank: int,
    output_folder: str | Path | None = None,
    # device = None,
    # threads = None,
    # EsmMod = None,
):
    score_function = PAIRWISEMATFUNCS.__getitem__(score_function_name)
    if output_folder is None:
        output_folder = Path(og.info_dict["analysis_folder"])
    output_folder = Path(output_folder)
    output_folder.mkdir(exist_ok=True, parents=True)
    # get idr sequence
    query_idr = og.query_sequence[og.idr_start : og.idr_end + 1]
    # get hit positions relative to the idr
    hit_st_idr = og.hit_start_position - og.idr_start
    hit_end_idr = og.hit_end_position - og.idr_start
    flanked_hit_st_idr, flanked_hit_end_idr, flanked_hit = tools.pad_hit(
        query_idr, hit_st_idr, hit_end_idr, lflank, rflank
    )
    orig_hit_st_in_flanked_hit = hit_st_idr - flanked_hit_st_idr
    orig_hit_end_in_flanked_hit = hit_end_idr - flanked_hit_st_idr
    k = len(flanked_hit)
    for level in levels:
        if level not in og.level_objects:
            continue
        lvlo = og.level_objects[level]
        aln_file = lvlo.alignment_file
        if score_key in og.info_dict["orthogroups"][level]["conservation_scores"]:
            print(f"WARNING: {score_key} already exists in {og.reference_index}-{level}")
            print(f"skipping {score_key} for {og.reference_index}-{level}")
            continue
        # output_file = output_folder / f'{og.reference_index}-{aln_file.stem}-{score_key}.json'
        output_file = output_folder / f'{og.reference_index}-{level}-{score_key}.json'
        score_function(
            input_alignment_file=aln_file,
            output_file=output_file,
            reference_id=og.query_gene_id,
            k=k,
            idr_aln_st=lvlo.idr_aln_start,
            idr_aln_end=lvlo.idr_aln_end,
            **score_params,
        )
        # if "conservation_scores" not in og.info_dict["orthogroups"][level]: # this should've been taken care of in step 1
        #     og.info_dict["orthogroups"][level]["conservation_scores"] = {}
        og.info_dict["orthogroups"][level]["conservation_scores"][
            f"{score_key}"
        ] = {}
        score_dict = og.info_dict["orthogroups"][level]["conservation_scores"][
            f"{score_key}"
        ]
        score_dict["matrix_file"] = str(output_file)
        score_dict["flanked_hit"] = flanked_hit
        score_dict["flanked_hit_start_position_in_idr"] = flanked_hit_st_idr
        score_dict["original_hit_st_in_flanked_hit"] = orig_hit_st_in_flanked_hit
        score_dict["original_hit_end_in_flanked_hit"] = orig_hit_end_in_flanked_hit
        score_dict["score_function_name"] = score_function_name
        score_dict["score_params"] = {k:v for k,v in score_params.items() if k != "mod"}
        score_dict["lflank"] = lflank
        score_dict["rflank"] = rflank
        og._overwrite_json()
    # gc.collect()


def compute_pairwise_matrices(
    json_file: str | Path,
    score_key: str,
    score_function_name: str,
    score_params: dict,
    level: str|None = None,
    lflank: int = 0,
    rflank: int = 0,
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
    if hasattr(PAIRWISEMATFUNCS, score_function_name):
        pairwise_matrices(
            og=og, 
            levels=levels,
            score_key=score_key, 
            score_function_name=score_function_name,
            score_params=score_params, 
            lflank=lflank, 
            rflank=rflank,
            # device=device,
            # threads=threads,
            # EsmMod=EsmMod,
        )
