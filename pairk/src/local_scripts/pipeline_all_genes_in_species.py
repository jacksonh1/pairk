import argparse
import multiprocessing
import shutil
import traceback
from pathlib import Path

import local_config.orthodb_pipeline_parameters as conf
import local_orthoDB_group_pipeline.sql_queries as sql_queries
# import local_scripts.create_filemap as create_filemap
import local_scripts.odb_group_pipeline as pipeline

SPECIES_ID = "9606_0"
N_CORES = multiprocessing.cpu_count() - 2


def multiple_levels(
    config: conf.PipelineParams, query_odb_gene_id: str, og_levels: list
):
    """
    run the pipeline for a single odb_gene_id for multiple og_levels
    """
    for og_level in og_levels:
        config.og_select_params.OG_level_name = og_level
        try:
            pipeline.main_pipeline(config, odb_gene_id=query_odb_gene_id)
        except ValueError as err:
            traceback.print_exc()
            # logger.error(f"{query_geneid} - {og_level} - {err}")
            print(f"{query_odb_gene_id} - {og_level} - {err}")


def main(
    config: conf.PipelineParams,
    og_levels: list,
    multiprocess=True,
    species_id=SPECIES_ID,
    n_cores=N_CORES,
    overwrite=False,
):
    odbgeneid_list = sql_queries.get_all_odb_gene_ids_from_species_id(species_id)
    if Path(config.main_output_folder).exists():
        if overwrite:
            shutil.rmtree(config.main_output_folder)
        else:
            raise FileExistsError(
                f"main_output_folder already exists: {config.main_output_folder}. Use -o flag to overwrite"
            )
    if multiprocess:
        p = multiprocessing.Pool(n_cores)
        f_args = [(config, i, og_levels) for i in odbgeneid_list]
        p.starmap(multiple_levels, f_args)
        p.close()
        p.join()
    else:
        for i in odbgeneid_list:
            multiple_levels(config, i, og_levels)


if __name__ == "__main__":
    og_levels = ["Eukaryota", "Mammalia", "Metazoa", "Tetrapoda", "Vertebrata"]
    parser = argparse.ArgumentParser(
        description="run the pipeline for all genes in an organism",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        metavar="<file>",
        default=None,
        help="""path to config file""",
    )
    parser.add_argument(
        "-n",
        "--n_cores",
        type=int,
        metavar="<int>",
        default=N_CORES,
        help=f"""number of cores to use""",
    )
    parser.add_argument(
        "-s",
        "--species_id",
        type=str,
        metavar="<str>",
        default=SPECIES_ID,
        help=f"""species id to use""",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        action="store_true",
        help="""if flag is provided and the main_output_folder exists, it will be removed and overwritten by the new files. Otherwise, an error will be raised if the folder exists""",
    )
    args = parser.parse_args()
    config = pipeline.load_config(args.config)
    main(
        config,
        og_levels,
        multiprocess=True,
        species_id=args.species_id,
        n_cores=args.n_cores,
        overwrite=args.overwrite,
    )
    # create_filemap.create_filemap(
    #     config.main_output_folder,
    #     output_file=Path(config.main_output_folder) / "filemap.json",
    # )
