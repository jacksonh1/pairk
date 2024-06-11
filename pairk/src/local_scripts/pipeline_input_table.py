import argparse
import multiprocessing
import shutil
import traceback
from pathlib import Path
from typing import Literal

import pandas as pd

import local_config.orthodb_pipeline_parameters as conf
import local_orthoDB_group_pipeline.sql_queries as sql_queries
# import local_scripts.create_filemap as create_filemap
import local_scripts.odb_group_pipeline as pipeline

N_CORES = multiprocessing.cpu_count() - 2
OG_LEVELS = ["Eukaryota", "Mammalia", "Metazoa", "Tetrapoda", "Vertebrata"]


def multiple_levels(
    config: conf.PipelineParams,
    gene_id: str,
    og_levels: list,
    id_type: Literal["odb_gene_id", "uniprot_id"],
):
    """
    run the pipeline for a single odb_gene_id for multiple og_levels
    """
    assert id_type in ["odb_gene_id", "uniprot_id"], f"id_type must be 'odb_gene_id' or 'uniprot_id', not {id_type}"
    for og_level in og_levels:
        config.og_select_params.OG_level_name = og_level
        if id_type == "odb_gene_id":
            try:
                pipeline.main_pipeline(config, odb_gene_id=gene_id)
            except ValueError as err:
                traceback.print_exc()
                # logger.error(f"{query_geneid} - {og_level} - {err}")
                print(f"{gene_id} - {og_level} - {err}")
                continue
            except FileExistsError as err:
                print(f"{gene_id} - {og_level} - {err}")
                print(f"skipping {gene_id} - {og_level}")
                continue
        else:
            try:
                pipeline.main_pipeline(config, uniprot_id=gene_id)
            except ValueError as err:
                traceback.print_exc()
                # logger.error(f"{query_geneid} - {og_level} - {err}")
                print(f"{gene_id} - {og_level} - {err}")
                continue
            except FileExistsError as err:
                print(f"{gene_id} - {og_level} - {err}")
                print(f"skipping {gene_id} - {og_level}")
                continue



def main(
    config: conf.PipelineParams,
    table_file: str,
    og_levels: list,
    odb_gene_id_column: str | None = None,
    uniprot_id_column: str | None = None,
    n_cores=N_CORES,
    clear_output_folder=False,
    multiprocess=True,
):
    table = pd.read_csv(table_file)
    if Path(config.main_output_folder).exists():
        if clear_output_folder:
            shutil.rmtree(config.main_output_folder)
        else:
            print("WARNING: OUTPUT FOLDER ALREADY EXISTS")
            print("any files that already exist will either be overwritten or skipped (depending on the overwrite flag in the config file)")
            print("set clear_output_folder=True to delete the folder and re-run the pipeline if you want to start fresh")

    if odb_gene_id_column is not None:
        table = table.dropna(subset=[odb_gene_id_column])
        id_list = list(table[odb_gene_id_column].unique())
        id_type = "odb_gene_id"
    elif uniprot_id_column is not None:
        table = table.dropna(subset=[uniprot_id_column])
        id_list = list(table[uniprot_id_column].unique())
        id_type = "uniprot_id"
    else:
        raise ValueError(
            "either odb_gene_id_column or uniprot_id_column must be provided"
        )
    if multiprocess:
        p = multiprocessing.Pool(n_cores)
        f_args = [(config, i, og_levels, id_type) for i in id_list]
        p.starmap(multiple_levels, f_args)
        p.close()
        p.join()
    else:
        for i in id_list:
            print(i)
            multiple_levels(config, i, og_levels, id_type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="run the pipeline for all odb_gene_ids in an input table",
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
        "-t",
        "--table",
        type=str,
        metavar="<file>",
        required=True,
        help=f"""input table file""",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--uniprot_id_column",
        type=str,
        metavar="<str>",
        help='the name of the column in the input table containing the uniprot ids of the genes of interest (e.g. "uniprot_id")',
    )
    group.add_argument(
        "--odb_gene_id_column",
        type=str,
        metavar="<str>",
        help="the name of the column in the input table containing the odb_gene_ids",
    )
    parser.add_argument(
        "--clear",
        action="store_true",
        help="""if flag is provided and the main_output_folder exists, it will be removed and overwritten by the new files""",
    )
    parser.add_argument(
        "-l",
        "--og_levels",
        nargs="*",
        metavar="<list>",
        default=OG_LEVELS,
        help=f"""list of orthologous group levels to run the pipeline for""",
    )
    args = parser.parse_args()
    # for arg in vars(args):
    #     print(f"{arg}: {getattr(args, arg)}")
    config = pipeline.load_config(args.config)
    main(
        config,
        table_file=args.table,
        og_levels=args.og_levels,
        odb_gene_id_column=args.odb_gene_id_column,
        uniprot_id_column=args.uniprot_id_column,
        n_cores=args.n_cores,
        clear_output_folder=args.clear,
        multiprocess=True,
    )
    # create_filemap.create_filemap(
    #     config.main_output_folder,
    #     output_file=Path(config.main_output_folder) / "filemap.json",
    # )
