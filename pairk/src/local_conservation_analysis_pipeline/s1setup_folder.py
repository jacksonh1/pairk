import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import local_env_variables.env_variables as env
import local_seqtools.general_utils as tools


def import_and_reindex_hits_df(hits_file):
    hits_df = pd.read_csv(hits_file)
    assert (
        "reference_index" not in hits_df.columns
    ), "hits_df already has a reference_index column, I need to use that name for reindexing. Please rename the column or set `new_index`=False."
    hits_df = hits_df.reset_index(names="reference_index")
    return hits_df


def aln_2_query_seq(aln_file: str | Path, query_id: str):
    fasta_importer = tools.FastaImporter(aln_file)
    aln_dict = fasta_importer.import_as_dict()
    query_seq = aln_dict[query_id]
    query_seq_stripped = tools.strip_dashes_from_str(str(query_seq.seq))
    return query_seq_stripped


def ortholog_analysis_setup(
    table_row: pd.Series,
    database_key: dict,
):
    output_dict = {}
    query_gene_id = table_row["gene_id"]
    ref_ind = int(table_row["reference_index"])
    output_dict["reference_index"] = ref_ind
    output_dict["query_gene_id"] = query_gene_id
    try:
        database_files = database_key[query_gene_id]
    except KeyError as ke:
        output_dict[
            "critical_error"
        ] = f"not found in database: {ref_ind} - {query_gene_id}: {ke}"
        return output_dict
    levels = [k for k in database_files.keys() if k != "query_uniprot_id"]
    if "query_uniprot_id" in database_files:
        query_uniprotid = database_files["query_uniprot_id"]
    else:
        query_uniprotid = None
    output_dict["query_uniprot_id"] = query_uniprotid
    output_dict["orthogroups"] = {}
    for level in levels:
        output_dict["orthogroups"][level] = copy.deepcopy(database_files[level])
        if "conservation_scores" not in output_dict["orthogroups"][level]:
            output_dict["orthogroups"][level]["conservation_scores"] = {}
    output_dict["query_sequence"] = aln_2_query_seq(
        output_dict["orthogroups"][levels[0]]["alignment_file"],
        query_gene_id,
    )
    return output_dict


def ortholog_analysis_setup_search(
    table_row: pd.Series,
    database_key: dict,
):
    output_dict = ortholog_analysis_setup(table_row, database_key)
    if "critical_error" in output_dict:
        return output_dict
    # make sure that hit sequence is uppercase
    hit_sequence = table_row["hit_sequence"].upper()
    output_dict["hit_sequence"] = hit_sequence
    # if hit_sequence not in output_dict["query_sequence"]:
    #     raise ValueError(f"hit sequence not found in query sequence: {output_dict['reference_index']} - {hit_sequence} not in {output_dict['query_sequence']}")
    return output_dict


def ortholog_analysis_setup_given_positions(
    table_row: pd.Series,
    database_key: dict,
):
    output_dict = ortholog_analysis_setup(table_row, database_key)
    if "critical_error" in output_dict:
        return output_dict
    hit_start_pos = int(table_row["hit start position"])
    hit_end_pos = int(table_row["hit end position"])
    output_dict["hit_start_position"] = hit_start_pos
    output_dict["hit_end_position"] = hit_end_pos
    output_dict["hit_sequence"] = output_dict["query_sequence"][
        hit_start_pos : hit_end_pos + 1
    ]

    if "hit_sequence" in table_row:
        if output_dict["hit_sequence"] != table_row["hit_sequence"]:
            raise ValueError(
                f"hit sequence from table does not match extracted sequence: {output_dict['reference_index']} - {output_dict['hit_sequence']} != {table_row['hit_sequence']}"
            )
    return output_dict


# def main(hits_file, database_key_file):
#     hits_df = import_and_reindex_hits_df(hits_file)
#     with open(database_key_file, 'r') as f:
#         database_key = json.load(f)
#     output_dicts = []
#     for i, row in hits_df.iterrows():
#         output_dict = ortholog_analysis_setup(row, database_key)
#         output_dicts.append(output_dict)
#     output_file = hits_file.replace(".csv", "_odb_search_setup.json")
#     with open(output_file, 'w') as f:
#         json.dump(output_dicts, f, indent=4)
#     return output_dicts


def main(hits_file, database_key_file, output_folder, hit_search_method, new_index):
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    if hit_search_method not in ["search", "given_positions"]:
        raise ValueError(
            f'hit_search_method must be "search" or "given_positions", not {hit_search_method}'
        )

    if new_index:
        hits_df = import_and_reindex_hits_df(hits_file)
        hits_df.to_csv(output_folder / Path(hits_file).name.replace(".csv", "_original_reindexed.csv"), index=False)
    else:
        hits_df = pd.read_csv(hits_file)
        assert (
            "reference_index" in hits_df.columns
        ), "hits_df should already have a reference_index column if new_index=False."

    with open(database_key_file, "r") as f:
        database_key = json.load(f)

    for i, row in hits_df.iterrows():
        if hit_search_method == "search":
            output_dict = ortholog_analysis_setup_search(row, database_key)
        else:
            output_dict = ortholog_analysis_setup_given_positions(row, database_key)
        # if 'critical_error' in output_dict:
        analysis_folder = output_folder / f'{output_dict["reference_index"]}-{str(output_dict["query_gene_id"]).replace(":","_")}'
        analysis_folder.mkdir(parents=True, exist_ok=True)
        output_dict["analysis_folder"] = str(analysis_folder.resolve())
        output_file = (
            analysis_folder / f'{output_dict["reference_index"]}-{str(output_dict["query_gene_id"]).replace(":","_")}.json'
        )
        if output_file.exists():
            print(f"File already exists: {output_file}")
            print("Skipping...")
            print("clear the folder and rerun if you want to re-do the analysis")
            continue
        with open(output_file, "w") as f:
            json.dump(output_dict, f, indent=4)


# if __name__ == "__main__":
#     main(HITS_FILE, DATABASE_KEY_FILE, OUTPUT_FOLDER, HIT_SEARCH_METHOD)

