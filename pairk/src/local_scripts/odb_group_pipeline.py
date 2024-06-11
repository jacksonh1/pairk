#!/usr/bin/env python

import argparse
import copy
import json
from pathlib import Path

import yaml
from attrs import asdict
from Bio import SeqIO

import local_env_variables.env_variables as env
import local_seqtools.cli_wrappers as cli_wrappers
from local_config import orthodb_pipeline_parameters
from local_orthoDB_group_pipeline import (cluster, filters, find_LDOs,
                                          og_selection, sql_queries,
                                          uniprotid_search)

ODB_DATABASE = env.orthoDB_database()

def load_config(config_file: str | None) -> orthodb_pipeline_parameters.PipelineParams:
    if config_file is None:
        config = orthodb_pipeline_parameters.PipelineParams()
    else:
        with open(config_file, 'r') as f:
            config_dict = yaml.safe_load(f)
        config = orthodb_pipeline_parameters.PipelineParams.from_dict(config_dict)
    return config

def generate_species_map(odb_gene_id_list: list[str]):
    species_map = {}
    for odb_gene_id in odb_gene_id_list:
        species_id = sql_queries.odb_gene_id_2_species_id(odb_gene_id)
        species_map[odb_gene_id] = ODB_DATABASE.data_species_dict[species_id]
    return species_map

def filter_sequences(min_fraction_shorter_than_query, query_seqrecord, sequence_dict):
    filtered_sequence_dict = filters.filter_seqs_with_nonaa_chars(
        sequence_dict,
    )
    min_length = min_fraction_shorter_than_query* len(query_seqrecord)
    filtered_sequence_dict = filters.filter_shorter_sequences(
        filtered_sequence_dict,
        min_length=min_length,
    )
    if query_seqrecord.id not in filtered_sequence_dict:
        filtered_sequence_dict[query_seqrecord.id] = copy.deepcopy(query_seqrecord)
    return filtered_sequence_dict


def save_info_json(output_dict: dict, output_file: str|Path):
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, 'w') as f:
        json.dump(output_dict, f, indent=4)


def _pipeline(config: orthodb_pipeline_parameters.PipelineParams, odb_gene_id: str) -> dict:
    """runs the pipeline for a odb_gene_id. This isn't meant to be called directly,
    Instead, use pipeline_from_uniprot_id or pipeline_from_odb_gene_id.

    Parameters
    ----------
    config : conf.PipelineParams
        pipeline parameters 
    odb_gene_id : str
        the orthoDB gene id of the gene of interest (e.g. "9606_0:001c7b")

    Returns
    -------
    dict
        the results of the pipeline in a dictionary
    """    
    results_dict = {}
    try:
        ogid, oglevel = og_selection.select_OG_by_level_name(
            odb_gene_id=odb_gene_id,
            level_name=config.og_select_params.OG_level_name
        )
    except ValueError as e:
        results_dict['critical error'] = str(e)
        return results_dict
    
    group_members = sql_queries.ogid_2_odb_gene_id_list(ogid)
    sequence_dict = ODB_DATABASE.get_sequences_from_list_of_seq_ids(group_members)
    query_seqrecord = sequence_dict[odb_gene_id]

    filtered_sequence_dict = filter_sequences(
        config.filter_params.min_fraction_shorter_than_query,
        query_seqrecord,
        sequence_dict,
    )

    pid_df, ldos = find_LDOs.find_LDOs_main(
        seqrecord_dict = filtered_sequence_dict,
        query_seqrecord = query_seqrecord,
        pid_method = config.ldo_select_params.LDO_selection_method,
        n_align_threads = config.ldo_select_params.LDO_mafft_threads,
        mafft_executable = config.ldo_select_params._LDO_mafft_exe,
        extra_args = config.ldo_select_params._LDO_mafft_additional_args,
    )
    ldo_seqrecord_dict = ODB_DATABASE.get_sequences_from_list_of_seq_ids(ldos)

    cdhit_command, clustered_ldo_seqrec_dict = cluster.cdhit_main(
        ldo_seqrecord_dict,
        odb_gene_id,
        cd_hit_executable=config._cd_hit_exe,
        extra_args=config._cd_hit_additional_args
    )

    results_dict['query_odb_gene_id'] = odb_gene_id
    results_dict['query_sequence_str'] = str(query_seqrecord.seq)
    results_dict['ogid'] = ogid
    results_dict['oglevel'] = oglevel
    results_dict['sequences'] = list(sequence_dict.keys())
    results_dict['sequences_filtered'] = list(filtered_sequence_dict.keys())
    results_dict['sequences_ldos'] = list(ldo_seqrecord_dict.keys())
    results_dict['sequences_clustered_ldos'] = clustered_ldo_seqrec_dict
    results_dict['cdhit_command'] = cdhit_command
    results_dict['species_map'] = generate_species_map(list(clustered_ldo_seqrec_dict.keys()))
    return results_dict


def pipeline_from_uniprot_id(config: orthodb_pipeline_parameters.PipelineParams, uniprot_id: str):
    try:
        odb_gene_id = uniprotid_search.uniprotid_2_odb_gene_id(uniprot_id)
    except ValueError as e:
        output_dict = {}
        output_dict['query_uniprot_id'] = uniprot_id
        output_dict['critical error'] = str(e)
        return output_dict
    output_dict = _pipeline(config, odb_gene_id)
    output_dict['query_uniprot_id'] = uniprot_id
    return output_dict


def pipeline_from_odb_gene_id(config: orthodb_pipeline_parameters.PipelineParams, odb_gene_id: str):
    query_uniprot_id = sql_queries.odb_gene_id_2_uniprotid(odb_gene_id)
    output_dict = _pipeline(config, odb_gene_id)
    output_dict['query_uniprot_id'] = query_uniprot_id
    return output_dict


def main_pipeline(config: orthodb_pipeline_parameters.PipelineParams, uniprot_id: str | None = None, odb_gene_id: str | None = None):
    """run the main pipeline for a single gene. Either uniprot_id or odb_gene_id must be provided

    Parameters
    ----------
    config : conf.PipelineParams
        pipeline parameters in a PipelineParams object
    uniprot_id : str | None, optional
        uniprot id of the query protein. If not provided, then `odb_gene_id` must be provided, by default None
    odb_gene_id : str | None, optional
        orthoDB gene id of the query protein. If not provided, then `uniprot_id` must be provided, by default None

    Returns
    -------
    dict
        the pipeline results in a dictionary

    Raises
    ------
    ValueError
        raises a ValueError if neither `uniprot_id` nor `odb_gene_id` are provided
    ValueError
        raises a ValueError if there is a "critical error" in the pipeline
        When the pipeline is run, errors are stored in the output dictionary under the key "critical error". This error is raised if it exists
    """    
    if odb_gene_id is not None:
        output_dict = pipeline_from_odb_gene_id(config, odb_gene_id) # type: ignore
    elif uniprot_id is not None:
        output_dict = pipeline_from_uniprot_id(config, uniprot_id) # type: ignore
    else:
        raise ValueError('either uniprot_id or odb_gene_id must be provided')
    
    output_dict['processing params'] = asdict(config)
    og_info_json_folder = Path(config.main_output_folder) / 'info_jsons'
    og_info_failure_folder = og_info_json_folder / 'failures'
    
    if 'critical error' in output_dict:
        if config.write_files:
            og_info_failure_folder.mkdir(parents=True, exist_ok=True)
            og_info_json_file = og_info_failure_folder / f'{uniprot_id}{odb_gene_id}_info.json'
            save_info_json(output_dict, og_info_json_file)
        raise ValueError(output_dict['critical error'])
    
    output_file_prefix = f'{output_dict["query_odb_gene_id"].replace(":", "_")}_{output_dict["oglevel"]}_{output_dict["ogid"]}'
    og_info_json_file = og_info_json_folder / f'{output_file_prefix}_info.json'

    if og_info_json_file.exists() and not config.overwrite:
        raise FileExistsError(f'info json file already exists: {og_info_json_file}')
    
    if config.align_params.align:
        mafft_command, aln = cli_wrappers.mafft_align_wrapper(
            list(output_dict['sequences_clustered_ldos'].values()),
            n_align_threads=config.align_params.n_align_threads,
            mafft_executable=config.align_params._mafft_exe,
            extra_args=config.align_params._mafft_additional_args,
            output_format = "list",
        )
        if config.write_files:
            alignment_folder = Path(config.main_output_folder) / 'alignments'
            alignment_folder.mkdir(parents=True, exist_ok=True)
            alignment_output_file = alignment_folder / f'{output_file_prefix}_clustered_ldos_aln.fasta'
            with open(alignment_output_file, 'w') as f:
                SeqIO.write(aln, f, 'fasta')
            output_dict['alignment_clustered_ldos_file'] = str(alignment_output_file.resolve())
            # output_dict['alignment_clustered_ldos_file_relative'] = str(alignment_output_file.resolve().relative_to(Path.cwd()))
            output_dict['alignment_clustered_ldos_command'] = mafft_command

    output_dict['sequences_clustered_ldos'] = list(output_dict['sequences_clustered_ldos'].keys())
    if config.write_files:
        save_info_json(output_dict, og_info_json_file)
    return og_info_json_file, output_dict


if __name__ == "__main__":
    # get the default parameters just to print them in the help message
    # this is a bit hacky but it works
    # filter out the private attributes (those that start with '_') because they
    # are more advanced and probably won't be used by most users
    d_params = ''
    for k,v in asdict(orthodb_pipeline_parameters.PipelineParams(), filter=lambda attr, value: not str(attr.name).startswith('_')).items():
        d_params += f'- {k}: {v}\n'

    parser = argparse.ArgumentParser(
        description=f'''run main orthoDB group generation pipeline for a single gene
    processing parameters should be provided in a config file. (-c/--config)
    if no config file is provided, default parameters will be used
    The default parameters are:
{d_params}''',
        formatter_class=argparse.RawTextHelpFormatter
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        '-unid',
        '--uniprot_id',
        type=str,
        metavar='<str>',
        help='the uniprot id of the gene of interest'
    )
    group.add_argument(
        '-odbid',
        '--odb_gene_id',
        type=str,
        metavar='<str>',
        help='the odb gene id of the gene of interest (e.g. "9606_0:001c7b")'
    )
    parser.add_argument(
        '-c',
        '--config',
        type=str,
        metavar='<file>',
        default=None,
        help='''path to config file, default=None'''
    )
    args = parser.parse_args()
    config = load_config(args.config)
    main_pipeline(config, args.uniprot_id, args.odb_gene_id)
