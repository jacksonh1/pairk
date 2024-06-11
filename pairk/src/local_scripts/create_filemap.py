'''
run this after the pipeline to create a json file that maps the gene ids to the files output by the pipeline

will create a file looking like:
{
    "odb_gene_id": {
        "Eukaryota": {
            "alignment_file": "path/to/file",
            "info_file": "path/to/file",
        },
        "Vertebrata": {
            "alignment_file": "path/to/file",
            "info_file": "path/to/file",
        },
        ...
    },
    ...
}
'''

import argparse
import json
from pathlib import Path


def get_json_file_list(json_dir: str|Path) -> list[Path]:
    json_dir = Path(json_dir)
    json_files = list(json_dir.glob('*.json'))
    return json_files


def get_map_info_from_json(json_file: str|Path) -> tuple[str, str, str]:
    with open(json_file, 'r') as f:
        json_data = json.load(f)
    aln_file = json_data['alignment_clustered_ldos_file']
    level = json_data['oglevel']
    odb_gene_id = json_data['query_odb_gene_id']
    return odb_gene_id, aln_file, level


def create_filemap(main_output_folder: str|Path, output_file: str|Path):
    json_dir = Path(main_output_folder) / 'info_jsons'
    json_files = get_json_file_list(json_dir)
    file_map = {}
    for json_file in json_files:
        odb_gene_id, aln_file, level = get_map_info_from_json(json_file)
        if odb_gene_id not in file_map:
            file_map[odb_gene_id] = {}
        if level in file_map[odb_gene_id]:
            raise ValueError(f'level {level} already in file map for gene {odb_gene_id}')
        file_map[odb_gene_id][level] = {}
        file_map[odb_gene_id][level]['alignment_file'] = aln_file
        file_map[odb_gene_id][level]['info_file'] = str(json_file.resolve())
    with open(output_file, 'w') as f:
        json.dump(file_map, f, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''
run this after the pipeline to create a single json file that maps the gene ids to the files output by the pipeline
creates a json file mapping odb gene ids to the files output by the pipeline
will create a file looking like:
{
    "odb_gene_id": {
        "Eukaryota": {
            "alignment_file": "path/to/file",
            "info_file": "path/to/file",
        },
        "Vertebrata": {
            "alignment_file": "path/to/file",
            "info_file": "path/to/file",
        },
        ...
    },
    ...
}''',
        formatter_class = argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '--main_output_folder',
        required=True,
        help='The main pipeline output folder. The folder should contain the folder `info_jsons/`, where it will look for the json files.'
    )
    parser.add_argument(
        '--output_file',
        required=True,
        help='path of the output json file'
    )
    args = parser.parse_args()
    create_filemap(args.main_output_folder, args.output_file)



