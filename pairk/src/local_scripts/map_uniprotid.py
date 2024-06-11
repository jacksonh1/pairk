import pandas as pd
from local_orthoDB_group_pipeline import uniprotid_search
import argparse
from pathlib import Path


def map_uniprot_id(table: pd.DataFrame, uniprotid_column_name: str = 'uniprot_id') -> pd.DataFrame:
    uniprot_ids = list(table[uniprotid_column_name].unique())
    id_map = {}
    for uniprot_id in uniprot_ids:
        try:
            id_map[uniprot_id] = uniprotid_search.uniprotid_2_odb_gene_id(uniprot_id)
        except ValueError:
            continue
    table['gene_id'] = table[uniprotid_column_name].map(id_map)
    return table


def main(input_file: str, uniprotid_column_name: str = 'uniprot_id', output_file: str|None = None):
    input_file = Path(input_file)
    if output_file is None:
        output_file = input_file.parent / f'{input_file.stem}_mapped_odbgeneid{input_file.suffix}'
    else:
        output_file = Path(output_file)
    table = pd.read_csv(input_file)
    assert uniprotid_column_name in table.columns, f'column {uniprotid_column_name} not found in input table'
    table = map_uniprot_id(table, uniprotid_column_name)
    table.to_csv(output_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='''maps a column of uniprot ids in a table to orthoDB gene ids.
exports a copy of the table with a new column containing the orthoDB gene ids.''',
        formatter_class = argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        required=True,
        metavar="<file>",
        help='''input table with a column of uniprot ids. In csv format.'''
    )
    parser.add_argument(
        '--uni_column',
        type=str,
        default='uniprot_id',
        metavar="<str>",
        help='''name of the column containing the uniprot ids. Default: uniprot_id'''
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        default=None,
        metavar="<file>",
        help='''output file name. Default: input file name + "_mapped_odbgeneid"'''
    )
    args = parser.parse_args()
    main(args.input, args.uni_column, args.output)


