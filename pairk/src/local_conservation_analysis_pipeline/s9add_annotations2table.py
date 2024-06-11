import json
from pathlib import Path

import numpy as np
import pandas as pd
from openpyxl import Workbook


def main(
    annotations_file,
    table_file,
    table_annotations: list[str],
    table_annotation_levels: list[str],
    output_table_file,
):
    with open(annotations_file, "r") as f:
        annotations = json.load(f)
    table_df = pd.read_csv(table_file)
    error_map = {}
    df = pd.DataFrame()
    for k, d in annotations.items():
        ref_ind = int(k)
        if d['critical_error'] is not None:
            error_map[ref_ind] = d['critical_error']
            continue
        if 'level_annotations' in d:
            lvl_annotations = d.pop('level_annotations')
            for level, lvl_d in lvl_annotations.items():
                if level not in table_annotation_levels:
                    continue
                for annotation in table_annotations:
                    if annotation not in lvl_d:
                        continue
                    if isinstance(lvl_d[annotation], list):
                        lvl_d[annotation] = str(lvl_d[annotation])
                    df.loc[ref_ind, f"{level}_{annotation}"] = lvl_d[annotation]
        for annotation in table_annotations:
            if annotation not in d:
                continue
            df.loc[ref_ind, annotation] = d[annotation]
    table_df['critical_error'] = table_df['reference_index'].map(error_map)
    df=df.reset_index(names='reference_index')
    table_df = pd.merge(table_df, df, on='reference_index', how='left')
    table_df.to_csv(output_table_file, index=False)
    # table_df.to_excel(output_table_file.replace(".csv", ".xlsx"), index=False)
    # workbook = Workbook()
    # sheet = workbook.active
    # for index, row in table_df.iterrows():
    #     sheet.append(row.to_dict())
    # workbook.save(output_table_file.replace(".csv", ".xlsx"))
    return table_df
















