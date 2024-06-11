# %%

import json
from pathlib import Path

import alv
from Bio import Align, AlignIO, SeqIO

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
import local_env_variables.env_variables as env
import local_seqtools.alignment_tools as aln_tools
import local_seqtools.general_utils as tools

# import local_seqtools.jch_alignment as jch_alignment

json_file = "./conservation_analysis/2-9606_0_002f40/2-9606_0_002f40.json"
n_flanking_aas = 20


def save_colored_protein_msa_html(
    alignment: Align.MultipleSeqAlignment,
    output_html_file,
    color_scheme="clustal",
    max_id_length=None,
    highlight_region=None,
):
    """
    Save a colored view of a protein multiple sequence alignment (MSA) as HTML.

    Parameters:
    - alignment_file (str): Path to the input protein MSA file (e.g., in FASTA format).
    - output_html_file (str): Path to save the output HTML file with colored protein MSA.
    - color_scheme (str): Color scheme for the alignment. Options: "default", "clustal", "mismatch".

    Returns:
    None
    """
    if color_scheme == "clustal":
        with open(env.COLOR_MAP_FILES["clustal"], "r") as f:
            colors = json.load(f)
    elif color_scheme == "mismatch":
        colors = {
            "A": "blue",
            "R": "red",
            "N": "green",
            "D": "orange",
            "C": "purple",
            "Q": "cyan",
            "E": "yellow",
            "G": "brown",
            "H": "pink",
            "I": "gray",
            "L": "olive",
            "K": "darkgreen",
            "M": "darkblue",
            "F": "darkred",
            "P": "darkorange",
            "S": "darkmagenta",
            "T": "darkyellow",
            "W": "darkcyan",
            "Y": "darkgray",
            "V": "lightgray",
            "-": "lightgray",
        }
    else:
        # Default color scheme
        colors = {
            "A": "green",
            "R": "blue",
            "N": "purple",
            "D": "red",
            "C": "orange",
            "Q": "magenta",
            "E": "yellow",
            "G": "cyan",
            "H": "lightblue",
            "I": "brown",
            "L": "pink",
            "K": "gray",
            "M": "olive",
            "F": "darkgreen",
            "P": "darkblue",
            "S": "darkred",
            "T": "darkorange",
            "W": "darkmagenta",
            "Y": "darkyellow",
            "V": "darkcyan",
            "-": "lightgray",
        }

    html_content = "<html><head><style>pre {font-family: 'Courier New', monospace;}</style></head><body><pre>\n"
    max_id_length = max_id_length or max(len(record.id) for record in alignment)
    for record in alignment:
        formatted_id = record.id[:max_id_length].ljust(max_id_length)
        sequence_line = ""
        for c, symbol in enumerate(record.seq):
            style = f'color: {colors.get(symbol, "black")};'
            if highlight_region is not None:
                if highlight_region[0] <= c + 1 <= highlight_region[1]:
                    style += "background-color: #f9f90262;"
            # if highlight_region and c in highlight_region:
            #     style += "font-weight: bold;"
            sequence_line += f'<span style="{style}">{symbol}</span>'
        html_content += f"{formatted_id} {sequence_line}\n"
    html_content += "</pre></body></html>"
    with open(output_html_file, "w") as html_file:
        html_file.write(html_content)


# def index2alnindex(aln, query_id, hit_start, hit_end, n_flanking_aas):
#     flaln = jch_alignment.jch_alignment(aln, query_id)
#     query_seq = flaln.query_unaligned_str
#     slice_start = max(0, hit_start - n_flanking_aas)
#     slice_end = min(len(query_seq) - 1, hit_end + n_flanking_aas)
#     flaln.slice_by_unaligned_positions_inclusive(slice_start, slice_end)


def index2alnindex_V2(
    lvlo: group_tools.ConserLevel, hit_start, hit_end, n_flanking_aas
):
    unaligned_query_seq, aln_index = tools.reindex_alignment_str(
        lvlo.query_aln_sequence
    )
    slice_start = max(0, hit_start - n_flanking_aas)
    slice_end = min(len(unaligned_query_seq) - 1, hit_end + n_flanking_aas)
    aln_slice_start = aln_index[slice_start]
    aln_slice_end = aln_index[slice_end]
    return aln_slice_start, aln_slice_end


def main(json_file, n_flanking_aas, whole_idr=False):
    og = group_tools.ConserGene(json_file)
    og.load_levels()
    og.query_gene_id
    hit_start = og.hit_start_position
    hit_end = og.hit_end_position
    for level in og.levels_passing_filters:
        lvl = og.level_objects[level]
        aln = lvl.aln
        # sort aln by percent identity to query
        aln.sort(
            key=lambda record: aln_tools.percent_identity(
                record, lvl.query_aln_sequence
            ),
            reverse=True,
        )
        if whole_idr:
            aln_slice_start, aln_slice_end = index2alnindex_V2(
                lvl, og.idr_start, og.idr_end, n_flanking_aas
            )
        else:
            aln_slice_start, aln_slice_end = index2alnindex_V2(
                lvl, hit_start, hit_end, n_flanking_aas
            )
        slice_file = (
            Path(og.analysis_folder)
            / f"{og.reference_index}-{og.query_gene_id.replace(':','')}-{level}_aln_slice.html"
        )
        og.add_item_to_lvl_orthogroup("aln_slice_file", str(slice_file), level)
        soffset = lvl.hit_aln_start - aln_slice_start
        eoffset = lvl.hit_aln_end - aln_slice_start
        save_colored_protein_msa_html(
            aln[:, aln_slice_start : aln_slice_end + 1],
            slice_file,
            color_scheme="clustal",
            highlight_region=(soffset, eoffset + 1),
        )
