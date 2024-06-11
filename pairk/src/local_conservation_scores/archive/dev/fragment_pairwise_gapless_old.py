# %%
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO, SeqRecord

from local_conservation_scores.tools import capra_singh_2007_scores, general
from local_env_variables import matrices
from local_seqtools import alignment_tools as aln_tools
from local_seqtools import general_utils as tools
from local_seqtools import jch_alignment as jch_aln


# def score_




def aln_kmer_2_seq_no_gaps(
    kmer: str,
    sequence: str,
    substitution_matrix_df: pd.DataFrame,
):
    positions = []
    subsequences = []
    scores = []
    if len(sequence) < len(kmer):
        return 0, "-"*len(kmer), -1
    for i in range(len(sequence) - len(kmer) + 1):
        subseq = sequence[i : i + len(kmer)]
        max_score = max(
            aln_tools.score_alignment(subseq, subseq, substitution_matrix_df),
            aln_tools.score_alignment(kmer, kmer, substitution_matrix_df),
        )
        score = aln_tools.score_alignment(kmer, subseq, substitution_matrix_df)/max_score
        positions.append(i)
        scores.append(score)
        subsequences.append(subseq)
    return np.max(scores), subsequences[np.argmax(scores)], positions[np.argmax(scores)]


def pad_hit(seq: str, st_pos: int, end_pos: int, l_flank: int = 0, r_flank: int = 0):
    st = max(0, st_pos - l_flank)
    end = min(len(seq)-1, end_pos + r_flank)
    return st, end, seq[st : end + 1]


def split_query_seq_for_kmer_gen(query_seq: str, hit_start: int, hit_end: int):
    pt1 = query_seq[:hit_start]
    pt2 = query_seq[hit_end + 1 :]
    return pt1, pt2


def get_kmers(query_seq, hit_st, hit_end):
    k = hit_end - hit_st + 1
    pt1, pt2 = split_query_seq_for_kmer_gen(query_seq, hit_st, hit_end)
    kmers = []
    kmers.extend(tools.gen_kmers(pt1, k))
    kmers.extend(tools.gen_kmers(pt2, k))
    kmers.append(query_seq[hit_st : hit_end + 1])
    return list(set(kmers))


def get_kmer_psuedo_alignment(kmer: str, idrs: list[SeqRecord]):


def main2(
    idr_aln: Align.MultipleSeqAlignment,
    k: int,
):





def main(
    input_alignment_file: str|Path,
    reference_id: str,
    idr_st: int,
    idr_end: int,
    mot_st: int,
    mot_end: int,
    scoring_matrix_name: str="EDSSMat50_max_off_diagonal_norm",
    l_flank: int = 0,
    r_flank: int = 0,
):
    fasta_importer = tools.FastaImporter(input_alignment_file)
    ref_seq_aln = fasta_importer.import_as_dict()[reference_id]
    ref_seq, alnindex = tools.reindex_alignment_str(str(ref_seq_aln.seq))
    mot = ref_seq[mot_st : mot_end + 1]
    idr_aln_st = alnindex[idr_st]
    idr_aln_end = alnindex[idr_end]
    aln = fasta_importer.import_as_alignment()
    idr_aln = aln[:, idr_aln_st : idr_aln_end+1]
    idrs = tools.strip_dashes_from_sequences(list(idr_aln))
    idrs = {i.id: i for i in idrs}
    ref_idr = idrs.pop(reference_id)
    
    mot_st_in_idr = mot_st - idr_st
    mot_end_in_idr = mot_end - idr_st
    hit_st_in_idr, hit_end_in_idr, hit_seq = pad_hit(
        str(ref_idr.seq), mot_st_in_idr, mot_end_in_idr, l_flank, r_flank
    )
    mot_st_in_hit = mot_st_in_idr - hit_st_in_idr
    mot_end_in_hit = mot_end_in_idr - hit_st_in_idr












    aln = fasta_importer.import_as_alignment()
    jchaln = jch_aln.jch_alignment(aln, reference_id)
    ref_seq = jchaln.query_unaligned_str
    
    stripped_fl_seqs = tools.strip_dashes_from_sequences(fasta_importer.import_as_list())

    hit_aln_st_idr = hit_aln_st - idr_aln_st
    hit_aln_end_idr = hit_aln_end - idr_aln_st
    hit_aln_seq = str(ref_seq_aln.seq)[hit_aln_st_idr : hit_aln_end_idr]
    ref_hit_seq = str(ref_idr.seq)[hit_aln_st_idr : hit_aln_end_idr]
    scoring_matrix_df = matrices.load_precomputed_matrix_df(scoring_matrix_name)
    print(ref_hit_seq)


# %%
# ==============================================================================
# // testing
# ==============================================================================    
import local_conservation_analysis_pipeline.group_conservation_objects as group_tools

json_file = "../../examples/table_annotation/conservation_analysis/3-9606_0_002f40/3-9606_0_002f40.json"
og = group_tools.ConserGene(json_file)
lvlo=og.get_level_obj('Vertebrata')
lvlo.alignment_file
main(
    lvlo.alignment_file,
    lvlo.idr_aln_start,
    lvlo.idr_aln_end,
    lvlo.hit_aln_start,
    lvlo.hit_aln_end,
    og.query_gene_id,
    "EDSSMat50_max_off_diagonal_norm"
)



og.hit_sequence





# %%


















    















