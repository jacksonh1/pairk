
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, AlignIO, Seq, SeqIO

import local_conservation_analysis_pipeline.group_conservation_objects as group_tools
from local_conservation_scores.tools import capra_singh_2007_scores, general
from local_env_variables import matrices
from local_seqtools import alignment_tools as aln_tools
from local_seqtools import general_utils as tools
from local_seqtools import jch_alignment as jch_aln


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

def get_kmer_psuedo_alignment(kmer: str, idr_seqrecords: list, matrix: pd.DataFrame):
    scores = []
    subseqs = []
    positions = []
    for i in idr_seqrecords:
        score, subseq, pos = aln_kmer_2_seq_no_gaps(kmer, str(i.seq), substitution_matrix_df=matrix)
        scores.append(score)
        subseqs.append(subseq)
        positions.append(pos)
    return scores, subseqs, positions


json_file = "../../examples/table_annotation/conservation_analysis/3-9606_0_002f40/3-9606_0_002f40.json"
og = group_tools.ConserGene(json_file)
lvlo=og.get_level_obj('Vertebrata')
lvlo.alignment_file

# inputs
input_alignment_file=lvlo.alignment_file
reference_id=og.query_gene_id
idr_st=og.idr_start
idr_end=og.idr_end
mot_st=og.hit_start_position
mot_end=og.hit_end_position
scoring_matrix_name="grantham_similarity_norm"
l_flank= 4
r_flank= 4


# from idr alignment
fasta_importer = tools.FastaImporter(input_alignment_file)
idraln 



# main function
refseq_aln = fasta_importer.import_as_dict()[reference_id]
refseq, alnindex = tools.reindex_alignment_str(str(refseq_aln.seq))
mot = refseq[mot_st : mot_end + 1]
idr_aln_st = alnindex[idr_st]
idr_aln_end = alnindex[idr_end]
refseq_aln_idr = refseq_aln[idr_aln_st : idr_aln_end+1]
aln = fasta_importer.import_as_alignment()
idr_aln = aln[:, idr_aln_st : idr_aln_end+1]
idr_aln.sort(
    key=lambda record: aln_tools.percent_identity(
        record, refseq_aln_idr
    ),
    reverse=True,
)
# sort alignment by percent identity to query
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


kmers=get_kmers(str(ref_idr.seq), hit_st_in_idr, hit_end_in_idr)

matrix = matrices.load_precomputed_matrix_df(scoring_matrix_name)

d = {}
for kmer in kmers:
    scores, seqs, positions = get_kmer_psuedo_alignment(
        kmer=kmer,
        idr_seqrecords=list(idrs.values()),
        matrix=matrix
    )
    d[kmer] = positions
scores, seqs, positions = get_kmer_psuedo_alignment(
    kmer=hit_seq,
    idr_seqrecords=list(idrs.values()),
    matrix=matrix
)
d[hit_seq] = positions

df = pd.DataFrame(d)

df.to_csv('temp.csv', index=False)


# for score, seq, pos in zip(scores, seqs, positions):
    # print(f'{score:.2f}', seq, pos)
