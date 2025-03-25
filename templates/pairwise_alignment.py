#!/usr/bin/env python3

import pandas as pd
import numpy as np
from Bio.Align import PairwiseAligner
from Bio import pairwise2
import re
from Bio import SeqIO
from Bio.Align import substitution_matrices

# ---------------------------------------
# Define functions
# ---------------------------------------


def compute_identity(alignment):
    matches = sum(a == b for a, b in zip(alignment[0], alignment[1]) if a != '-' and b != '-')
    alignment_length = alignment.length  # Total aligned length
    return (matches / alignment_length) * 100 if alignment_length > 0 else 0

def compute_alignment_all_promoters(promoters_name, promoters_name_mouse):
    # """Computes pairwise alignments between human and mouse promoter sequences."""
    # Initialize the aligner
    aligner = PairwiseAligner()
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"

    results = []
    not_in_fasta = 0

    for gene_human in promoters_name.keys():
        for gene_mouse in promoters_name_mouse.keys():
            print(f"aligning {gene_human} and {gene_mouse}")
            alignments = aligner.align(promoters_name[gene_human], promoters_name_mouse[gene_mouse])
            alignment = alignments[0]
            score = alignment.score
            identity = compute_identity(alignment)
            results.append((gene_human, gene_mouse, score, identity))
        else:
            not_in_fasta += 1
        

    df = pd.DataFrame(results, columns=["Gene_human", "Gene_mouse", "alignment_score", "promoter_identity"])
    print(f"{not_in_fasta} gene pairs not found in the FASTA file")
    print(f"{len(results)} alignments computed")
    return df

def compute_alignment_all_domains(proteins, proteins_mouse):
    # """Computes pairwise alignments between human and mouse protein sequences."""

    # Initialize the aligner
    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10  # Opening a gap penalty
    aligner.extend_gap_score = -1  # Extending a gap penalty
    aligner.mode = "global"

    results = []
    not_in_fasta = 0

    for prot in proteins.keys():
        for prot_mouse in proteins_mouse.keys():
            print(f"aligning {prot} and {prot_mouse}")
            alignments = aligner.align(proteins[prot], proteins_mouse[prot_mouse])
            alignment = alignments[0]
            score = alignment.score
            identity = compute_identity(alignment)
            results.append((prot, prot_mouse, score, identity))
        else:
            not_in_fasta += 1
        

    df = pd.DataFrame(results, columns=["Gene_human", "Gene_mouse", "alignment_score", "domain_identity"])
    print(f"{not_in_fasta} gene pairs not found in the FASTA file")
    print(f"{len(results)} alignments computed")
    return df

def clean_sequence(seq):
    # """Removes non-standard amino acids"""
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq.upper())


def compute_alignment_all(sequence, sequence_mouse, sequence_type):
    if sequence_type == "protein":
        sequence = {gene: clean_sequence(seq) for gene, seq in sequence.items()}
        sequence_mouse = {gene: clean_sequence(seq) for gene, seq in sequence_mouse.items()}
        return compute_alignment_all_domains(sequence, sequence_mouse)
    else:
        return compute_alignment_all_promoters(sequence, sequence_mouse)
    

# ---------------------------------------
# Perform pairwise alignments
# ---------------------------------------

sequences = {}
for record in SeqIO.parse("${human_fasta}", "fasta"):
    gene_id = record.id
    print(gene_id)
    sequences[gene_id] = str(record.seq)



sequences_mouse = {}
for record in SeqIO.parse("${mouse_fasta}", "fasta"):
    gene_id = record.id
    print(gene_id)
    sequences_mouse[gene_id] = str(record.seq)

# Check if both dictionaries have the same keys
print(f"Sequence keys match: {sequences.keys() == sequences_mouse.keys()}")

# Perform pairwise alignments
alignments = compute_alignment_all(sequences, sequences_mouse, "${sequence_type}")

family_id = "${params.family_id}"
sequence_type = "${sequence_type}"

alignments.to_csv(f"{family_id}_{sequence_type}_alignment.tsv", sep="\\t", index=False)