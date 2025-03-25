# This script is used to retrieve protein sequences from UniprotKB/Swissprot databases. 
# Input: orthologs + mapping files
# Output: 
# - Fasta file with protein sequences (Uniprot ID as identifier)
# - Fasta file with protein sequences (Ensembl Gene ID as identifier)

# %%
# -----------------------------
# Import libraries
# -----------------------------
import pandas as pd
import numpy as np
import requests
import re


# %%
#---------------------------------------
# Define functions
#---------------------------------------

# Function to get the sequence from the uniprot id 
def get_swissprot_sequence(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)

    if response.ok:
        return response.text
    else:
        return f"Error {response.status_code}: {response.text}"

#%%
# ---------------------------------------
# Define input, output files and parameters
# ---------------------------------------
family = "kinase"  # Family of interest

# Input file
orthologs_mapped_file = "../../data/db_ids/human_mouse_one2one_filtered_mapped.txt"
orthologs_mapped = pd.read_csv(orthologs_mapped_file, sep="\t")

# Output files
raw_fasta_file = f"../../data/sequence/human/raw_protein_{family}_sequences.fasta"
renamed_fasta_file = f"../../data/sequence/human/renamed_protein_{family}_sequences.fasta"

raw_fasta_mouse_file = f"../../data/sequence/mouse/mouse_raw_protein_{family}_sequences.fasta"
renamed_fasta_mouse_file = f"../../data/sequence/mouse/mouse_renamed_protein_{family}_sequences.fasta"

# %%
# ---------------------------------------
# Extract protein sequences
# ---------------------------------------

# Extract tuples of [gene_id, uniprot_id] from the orthologs file
gene_prot_human_tuples = orthologs_mapped[["Gene stable ID", "Entry_human"]].values
gene_prot_mouse_tuples = orthologs_mapped[["Mouse gene stable ID", "Entry_mouse"]].values


# Open human FASTA files for writing
with open(raw_fasta_file, "w") as raw_fasta, open(renamed_fasta_file, "w"
) as renamed_fasta:
    for gene_id, protein_id in gene_prot_human_tuples:
        fasta_entry = get_swissprot_sequence(protein_id)

        if fasta_entry:
            # Write raw FASTA entry to file
            raw_fasta.write(fasta_entry)

            # Modify header with Ensembl Gene ID
            fasta_lines = fasta_entry.split("\n")
            fasta_lines[0] = f">{gene_id}"  # Replace header with Ensembl Gene ID
            renamed_entry = "\n".join(fasta_lines)

            # Write modified FASTA entry to second file
            renamed_fasta.write(renamed_entry)

print("FASTA files saved:")
print(f"- {raw_fasta_file} (original UniProt headers)")
print(f"- {renamed_fasta_file} (headers replaced with Ensembl Gene IDs)")


# Open mouse FASTA files for writing
with open(raw_fasta_mouse_file, "w") as raw_fasta, open(renamed_fasta_mouse_file, "w"
) as renamed_fasta:
    for gene_id, protein_id in gene_prot_mouse_tuples:
        fasta_entry = get_swissprot_sequence(protein_id)  # Retrieve full FASTA entry

        if fasta_entry:
            # Write raw FASTA entry to file
            raw_fasta.write(fasta_entry)

            # Modify header with Ensembl Gene ID
            fasta_lines = fasta_entry.split("\n")
            fasta_lines[0] = f">{gene_id}"  # Replace header with Ensembl Gene ID
            renamed_entry = "\n".join(fasta_lines)

            # Write modified FASTA entry to second file
            renamed_fasta.write(renamed_entry)

print("FASTA files saved:")
print(f"- {raw_fasta_mouse_file} (original UniProt headers)")
print(f"- {renamed_fasta_mouse_file} (headers replaced with Ensembl Gene IDs)")
