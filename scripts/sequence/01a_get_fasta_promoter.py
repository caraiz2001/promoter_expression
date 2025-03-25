# This script is used to retrieve promoter sequences from Ensembl databases. 
# Input: 
# - orthologs_file
# - mapping_file

# Output: Fasta file with promoter sequences (ENSEMBL ID as identifier)
# - human_promoter_file
# - mouse_promoter_file


# ---------------------------------------
# Import libraries
# ---------------------------------------
#%%
import pandas as pd
import numpy as np
import requests
import re

# ---------------------------------------
# Load the orthologs and mapping information
# ---------------------------------------

#%%
orthologs_file = "../../data/db_ids/human_mouse_one2one_filtered.txt"
orthologs = pd.read_csv(orthologs_file, sep="\t")

mapping_human = "../../data/db_ids/human_idmapping_2025_01_24.tsv"
mapping_mouse = "../../data/db_ids/mouse_idmapping_2025_02_11.tsv"

mapping_human = pd.read_csv(mapping_human, sep="\t")
mapping_mouse = pd.read_csv(mapping_mouse, sep="\t")

# Combine the orthologs file with the mapping information
mapping_human.columns = [col + "_human" for col in mapping_human.columns]
mapping_mouse.columns = [col + "_mouse" for col in mapping_mouse.columns]
orthologs_mapped = orthologs.merge(mapping_human, left_on="Gene stable ID", right_on="From_human")
orthologs_mapped = orthologs_mapped.merge(mapping_mouse, left_on="Mouse gene stable ID", right_on="From_mouse")

# Save the mapped orthologs file
orthologs_mapped.to_csv("../../data/db_ids/human_mouse_one2one_filtered_mapped.txt", sep="\t", index=False)


# %%
# ---------------------------------------
# Define functions
# ---------------------------------------

# Function to get the sequence from the ensembl id
# This returns the 500 upstream promoter sequence + the gene sequence
def get_promoter_sequence(ensembl_id, upstream=500):
    url = f"https://rest.ensembl.org/sequence/id/{ensembl_id}?"
    params = {
        "type": "genomic",
        "expand_5prime": upstream,  # Fetch upstream promoter region
    }
    
    headers = {"Content-Type": "text/plain"}
    response = requests.get(url, params=params, headers=headers)
    
    if response.ok:
        return response.text.strip()
    else:
        print(f"Error fetching {ensembl_id}: {response.status_code}")
        return None

# %%
# ---------------------------------------
# Retrieve the promoter sequences
# ---------------------------------------
family = "kinase"

# Extract unique gene ids from the orthologs file
gene_ids = orthologs_mapped['Gene stable ID'].unique() # Human gene ids (315)
gene_ids_mouse = orthologs_mapped['Mouse gene stable ID'].unique() # Mouse gene ids (315)

# Print the number of genes
print(f"Number of human genes: {len(gene_ids)}")
print(f"Number of mouse genes: {len(gene_ids_mouse)}")

# Retrieve the promoter sequences for each gene
# Fetch promoter sequences
promoters = {gene: get_promoter_sequence(gene, upstream=500) for gene in gene_ids}
promoters_mouse = {gene: get_promoter_sequence(gene, upstream=500) for gene in gene_ids_mouse}

# Filter the promoter sequence to only include the first 500 nucleotides
promoters = {gene: seq[:500] for gene, seq in promoters.items()}
promoters_mouse = {gene: seq[:500] for gene, seq in promoters_mouse.items()}

human_promoter_file = f"../../data/sequence/human/only_promoter_{family}_human_sequences.fasta"
mouse_promoter_file = f"../../data/sequence/mouse/only_promoter_{family}_mouse_sequences.fasta"

# Write output to file in the data/sequence folder
with open(human_promoter_file, "w") as f:
    for gene, seq in promoters.items():
        f.write(f">{gene}\n{seq}\n")

# Write output to file
with open(mouse_promoter_file, "w") as f:
    for gene, seq in promoters_mouse.items():
        f.write(f">{gene}\n{seq}\n")
