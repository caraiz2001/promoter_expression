# This script is used to extract kinase domains from protein sequences in a FASTA file.
# The script fetches domain information from UniProt for each protein sequence and extracts the kinase domain.
# The extracted kinase domains are saved to a new FASTA file.

# Input: A FASTA file containing protein sequences. The header of each sequence should contain the UniProt ID.
# Output: A new FASTA file containing the extracted kinase domains. (Both human and mouse have 315 sequences)

#%%
# ---------------------------------------
# Import libraries
# ---------------------------------------
import requests
from Bio import SeqIO
from io import StringIO
import pandas as pd


# ---------------------------------------
# Define functions
# ---------------------------------------

def get_domain_coordinates(uniprot_id):
    """Fetches kinase domain coordinates from UniProt for a given UniProt ID."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
    response = requests.get(url)

    if response.ok:
        # print(response.text)
        domain_coordinates = []
        lines = response.text.split("\n")

        for i, line in enumerate(lines):
            # Check for FT DOMAIN and ensure the next line contains "Protein kinase"
            if line.startswith("FT   DOMAIN"):
                print(f"The important line is {line}")
                parts = line.split()

                if len(parts) > 2:
                    range_str = parts[2]
                    if ".." in range_str:
                        start, end = map(int, range_str.split(".."))

                    # print(start, end)

                    # Check the next line for the "Protein kinase" keyword
                    if i + 1 < len(lines) and domain_name in lines[i + 1]:
                        print(f"The next line is {lines[i + 1]}")
                        domain_coordinates.append((start, end))

                        # print(f"The domain coordinates are {domain_coordinates}")
                    else:
                        print(f"Skipping domain {start}-{end} for {uniprot_id}")

        return domain_coordinates
    else:
        print(f"Error fetching data for {uniprot_id}: {response.status_code}")
        return []

def extract_domains_from_fasta(fasta_file, output_file):
    """Extracts kinase domains from sequences in a FASTA file."""
    num_searches = 0
    with open(output_file, "w") as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Extract UniProt ID from FASTA header (e.g., sp|O88866|HUNK_MOUSE)
            uniprot_id = record.id.split("|")[1] if "|" in record.id else record.id
            num_searches += 1
            # Get kinase domain coordinates
            domains = get_domain_coordinates(uniprot_id)
            print(f"the uniprot id is {uniprot_id} and the coordinates {domains}")
            # Extract and save each kinase domain
            for start, end in domains:
                kinase_seq = record.seq[start - 1:end]  # UniProt is 1-indexed
                out_fasta.write(f">{record.id}_domain_{start}_{end}\n{kinase_seq}\n")
    print(f"Processed {num_searches} sequences")


# ---------------------------------------
# Extract domains from protein sequences
# ---------------------------------------
family = "kinase"  # Family of interest
# Define the domain name to search for in UniProt entries
domain_name = "Protein kinase"

# input files
protein_fasta_uniprot_human = f"../../data/sequence/human/raw_protein_{family}_sequences.fasta"
protein_fasta_uniprot_mouse = f"../../data/sequence/mouse/mouse_raw_protein_{family}_sequences.fasta"

# output files
domain_fasta_uniprot_human = f"../../data/sequence/human/{family}_domains.fasta"
domain_fasta_uniprot_mouse = f"../../data/sequence/mouse/mouse_{family}_domains.fasta"

domain_fasta_unique_human = f"../../data/sequence/human/{family}_domains_unique.fasta"
domain_fasta_unique_mouse = f"../../data/sequence/mouse/mouse_{family}_domains_unique.fasta"

domain_fasta_ensemble_human = f"../../data/sequence/human/{family}_domains_ensemblID.fasta"
domain_fasta_ensemble_mouse = f"../../data/sequence/mouse/mouse_{family}_domains_ensemblID.fasta"



extract_domains_from_fasta(protein_fasta_uniprot_human, domain_fasta_uniprot_human)
extract_domains_from_fasta(protein_fasta_uniprot_mouse, domain_fasta_uniprot_mouse)

# For entries with the same uniprot ID, keep only the first one
# This is because some proteins have multiple kinase domains, we will keep the first one for simplicity

# Check that the number of sequences in the output file is correct
human = list(SeqIO.parse(domain_fasta_uniprot_human, "fasta"))
human_raw = list(SeqIO.parse(protein_fasta_uniprot_human, "fasta"))

unique_human = []
seen = set()
for seq in human:
    uniprot_id = seq.id.split("|")[1]
    if uniprot_id not in seen:
        unique_human.append(seq)
        seen.add(uniprot_id)

# Save the unique sequences to a new file
SeqIO.write(unique_human, domain_fasta_unique_human, "fasta")

# Repeat the same process for mouse
mouse = list(SeqIO.parse(domain_fasta_uniprot_mouse, "fasta"))
mouse_raw = list(SeqIO.parse(protein_fasta_uniprot_mouse, "fasta"))

# For entries with the same uniprot ID, keep only the first one
unique_mouse = []
seen = set()
for seq in mouse:
    uniprot_id = seq.id.split("|")[1]
    if uniprot_id not in seen:
        unique_mouse.append(seq)
        seen.add(uniprot_id)

# Save the unique sequences to a new file
SeqIO.write(unique_mouse, domain_fasta_unique_mouse, "fasta")

# Print number of sequences in the unique files
print(f"Number of unique human sequences: {len(unique_human)}")
print(f"Number of unique mouse sequences: {len(unique_mouse)}")

# Load orthologs_mapped file
orthologs_mapped = pd.read_csv("../../data/db_ids/human_mouse_one2one_filtered_mapped.txt", sep="\t")

# Write another fasta file in which the headers are ensemble gene ids instead of uniprot ids
with open(domain_fasta_ensemble_human, "w") as f:
    for seq in unique_human:
        uniprot_id = seq.id.split("|")[1]
        gene_id = orthologs_mapped[orthologs_mapped["Entry_human"] == uniprot_id]["Gene stable ID"].values[0]
        f.write(f">{gene_id}\n{seq.seq}\n")

# Repeat the same process for mouse
with open(domain_fasta_ensemble_mouse, "w") as f:
    for seq in unique_mouse:
        uniprot_id = seq.id.split("|")[1]
        gene_id = orthologs_mapped[orthologs_mapped["Entry_mouse"] == uniprot_id]["Mouse gene stable ID"].values[0]
        f.write(f">{gene_id}\n{seq.seq}\n")

# --------------------------
# Final notes
# --------------------------
# The script extracts domains from protein sequences in a FASTA file.
# The extracted domains are saved to a new FASTA file.
# The uniprot identifiers are taken from the headers of the raw fasta file (315 seq in mouse and human)
# The resulting domain files have more sequences, because some proteins have multiple kinase domains.
# We have kept only the first one. 
# %%
