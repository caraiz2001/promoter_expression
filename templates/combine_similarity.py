#!/usr/bin/env python3

# This script is used to combine the sequence similarity with the expression similarity.
# The input files are:
# - promoter alignmnet
# - domain alignment
# - expression similarity

# The output file is:
# A single df with Gene name mouse, Gene name human, Gene ID mouse, Gene ID human, Domain similarity, Promoter similarity, Expression similarity (cosine) and expression similarity (pearson)

# Load the data
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
print(os.getcwd()) 

family_id = "${family_id}"
transform_method = "${transform_method}"

# Load the data
#%%
promoter_sim_file = "${promoter_similarity}"
domain_sim_file = "${protein_similarity}"
expression_sim_file = "${expression_similarity}"
orthologs_mapped_file = "${orthologs_mapped}"

promoter_similarity = pd.read_csv(promoter_sim_file, sep="\\t")
domain_similarity = pd.read_csv(domain_sim_file, sep="\\t")
expression_similarity = pd.read_csv(expression_sim_file, sep="\\t")
orthologs_mapped = pd.read_csv(orthologs_mapped_file, sep="\\t")

# %%
# Merge the data

# Rename the promoter similarity columns to gene_id_human and gene_id_mouse
promoter_similarity = promoter_similarity.rename(columns={"Gene_human": "gene_id_human", "Gene_mouse": "gene_id_mouse"})

# Remove the column alignment_score
promoter_similarity = promoter_similarity.drop(columns=["alignment_score"])

# Same for domain similarity
domain_similarity = domain_similarity.rename(columns={"Gene_human": "gene_id_human", "Gene_mouse": "gene_id_mouse"})
domain_similarity = domain_similarity.drop(columns=["alignment_score"])

# Merge the data frames
sequence_similarity = promoter_similarity.merge(domain_similarity, on=["gene_id_human", "gene_id_mouse"], how="inner")

# %%

# Map each gene name (columns gene_name_human and gene_name_mouse in the expression similarity df) to the gene id using orthologs file
# Do a dictionary with the gene name as key and the gene id as value
orthologs_mapped = orthologs_mapped.set_index("Gene name")
gene_id_human_dict = orthologs_mapped["Gene stable ID"].to_dict()
gene_id_mouse_dict = orthologs_mapped["Mouse gene stable ID"].to_dict()

# Map the gene names to the gene ids
expression_similarity["gene_id_human"] = expression_similarity["gene_name_human"].map(gene_id_human_dict)
expression_similarity["gene_id_mouse"] = expression_similarity["gene_name_mouse"].map(gene_id_mouse_dict)

# Add one extra column named "orthologs" to the expression similarity df
expression_similarity["orthologs"] = expression_similarity["gene_name_human"] == expression_similarity["gene_name_mouse"]

# Reorder the columns to have gene ids, then gene names, then orthologs and then the similarity values
expression_similarity = expression_similarity[["gene_id_human", "gene_id_mouse", "gene_name_human", "gene_name_mouse", "orthologs", "exp_sim_cosine", "exp_sim_pearson"]]

# Merge the expression similarity with the sequence similarity
combined_df = expression_similarity.merge(sequence_similarity, on=["gene_id_human", "gene_id_mouse"], how="inner")


# Save the combined df
combined_df.to_csv(f"{family_id}_{transform_method}_combined_similarity.tsv", sep="\\t", index=False)

#%%

# Separte the df into orthologs = True and orthologs = False
orthologs_true = combined_df[combined_df['orthologs'] == True]
orthologs_false = combined_df[combined_df['orthologs'] == False]


# Plot promoter similarity against expression similarity (scatter plot)
plt.figure()
plt.scatter(orthologs_false['promoter_identity'], orthologs_false['exp_sim_cosine'], s=2, alpha=0.5, c = 'coral')
plt.scatter(orthologs_true['promoter_identity'], orthologs_true['exp_sim_cosine'], s=3, alpha=1, c='blue')
plt.legend(['Orthologs = False', 'Orthologs = True'])
plt.xlabel('Promoter Similarity')
plt.ylabel('Expression Similarity (cosine)')
plt.title('Promoter Similarity vs Expression Similarity')

plt.savefig(f'{family_id}_{transform_method}_promoter_vs_exp_sim_cosine.png', dpi=300, bbox_inches='tight')

# Plot promoter similarity against expression similarity (scatter plot)
plt.figure()
plt.scatter(orthologs_false['promoter_identity'], orthologs_false['exp_sim_pearson'], s=2, alpha=0.5, c = 'coral')
plt.scatter(orthologs_true['promoter_identity'], orthologs_true['exp_sim_pearson'], s=3, alpha=1, c='blue')
plt.legend(['Orthologs = False', 'Orthologs = True'])
plt.xlabel('Promoter Similarity')
plt.ylabel('Expression Similarity (pearson)')
plt.title('Promoter Similarity vs Expression Similarity')

plt.savefig(f'{family_id}_{transform_method}_promoter_vs_exp_sim_pearson.png', dpi=300, bbox_inches='tight')