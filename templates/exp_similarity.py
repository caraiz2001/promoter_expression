#!/usr/bin/env python3

# ---------------------------------------
# Import packages + define parameters + load data
# ---------------------------------------

import pandas as pd
from sklearn.metrics.pairwise import cosine_similarity
import itertools
from scipy.stats import pearsonr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

human_tpm_file = "${human_transformed}"
mouse_tpm_file = "${mouse_transformed}"
orthologs_mapped_file = "${orthologs_mapped}"
family_id = "${family_id}"
transform_method = "${transform_method}"

human_tpm = pd.read_csv(human_tpm_file, sep="\\t")
mouse_tpm = pd.read_csv(mouse_tpm_file, sep="\\t")
orthologs_mapped = pd.read_csv(orthologs_mapped_file, sep="\\t")

# Filter the TPM data to keep only the genes ids present in the orthologs file
#%%
human_tpm_orthologs = human_tpm.loc[human_tpm["gene_id"].isin(orthologs_mapped["Gene stable ID"])]
mouse_tpm_orthologs = mouse_tpm.loc[mouse_tpm["gene_id"].isin(orthologs_mapped["Mouse gene stable ID"])]

human_tpm_orthologs = human_tpm_orthologs.copy()
mouse_tpm_orthologs = mouse_tpm_orthologs.copy()

# Add a column at the beginning of the df with the gene names instead of gene ids
human_tpm_orthologs["gene_name"] = human_tpm_orthologs["gene_id"].map(orthologs_mapped.set_index("Gene stable ID")["Gene name"])
mouse_tpm_orthologs["gene_name"] = mouse_tpm_orthologs["gene_id"].map(orthologs_mapped.set_index("Mouse gene stable ID")["Gene name"])

# Set the order of the columns to gene_id, gene_name and then the rest of the columns
human_tpm_orthologs = human_tpm_orthologs[["gene_id", "gene_name"] + [col for col in human_tpm_orthologs.columns if col not in ["gene_id", "gene_name"]]]
mouse_tpm_orthologs = mouse_tpm_orthologs[["gene_id", "gene_name"] + [col for col in mouse_tpm_orthologs.columns if col not in ["gene_id", "gene_name"]]]

# Print the number of genes in each data frame
print(f"Number of genes in human orthologs data frame: {len(human_tpm_orthologs)}")
print(f"Number of genes in mouse orthologs data frame: {len(mouse_tpm_orthologs)}")

human_filtered = human_tpm_orthologs.copy()
mouse_filtered = mouse_tpm_orthologs.copy()

# Compute the mean expression for each gene across all samples
human_mean_expression = human_filtered.iloc[:, 2:].mean(axis=1)
mouse_mean_expression = mouse_filtered.iloc[:, 2:].mean(axis=1)

# Print summary of the mean expression
print(f"Mean expression in human data frame: {human_mean_expression.describe()}")
print(f"Mean expression in mouse data frame: {mouse_mean_expression.describe()}")


# Keep only the columns that have a mean expression above the threshold. Ignore the columns gene_id and gene_name
numeric_cols = human_tpm_orthologs.iloc[3].apply(type) != str # True if the value is a string



# Which genes are in both data frames
common_genes = set(human_filtered["gene_name"]).intersection(set(mouse_filtered["gene_name"]))
print(f"Number of common genes: {len(common_genes)}")

human_filtered = human_filtered.loc[human_filtered["gene_name"].isin(common_genes)]
mouse_filtered = mouse_filtered.loc[mouse_filtered["gene_name"].isin(common_genes)]

# Order rows ion mouse according to human df
mouse_filtered = mouse_filtered.set_index("gene_name").loc[human_filtered["gene_name"]].reset_index()
# Rename column gene name to gene_name_human and gene_name_mouse
human_filtered = human_filtered.rename(columns={"gene_name": "gene_name_human"})
mouse_filtered = mouse_filtered.rename(columns={"gene_name": "gene_name_mouse"})

# Compute the similarity of expression. There are 2 options
# 1. Cosine similarity 
# 2. Pearson correlation



cosine_sim_matrix = cosine_similarity(human_filtered.loc[:, human_filtered.columns[numeric_cols]], mouse_filtered.loc[:, mouse_filtered.columns[numeric_cols]])
cosine_sim_df = pd.DataFrame(cosine_sim_matrix, index=human_filtered["gene_name_human"], columns=mouse_filtered["gene_name_mouse"])

cosine_sim_long = cosine_sim_df.reset_index().melt(id_vars='gene_name_human', var_name='gene_name_mouse', value_name='exp_sim_cosine')

human_num = human_filtered.loc[:, human_filtered.columns[numeric_cols]]
mouse_num = mouse_filtered.loc[:, mouse_filtered.columns[numeric_cols]]

# Set the index to be the gene names
human_num.index = human_filtered["gene_name_human"]
mouse_num.index = mouse_filtered["gene_name_mouse"]


correlation_matrix = np.corrcoef(human_num, mouse_num, rowvar=True)

# Extract the relevant portion (human vs. mouse correlations)
human_count = human_num.shape[0]
mouse_count = mouse_num.shape[0]

# The upper-left block of the correlation matrix is human-human, bottom-right is mouse-mouse
# The relevant block (human vs. mouse) is in the top-right (or bottom-left since it's symmetric)
human_mouse_correlations = correlation_matrix[:human_count, human_count:]

# Convert matrix to DataFrame
correlation_df = pd.DataFrame(human_mouse_correlations, 
                              index=human_num.index, 
                              columns=mouse_num.index)
# Display results
print(correlation_df.head())

# Melt the DataFrame to long format
correlation_df_long = correlation_df.reset_index().melt(id_vars="gene_name_human", 
                                                   var_name="gene_name_mouse", 
                                                   value_name="exp_sim_pearson")

# Merge the two similarity DataFrames
exp_sim = cosine_sim_long.merge(correlation_df_long, on=["gene_name_human", "gene_name_mouse"])

# Filter for the rows with gene_name_human == gene_name_mouse
sim_orthologs = exp_sim[exp_sim["gene_name_human"] == exp_sim["gene_name_mouse"]]
sim_non_orthologs = exp_sim[exp_sim["gene_name_human"] != exp_sim["gene_name_mouse"]]

# Plot density plot

plt.figure()
sns.kdeplot(sim_non_orthologs['exp_sim_cosine'], color = "coral", label = "non-orthologs")
sns.kdeplot(sim_orthologs['exp_sim_cosine'], color = "blue", label = "orthologs")
plt.legend(['Orthologs = False', 'Orthologs = True'])
plt.title('Expression similarity (Cosine)')
plt.savefig(f'{family_id}_{transform_method}_exp_sim_cosine_density.png', dpi=300, bbox_inches='tight')

plt.figure()
sns.kdeplot(sim_non_orthologs['exp_sim_pearson'], color = "coral", label = "non-orthologs")
sns.kdeplot(sim_orthologs['exp_sim_pearson'], color = "blue", label = "orthologs")
plt.legend(['Orthologs = False', 'Orthologs = True'])
plt.title('Expression similarity (Pearson)')
plt.savefig(f'{family_id}_{transform_method}_exp_sim_cosine_Pearson.png', dpi=300, bbox_inches='tight')



# Save the similarity matrix to a file
exp_sim.to_csv(f"{family_id}_{transform_method}_exp_similarity.tsv", index=False, sep="\\t")


# ---------------------------
# Write report to a file
# ---------------------------

# - Number of genes in the input df (transformed and filtered for removing non expressed)
# - Number of orthologs
# - Avereage expression gene-wise -> describe
# - Exp similarity -> describe
# - Plot orthologs vs non-orthologs (cos/pearson)
#


text = f"""

    The number of genes in the input df (transformed and filtered for non-expressed genes) is:
    {human_tpm.shape[0]} (human) and {mouse_tpm.shape[0]} (mouse)

    The number of orthologs:
    {human_tpm_orthologs.shape[0]} (human) and {mouse_tpm_orthologs.shape[0]} (mouse)

    Mean expression in human data frame:
    {human_mean_expression.describe()} (human)
    {mouse_mean_expression.describe()} (mouse)

    Expression similarity across genes pairs:
    - Cosine: {exp_sim['exp_sim_cosine'].describe()}
    - Pearson: {exp_sim['exp_sim_cosine'].describe()}

"""

# Write the report to a file
report_file = f"{family_id}_{transform_method}_report.txt"
with open(report_file, "w") as file:
    file.write(text)

print(f"Report written to {report_file}")