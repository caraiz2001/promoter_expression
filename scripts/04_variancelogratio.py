#%%
# # Load expression data
import pandas as pd
import numpy as np
import os
from itertools import combinations
from scipy.stats import linregress

#%%

# Load the TPM data
human_tpm = pd.read_csv("data/expression/human/merged_human_data.tsv", sep="\t")
mouse_tpm = pd.read_csv("data/expression/mouse/merged_mouse_data.tsv", sep="\t")
# # If the first column is not numeric, set it as row names
# if not np.issubdtype(human_tpm.iloc[:, 0].dtype, np.number):
#     human_tpm.set_index(human_tpm.columns[0], inplace=True)
#     mouse_tpm.set_index(mouse_tpm.columns[0], inplace=True)
#     # human_tpm = human_tpm.iloc[:, 1:]
#     # mouse_tpm = mouse_tpm.iloc[:, 1:]

# Define non-expressed genes as those in which all samples have expression below a threshold
expression_threshold = 1

# Identify non-expressed genes
human_non_expressed = (human_tpm.iloc[:, 1:] <= expression_threshold).all(axis=1)
mouse_non_expressed = (mouse_tpm.iloc[:, 1:] <= expression_threshold).all(axis=1)

# Filter out non-expressed genes
human_tpm_filtered = human_tpm[~human_non_expressed]
mouse_tpm_filtered = mouse_tpm[~mouse_non_expressed]



# %%
# Load the orthologs file
orthologs_mapped_file = "data/db_ids/human_mouse_one2one_filtered_mapped.txt"
orthologs_mapped = pd.read_csv(orthologs_mapped_file, sep="\t")

# Filter the TPM data to keep only the genes ids present in the orthologs file
human_tpm_orthologs = human_tpm_filtered.loc[human_tpm['gene_id'].isin(orthologs_mapped["Gene stable ID"])]
mouse_tpm_orthologs = mouse_tpm_filtered.loc[mouse_tpm['gene_id'].isin(orthologs_mapped["Mouse gene stable ID"])]

# Map the gene ids to gene names
human_tpm_orthologs["gene_name_human"] = human_tpm_orthologs['gene_id'].map(orthologs_mapped.set_index("Gene stable ID")["Gene name"])
mouse_tpm_orthologs["gene_name_mouse"] = mouse_tpm_orthologs['gene_id'].map(orthologs_mapped.set_index("Mouse gene stable ID")["Gene name"])

# Set the order of the columns to gene_id, gene_name and then the rest of the columns
human_tpm_orthologs = human_tpm_orthologs[["gene_name_human"] + [col for col in human_tpm_orthologs.columns if col != "gene_name_human"]]
mouse_tpm_orthologs = mouse_tpm_orthologs[["gene_name_mouse"] + [col for col in mouse_tpm_orthologs.columns if col != "gene_name_mouse"]]
# Print the number of genes in each data frame
print(f"Number of genes in human orthologs data frame: {len(human_tpm_orthologs)}")
print(f"Number of genes in mouse orthologs data frame: {len(mouse_tpm_orthologs)}")


# Remove the gene_id column
human_tpm_orthologs = human_tpm_orthologs.drop(columns=["gene_id"])
mouse_tpm_orthologs = mouse_tpm_orthologs.drop(columns=["gene_id"])


# # Set the index to gene names
# human_tpm_orthologs.set_index("gene_name_human", inplace=True)
# mouse_tpm_orthologs.set_index("gene_name_mouse", inplace=True)

# Merge the 2 df on the gene names 
merged_tpm_orthologs = human_tpm_orthologs.merge(mouse_tpm_orthologs, left_on="gene_name_human", right_on="gene_name_mouse", how="inner", suffixes=("_human", "_mouse"))

# Remove the gene_name_mouse column
merged_tpm_orthologs = merged_tpm_orthologs.drop(columns=["gene_name_mouse"])
# Rename the gene_name_human column to gene_name
merged_tpm_orthologs = merged_tpm_orthologs.rename(columns={"gene_name_human": "gene_name"})

# Set the index to gene names
merged_tpm_orthologs.set_index("gene_name", inplace=True)

# # Intersection of the rownames of both data frames
# common_genes = set(human_tpm_orthologs.index).intersection(set(mouse_tpm_orthologs.index))

#%%
# Calculate log ratios for each pair of genes

# human_tpm_orthologs = merged_tpm_orthologs.iloc[:, :5]
# mouse_tpm_orthologs = merged_tpm_orthologs.iloc[:, 5:]

def calculate_log_ratios(df):
    gene_pairs = list(combinations(df.index, 2))
    log_ratios = []

    for gene1, gene2 in gene_pairs:
        row = {"Gene1": gene1, "Gene2": gene2}
        for sample in df.columns:
            row[f"ratio_{sample}"] = np.log((df.loc[gene1, sample] + 1) / (df.loc[gene2, sample] + 1))
            # error message if the log ratio is NaN
            if np.isnan(row[f"ratio_{sample}"]):
                print(f"NaN value encountered for genes {gene1} and {gene2} in sample {sample}.")
                break
        log_ratios.append(row)

    return pd.DataFrame(log_ratios)

merged_log_ratios = calculate_log_ratios(merged_tpm_orthologs)


# Plot the first row (first gene pair)
import matplotlib.pyplot as plt

# Extract the first row
first_human_row = merged_log_ratios.iloc[100, 2:7]  # Skip Gene1 and Gene2 columns
first_mouse_row = merged_log_ratios.iloc[100, 7:]  # Skip Gene1 and Gene2 columns

# Create the plot
plt.figure(figsize=(8, 6))
plt.scatter(first_human_row, first_mouse_row, color='blue', alpha=0.7)
plt.title(f"Log Ratios for Gene Pair: {merged_log_ratios.iloc[0]['Gene1']} vs {merged_log_ratios.iloc[0]['Gene2']}")
plt.xlabel("Human Log Ratios")
plt.ylabel("Mouse Log Ratios")
plt.grid(True)
plt.show()

correlation_data = []

for i in range(len(merged_log_ratios)):
    human_log_ratios = merged_log_ratios.iloc[i, 2:7]  # Skip Gene1 and Gene2 columns
    mouse_log_ratios = merged_log_ratios.iloc[i, 7:]  # Skip Gene1 and Gene2 columns
    gene1 = merged_log_ratios.iloc[i]["Gene1"]
    gene2 = merged_log_ratios.iloc[i]["Gene2"]

    human_values = pd.to_numeric(human_log_ratios)  # Skip Gene1 and Gene2 columns
    mouse_values = pd.to_numeric(mouse_log_ratios) # Skip Gene1 and Gene2 columns
    
    # Calculate Pearson correlation and slope
    slope, intercept, r_value, p_value, std_err = linregress(human_values, mouse_values)
    
    correlation_data.append({
        "Gene1": gene1,
        "Gene2": gene2,
        "Pearson_Correlation_R2": r_value**2,
        "Slope": slope
    })

# Create a new DataFrame with the results
correlation_df = pd.DataFrame(correlation_data)

# Load the combined similarity data
combined_df_file = "results/combined_similarity/kinase_clr_combined_similarity.tsv"
combined_df = pd.read_csv(combined_df_file, sep="\t")

model_results_df = pd.read_csv("./results/model_results/kinase_vlr/all_pairs_all_metrics_vlr_rho_rank_5000.tsv", sep = "\t")

# KEep only the columns of interest
model_results_small = model_results_df[["gene_name_human", "gene_name_mouse", "exp_sim_cosine", "promoter_identity", "emb_sim_cos", "Rho", "label"]]

# 
model_results_forward = model_results_small.iloc[::2, :]
model_results_reverse = model_results_small.iloc[1::2, :]

# Merge the df on the gene names with the suffixes _forward and _reverse


model_results_small = model_results_forward.merge(model_results_reverse, left_on=["gene_name_human", "gene_name_mouse"], right_on=["gene_name_mouse", "gene_name_human"], how="inner", suffixes=("_forward", "_reverse"))

# Calculate the avergae of emb_sim_cos_forward and emb_sim_cos_reverse
model_results_small["emb_sim_cos"] = (model_results_small["emb_sim_cos_forward"] + model_results_small["emb_sim_cos_reverse"]) / 2
# Calculate the average of exp_sim_cosine_forward and exp_sim_cosine_reverse
model_results_small["exp_sim_cosine"] = (model_results_small["exp_sim_cosine_forward"] + model_results_small["exp_sim_cosine_reverse"]) / 2

model_results_small["Rho"] = (model_results_small["Rho_forward"] + model_results_small["Rho_reverse"]) / 2

# PLot emb_sim_cos vs Rho
plt.figure() 
plt.scatter(model_results_small['emb_sim_cos'], model_results_small['Rho'], s=2, alpha=0.5, c = 'coral')
plt.xlabel('Emb Sim Cos')
plt.ylabel('Rho')
plt.title('Emb Sim Cos vs Rho')

# Same plot but different color for label == P/N and U
plt.figure()
plt.scatter(model_results_small[model_results_small['label_forward'] == 'P']['emb_sim_cos'], model_results_small[model_results_small['label_forward'] == 'P']['Rho'], s=2, alpha=0.5, c = 'green')
plt.scatter(model_results_small[model_results_small['label_forward'] == 'N']['emb_sim_cos'], model_results_small[model_results_small['label_forward'] == 'N']['Rho'], s=2, alpha=0.5, c = 'green')
plt.scatter(model_results_small[model_results_small['label_forward'] == 'U']['emb_sim_cos'], model_results_small[model_results_small['label_forward'] == 'U']['Rho'], s=2, alpha=0.5, c = 'blue')
plt.xlabel('Emb Sim Cos')
plt.ylabel('Rho')
plt.title('Emb Sim Cos vs Rho')

# COmpute correlation between emb_sim_cos and Rho
correlation = model_results_small[model_results_small['label_forward'] == 'U'][['emb_sim_cos', 'Rho']].corr().iloc[0, 1]
print(f"Correlation between emb_sim_cos and Rho: {correlation}")



# Plot exp_sim_cosine vs Rho
plt.figure()
plt.scatter(model_results_small['exp_sim_cosine'], model_results_small['Rho'], s=2, alpha=0.5, c = 'coral')
plt.xlabel('Exp Sim Cosine')
plt.ylabel('Rho')
plt.title('Exp Sim Cosine vs Rho')

# Plot exp_sim_cosine vs Rho
plt.figure()
plt.scatter(model_results_small['exp_sim_cosine_forward'], model_results_small['Rho_forward'], s=2, alpha=0.5, c = 'coral')
plt.scatter(model_results_small['exp_sim_cosine_reverse'], model_results_small['Rho_forward'], s=2, alpha=0.5, c = 'coral')
plt.xlabel('Exp Sim Cosine')
plt.ylabel('Rho')
plt.title('Exp Sim Cosine vs Rho')


# Filter to keep the rows with odd index
model_results_small = model_results_small.iloc[::2, :]

# Merge the correlation data with the combined similarity data on Gene1 and Gene2
merged_df = combined_df.merge(correlation_df, left_on=["gene_name_human", "gene_name_mouse"], right_on=["Gene1", "Gene2"], how="inner")

# Plot domain identity vs Pearson correlation
import matplotlib.pyplot as plt


orthologs_true = merged_df[merged_df['orthologs'] == True]
orthologs_false = merged_df[merged_df['orthologs'] == False]


# Plot promoter similarity against expression similarity (scatter plot)
plt.figure()
plt.scatter(orthologs_false['promoter_identity'], orthologs_false['Pearson_Correlation_R2'], s=2, alpha=0.5, c = 'coral')
plt.scatter(orthologs_true['promoter_identity'], orthologs_true['Pearson_Correlation_R2'], s=3, alpha=1, c='blue')
plt.legend(['Orthologs = False', 'Orthologs = True'])
plt.xlabel('Promoter Similarity')
plt.ylabel('Expression Similarity (cosine)')
plt.title('Promoter Similarity vs Expression Similarity')


# Plot promoter similarity against expression similarity (scatter plot)
plt.figure()
plt.scatter(orthologs_false['promoter_identity'], orthologs_false['Slope'], s=2, alpha=0.5, c = 'coral')
plt.scatter(orthologs_true['promoter_identity'], orthologs_true['Slope'], s=3, alpha=1, c='blue')
plt.legend(['Orthologs = False', 'Orthologs = True'])
plt.xlabel('Promoter Similarity')
plt.ylabel('Expression Similarity (cosine)')
plt.title('Promoter Similarity vs Expression Similarity')

#rho(h,m)=2cov(h,m)/(var(h)+var(m))

rho_data = []

for i in range(len(merged_log_ratios)):
    human_log_ratios = merged_log_ratios.iloc[i, 2:7]  # Skip Gene1 and Gene2 columns
    mouse_log_ratios = merged_log_ratios.iloc[i, 7:]  # Skip Gene1 and Gene2 columns
    gene1 = merged_log_ratios.iloc[i]["Gene1"]
    gene2 = merged_log_ratios.iloc[i]["Gene2"]

    human_values = pd.to_numeric(human_log_ratios)  # Skip Gene1 and Gene2 columns
    mouse_values = pd.to_numeric(mouse_log_ratios) # Skip Gene1 and Gene2 columns
    
    # Calculate Pearson correlation and slope
    covariance = np.cov(human_values, mouse_values)[0][1]
    var_human = np.var(human_values)
    var_mouse = np.var(mouse_values)
    rho = 2 * covariance / (var_human + var_mouse)

    rho_data.append({
        "Gene1": gene1,
        "Gene2": gene2,
        "Rho": rho
    })

# Create a new DataFrame with the results
rho_df = pd.DataFrame(rho_data)

# Add one column with the gene
rho_df = rho_df.rename(columns={"Gene1": "Gene_A", "Gene2": "Gene_B"})
# Add one column with the Gene_A_id_human and another one with the Gene_A_id_mouse
rho_df["Gene_A_id_human"] = rho_df["Gene_A"].map(orthologs_mapped.set_index("Gene name")["Gene stable ID"])
rho_df["Gene_A_id_mouse"] = rho_df["Gene_A"].map(orthologs_mapped.set_index("Gene name")["Mouse gene stable ID"])

# Add one column with the Gene_B_id_human and another one with the Gene_B_id_mouse
rho_df["Gene_B_id_human"] = rho_df["Gene_B"].map(orthologs_mapped.set_index("Gene name")["Gene stable ID"])
rho_df["Gene_B_id_mouse"] = rho_df["Gene_B"].map(orthologs_mapped.set_index("Gene name")["Mouse gene stable ID"])


# Sort the data frame by rho values
rho_df = rho_df.sort_values(by="Rho", ascending= False)

rho_df['label'] = "U"
# Write P for the top 3000 pairs
# Write N for the bottom 3000 pairs
rho_df.loc[rho_df.index[:5000], "label"] = "P"
rho_df.loc[rho_df.index[-5000:], "label"] = "N"

# Save to a file
rho_df.to_csv("results/pairs_labeled/kinase_vlr_com_similarity_rho_rank_3000.tsv", sep="\t", index=False)






# New df with the Gene 1 and Gene 2 columns inverted
rho_df_inverted = rho_df[['Gene2', 'Gene1', 'Rho']]
# Rename the columns of the data frame to reference gene and ortholog
rho_df_inverted = rho_df_inverted.rename(columns={"Gene2": "Gene1", "Gene1": "Gene2"})


# Combine the two data frames
rho_df_combined = pd.concat([rho_df, rho_df_inverted], ignore_index=True)


    
# Merge the rho data with the combined similarity data on Gene1 and Gene2
merged_df = combined_df.merge(rho_df_combined, left_on=["gene_name_human", "gene_name_mouse"], right_on=["Gene1", "Gene2"], how="inner")



orthologs_true = merged_df[merged_df['orthologs'] == True]
orthologs_false = merged_df[merged_df['orthologs'] == False]

# Plot promoter similarity against expression similarity (scatter plot)
plt.figure()
plt.scatter(orthologs_false['promoter_identity'], orthologs_false['Rho'], s=2, alpha=0.5, c = 'coral')
plt.scatter(orthologs_true['promoter_identity'], orthologs_true['Rho'], s=3, alpha=1, c='blue')
plt.legend(['Orthologs = False', 'Orthologs = True'])
plt.xlabel('Promoter Similarity')
plt.ylabel('Rho')
plt.title('Promoter Similarity vs Expression Similarity')
# Save the correlation data


correlation_df.to_csv("data/expression/human_mouse_correlation.tsv", sep="\t", index=False)






# Sort the data frame by rho values
sorted_df = merged_df.sort_values(by="Rho", ascending= False)
sorted_df['label'] = "U"
# Write P for the top 6000 pairs
# Write N for the bottom 6000 pairs

sorted_df.loc[sorted_df.index[:5000], "label"] = "P"
sorted_df.loc[sorted_df.index[-5000:], "label"] = "N"

# Save to a file
sorted_df.to_csv("results/pairs_labeled/kinase_vlr_similarity_rho_rank_5000.tsv", sep="\t", index=False)


# Save the top 5000 pairs to a list of tuples
top_6000_pairs = sorted_df.head(6000)
# Append to the list all pairs inverted
# top_6000_pairs_inverted = sorted_df.head(3000)[["Gene2", "Gene1"]].values.tolist()
# Save the bottom 5000 pairs to a list of tuples
bottom_6000_pairs = sorted_df.tail(6000)
# bottom_6000_pairs_inverted = sorted_df.tail(3000)[["Gene2", "Gene1"]].values.tolist()



# join the two lists
top_3000_pairs.extend(top_3000_pairs_inverted)
bottom_3000_pairs.extend(bottom_3000_pairs_inverted)


# Select first 5 rows of df
df = merged_log_ratios.iloc[:5, :]

# Save the results to files
human_log_ratios.to_csv("data/expression/human/human_log_ratios.tsv", sep="\t")
mouse_log_ratios.to_csv("data/expression/mouse/mouse_log_ratios.tsv", sep="\t")

# New column called label to the combined_df
merged_df["label"] = "U"
# Assign P if the pair is in the top 3000 pairs
for pair in top_6000_pairs:
    combined_df.loc[(combined_df["gene_name_human"] == pair[0]) & (combined_df["gene_name_mouse"] == pair[1]), "label"] = "P"
# Assign N if the pair is in the bottom 3000 pairs
for pair in bottom_6000_pairs:
    combined_df.loc[(combined_df["gene_name_human"] == pair[0]) & (combined_df["gene_name_mouse"] == pair[1]), "label"] = "N"