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
import seaborn as sns
from scipy.stats import gaussian_kde
import numpy as np
from scipy.special import rel_entr
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

number_orthologs_true = len(orthologs_true)
number_orthologs_false = len(orthologs_false)

# Plot promoter similarity against expression similarity (scatter plot)
plt.figure()
plt.scatter(orthologs_false['promoter_identity'], orthologs_false['exp_sim_cosine'], s=2, alpha=0.5, c = 'coral')
plt.scatter(orthologs_true['promoter_identity'], orthologs_true['exp_sim_cosine'], s=3, alpha=1, c='blue')
plt.legend([f'Orthologs = False ({number_orthologs_false})', f'Orthologs = True ({number_orthologs_true})'])
plt.xlabel('Promoter Similarity')
plt.ylabel('Expression Similarity (cosine)')
plt.title('Promoter Similarity vs Expression Similarity')

plt.savefig(f'{family_id}_{transform_method}_promoter_vs_exp_sim_cosine.png', dpi=300, bbox_inches='tight')

# Plot promoter similarity against expression similarity (scatter plot)
plt.figure()
plt.scatter(orthologs_false['promoter_identity'], orthologs_false['exp_sim_pearson'], s=2, alpha=0.5, c = 'coral')
plt.scatter(orthologs_true['promoter_identity'], orthologs_true['exp_sim_pearson'], s=3, alpha=1, c='blue')
plt.legend([f'Orthologs = False ({number_orthologs_false})', f'Orthologs = True ({number_orthologs_true})'])
plt.xlabel('Promoter Similarity')
plt.ylabel('Expression Similarity (pearson)')
plt.title('Promoter Similarity vs Expression Similarity')

plt.savefig(f'{family_id}_{transform_method}_promoter_vs_exp_sim_pearson.png', dpi=300, bbox_inches='tight')

# Plot scatter + density plot combined
# Create JointGrid manually
g = sns.JointGrid(data=combined_df, x="promoter_identity", y="exp_sim_cosine", height=8)

sns.scatterplot(data=orthologs_false, x="promoter_identity", y="exp_sim_cosine",
                color="coral", alpha=0.4, s=10, label= f"Orthologs = False ({number_orthologs_false}", ax=g.ax_joint, edgecolor = 'white')

sns.scatterplot(data=orthologs_true, x="promoter_identity", y="exp_sim_cosine",
                color="blue", alpha=0.7, s=15, label= f"Orthologs = True ({number_orthologs_true})", ax=g.ax_joint)

# KDEs — normalized *independently*
sns.kdeplot(data=orthologs_false, x="promoter_identity", ax=g.ax_marg_x,
            fill=True, color="coral", common_norm=False)

sns.kdeplot(data=orthologs_true, x="promoter_identity", ax=g.ax_marg_x,
            fill=True, color="blue", common_norm=False)

sns.kdeplot(data=orthologs_false, y="exp_sim_cosine", ax=g.ax_marg_y,
            fill=True, color="coral", common_norm=False)

sns.kdeplot(data=orthologs_true, y="exp_sim_cosine", ax=g.ax_marg_y,
            fill=True, color="blue", common_norm=False)

# Labels & Legend
g.ax_joint.set_xlabel("Promoter Similarity")
g.ax_joint.set_ylabel("Expression Similarity (cosine)")
g.ax_joint.legend(loc="upper right")

plt.suptitle("Promoter vs Expression Similarity with density plot", fontsize=16)
plt.tight_layout()
plt.subplots_adjust(top=0.92)
plt.savefig(f'{family_id}_{transform_method}_promoter_vs_exp_sim_cosine_jointgrid.png', dpi=300, bbox_inches='tight')


#%%
# Compute some statistics to compare both distributions (orthologs vs non-orthologs)
# For expression similarity (cosine)

# mean and std and iqr
orthologs_false_mean = orthologs_false['exp_sim_cosine'].mean()
orthologs_false_std = orthologs_false['exp_sim_cosine'].std()
orthologs_true_std = orthologs_true['exp_sim_cosine'].std()
orthologs_true_mean = orthologs_true['exp_sim_cosine'].mean()
# Interquartile range
orhtologs_true_iqr = orthologs_true["promoter_identity"].quantile(0.75) - orthologs_true["promoter_identity"].quantile(0.25)
orthologs_false_iqr = orthologs_false["promoter_identity"].quantile(0.75) - orthologs_false["promoter_identity"].quantile(0.25)

# median
orthologs_false_median = orthologs_false['exp_sim_cosine'].median()
orthologs_true_median = orthologs_true['exp_sim_cosine'].median()

# mode: most frquent discrete value (after rounding to 3 decimals), essentially mode of a histogram with a fixed bin size
orthologs_false_mode = orthologs_false['exp_sim_cosine'].mode(3)[0]
orthologs_true_mode = orthologs_true['exp_sim_cosine'].mode(3)[0]

# mode kde: peak of the estimated probability density function. To quantify where each distribution is centered most densely.
# Compute KDE and mode for the y-axis (Expression Similarity)
y_vals = np.linspace(-1.1, 1.1, 500)

kde_true_y = gaussian_kde(orthologs_true["exp_sim_cosine"])
density_true_y = kde_true_y(y_vals)
orthologs_true_kdemode = y_vals[np.argmax(density_true_y)]

kde_false_y = gaussian_kde(orthologs_false["exp_sim_cosine"])
density_false_y = kde_false_y(y_vals)
orthologs_false_kdemode = y_vals[np.argmax(density_false_y)]


# KL Divergence: measure of how one probability distribution diverges from a second, expected probability distribution.
p = density_true_y / density_true_y.sum()
q = density_false_y / density_false_y.sum()

kl_div = np.sum(rel_entr(p, q))
print(f"KL Divergence (True || False): {kl_div:.4f}")

# Cohen's d: effect size to quantify the difference between the two distributions
pooled_std = np.sqrt((orthologs_true_std**2 + orthologs_false_std**2) / 2)

cohens_d = (orthologs_true_mean - orthologs_false_mean) / pooled_std

# Print the results and write to a file (with the description of each statistic)
stats_file = f"{family_id}_{transform_method}_expression_similarity_stats.txt"
text = f"""
Expression Similarity (cosine) Statistics:
Orthologs = False:
- Mean: {orthologs_false_mean:.3f}
- Std: {orthologs_false_std:.3f}
- Median: {orthologs_false_median:.3f}
- Mode: {orthologs_false_mode:.3f}
- IQR: {orthologs_false_iqr:.3f}
- Mode KDE: {orthologs_false_kdemode:.3f}


Orthologs = True:
- Mean: {orthologs_true_mean:.3f}
- Std: {orthologs_true_std:.3f}
- Median: {orthologs_true_median:.3f}
- Mode: {orthologs_true_mode:.3f}
- IQR: {orhtologs_true_iqr:.3f}
- Mode KDE: {orthologs_true_kdemode:.3f}

Cohen's d: {cohens_d:.3f}
KL Divergence (True || False): {kl_div:.4f}
"""

print(text)
with open(stats_file, "w") as f:
    f.write(text)

print(f"Report written to {stats_file}")