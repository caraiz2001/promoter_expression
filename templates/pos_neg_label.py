#!/usr/bin/env python3

# ---------------------------------------
# Import libraries
# ---------------------------------------
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


#%% 

# ---------------------------------------
# Define parameters
# ---------------------------------------
label_method = "${label_method}"
family_id = "${family_id}"
transform_method = "${transform_method}"
combined_df_file = "${combined_df}"
output_file = f"{family_id}_{transform_method}_similarity_{label_method}.tsv"

# ---------------------------------------
# Load the data
# ---------------------------------------
combined_df = pd.read_csv(combined_df_file, sep="\\t")

#%%
# ---------------------------------------
# Label the pairs
# ---------------------------------------

combined_df_labeled = combined_df.copy()

if label_method == "cos_pear_average_rank_5000":
    # Compute the average of the cosine and pearson correlation
    combined_df_labeled["average_sim"] = combined_df_labeled[["exp_sim_cosine", "exp_sim_pearson"]].mean(axis=1)
    # Sort the values by the average similarity
    combined_df_labeled = combined_df_labeled.sort_values("average_sim", ascending=False)
    # Label the first 5000 pairs as positive, the last 5000 as negative and the rest as unknown
    combined_df_labeled["label"] = "U"
    combined_df_labeled.loc[combined_df_labeled.index[:5000], "label"] = "P"
    combined_df_labeled.loc[combined_df_labeled.index[-5000:], "label"] = "N"

elif label_method == "cos_rank_5000":
    # Sort the values by the cosine similarity
    combined_df_labeled = combined_df_labeled.sort_values("exp_sim_cosine", ascending=False)
    # Label the first 5000 pairs as positive, the last 5000 as negative and the rest as unknown
    combined_df_labeled["label"] = "U"
    combined_df_labeled.loc[combined_df_labeled.index[:5000], "label"] = "P"
    combined_df_labeled.loc[combined_df_labeled.index[-5000:], "label"] = "N"


elif label_method == "vlr_5000":
    # compute the variance of the log ratio
    print("Computing the variance of the log ratio")

elif label_method == "cos_pear_threshold":
    if transform_method == "clr":
        cosine_upper = 0.94
        cosine_lower = 0.25
        pearson_upper = 0.5
        pearson_lower = 0.25
    else:
        cosine_upper = 0.94
        cosine_lower = 0.5
        pearson_upper = 0.5
        pearson_lower = 0.25
    combined_df_labeled['label'] = "U"
    # Set "P" to those pairs with exp_sim_cosine > 0.94 AND exp_sim_pearson > 0.5
    combined_df_labeled.loc[(combined_df_labeled['exp_sim_cosine'] > cosine_upper) & (combined_df_labeled['exp_sim_pearson'] > pearson_upper), 'label'] = "P"
    # Set "N" to those pairs with exp_sim_cosine < 0.5 AND exp_sim_pearson < 0.25
    combined_df_labeled.loc[(combined_df_labeled['exp_sim_cosine'] < cosine_lower) & (combined_df_labeled['exp_sim_pearson'] < pearson_lower), 'label'] = "N"

else:
    raise ValueError(f"Label method {label_method} not recognized. Please choose one of the following: cos_pear_average_rank_5000, cos_rank_5000, vlr_5000, cos_pear_threshold")
    

# %%
combined_df_labeled.to_csv(output_file, sep="\\t", index=False)

pos_pairs = combined_df_labeled[combined_df_labeled['label'] == 'P']
neg_pairs = combined_df_labeled[combined_df_labeled['label'] == 'N']
unknown_pairs = combined_df_labeled[combined_df_labeled['label'] == 'U']
# %%
# Same for promoter identity
plt.figure()
sns.kdeplot(pos_pairs['promoter_identity'], color = "coral", label = "positive pairs")
sns.kdeplot(neg_pairs['promoter_identity'], color = "blue", label = "negative pairs")
sns.kdeplot(unknown_pairs['promoter_identity'], color = "green", label = "unknown pairs")
plt.xlabel('Promoter identity')
plt.legend([f'Positive pair {len(pos_pairs)}', f'Negative pairs {len(neg_pairs)}', f'Unknown pairs {len(unknown_pairs)}'])
plt.title('Promoter identity')

plt.savefig(f'{family_id}_{transform_method}_{label_method}_promoter_identity_density.png', dpi=300, bbox_inches='tight')

plt.figure()
sns.kdeplot(pos_pairs['exp_sim_cosine'], color = "coral", label = "positive pairs")
sns.kdeplot(neg_pairs['exp_sim_cosine'], color = "blue", label = "negative pairs")
sns.kdeplot(unknown_pairs['exp_sim_cosine'], color = "green", label = "unknown pairs")
plt.xlabel('Expression similarity (cosine)')
plt.legend([f'Positive pair {len(pos_pairs)}', f'Negative pairs {len(neg_pairs)}', f'Unknown pairs {len(unknown_pairs)}'])
#plt.legend(['Positive pairs', 'Negative pairs', 'Unknown pairs'])
plt.title('Expression similarity density (cosine)')

plt.savefig(f'{family_id}_{transform_method}_{label_method}_exp_similarity_density.png', dpi=300, bbox_inches='tight')



print(combined_df_labeled['label'].value_counts())
# %%

