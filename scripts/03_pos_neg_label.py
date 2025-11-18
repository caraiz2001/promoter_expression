# This script is used to label pairs are positive and negative for training
# It can be done using different criteria
# Input is combined similarity df 
# Output is the same df with an extra column "label" with values P, N or U


#%%
# ---------------------------------------
# Import libraries
# ---------------------------------------
import pandas as pd
import os


#%% 

# ---------------------------------------
# Define parameters
# ---------------------------------------
label_method = "cos_rank_5000"
family_id = "kinase"
transform_method = "clr"
combined_df_file = f"./results/combined_similarity/{family_id}_{transform_method}_combined_similarity.tsv"
output_file = f"./results/pairs_labeled/{family_id}_{transform_method}_similarity_{label_method}.tsv"
path = "/Users/crsitina/Documents/phd/promoter_expression"
os.chdir(path)

# ---------------------------------------
# Load the data
# ---------------------------------------
combined_df = pd.read_csv(combined_df_file, sep="\t")

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
    combined_df_labeled['label'] = "U"
    # Set "P" to those pairs with exp_sim_cosine > 0.94 AND exp_sim_pearson > 0.5
    combined_df_labeled.loc[(combined_df_labeled['exp_sim_cosine'] > 0.94) & (combined_df_labeled['exp_sim_pearson'] > 0.5), 'label'] = "P"
    # Set "N" to those pairs with exp_sim_cosine < 0.5 AND exp_sim_pearson < 0.25
    combined_df_labeled.loc[(combined_df_labeled['exp_sim_cosine'] < 0.5) & (combined_df_labeled['exp_sim_pearson'] < 0.25), 'label'] = "N"

else:
    raise ValueError(f"Label method {label_method} not recognized. Please choose one of the following: cos_pear_average_rank_5000, cos_rank_5000, vlr_5000, cos_pear_threshold")
    

# %%
combined_df_labeled.to_csv(output_file, sep="\t", index=False)
# %%
print(combined_df_labeled['label'].value_counts())
# %%
