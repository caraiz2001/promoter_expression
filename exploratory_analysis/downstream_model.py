#%%

# ------------------------------------
# Import libraries
# ------------------------------------
import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression
from scipy.stats import kendalltau
from scipy.stats import pearsonr
import seaborn as sns
from sklearn.metrics import mean_squared_error, r2_score


#%%
# --------------------
# Define parameters
# --------------------

family_id = "kinase"
transform_method = "ratios"
label_method = "var_rank_5000"

# Load orthologs df
orthologs_df = pd.read_csv("./data/db_ids/human_mouse_one2one_filtered_mapped.txt", sep="\t")
folder_id = f"{family_id}_{transform_method}"

model_results_file = f"./results/model_results/{folder_id}/all_pairs_all_metrics_{transform_method}_{label_method}.tsv"
set_id_pairs_file = f"./results/model_results/{folder_id}/{family_id}_{transform_method}_id_pairs_set_{label_method}.tsv"
path = '/Users/crsitina/Documents/phd/promoter_expression/'
print(path)
os.chdir(path)

# Load the Rho data


#%% 

# ------------------------------------
# Define functions
# ------------------------------------

def analyze_intervals(dataframe, reference_column_exp, reference_column_sequence):
    # interval_df = dataframe.copy()
    interval_df = dataframe.sort_values(by=reference_column_sequence, ascending=False)
    # Ensure the reference column exists in the dataframe
    if reference_column_exp not in interval_df.columns or reference_column_sequence not in interval_df.columns:
        raise ValueError("The specified reference columns do not exist in the dataframe.")
    else:
        print(f"Reference columns {reference_column_exp} and {reference_column_sequence} exist in the dataframe.")
    # Initialize a list to store the results
    results = []
    
    # Calculate the range of the reference_column_sequence
    min_value = interval_df[reference_column_sequence].min().round(2)
    max_value = interval_df[reference_column_sequence].max().round(2)

    steps = max_value / 40

    print(f"Min value: {min_value}, Max value: {max_value}, Steps: {steps}")

    # Iterate over the intervals
    for i in np.arange(max_value, min_value, -steps):
        print(f"Interval: {i:.2f}-{i-steps:.2f}")
        # Define the interval
        df = interval_df[(interval_df[reference_column_sequence] <= i) & (interval_df[reference_column_sequence] > i - steps)]
        print(df.shape)
        # Calculate the number of pairs (rows) in the interval
        num_pairs = df.shape[0]
        print(num_pairs)

        # Extract maximum and minimum values of the reference column
        max_value = df[reference_column_exp].max()
        min_value = df[reference_column_exp].min()

        # Difference between maximum and minimum
        diff_value = max_value - min_value

        # Extract 1st and 3rd quartile of the reference column
        q1_value = df[reference_column_exp].quantile(0.05)
        q3_value = df[reference_column_exp].quantile(0.95)

        # Compute difference of q3 and q1
        iqr_value = q3_value - q1_value

        # Append the results to the list
        results.append({
            f"Interval of {reference_column_sequence}": f"{i:.2f}-{i-steps:.2f}",
            "num_pairs": num_pairs,
            f"max_{reference_column_exp}": max_value,
            f"min_{reference_column_exp}": min_value,
            f"diff_{reference_column_exp}": diff_value,
            f"q1_{reference_column_exp}": q1_value,
            f"q3_{reference_column_exp}": q3_value,
            f"iqr_{reference_column_exp}": iqr_value
        })

        # Convert the results to a dataframe
        results_df = pd.DataFrame(results)

        # Round the values to 2 decimal places
        results_df[f"max_{reference_column_exp}"] = results_df[f"max_{reference_column_exp}"].round(2)
        results_df[f"min_{reference_column_exp}"] = results_df[f"min_{reference_column_exp}"].round(2)
        results_df[f"q1_{reference_column_exp}"] = results_df[f"q1_{reference_column_exp}"].round(2)
        results_df[f"q3_{reference_column_exp}"] = results_df[f"q3_{reference_column_exp}"].round(2)
        results_df[f"iqr_{reference_column_exp}"] = results_df[f"iqr_{reference_column_exp}"].round(2)
        results_df["num_pairs"] = results_df["num_pairs"].astype(int)
    
    print(f"Number of intervals: {len(results)}")
    return results_df

def calculate_statistics(results_df, reference_column, reference_column_sequence):
    # Average of the "diff_{reference_column}"
    average_diff = results_df[f"diff_{reference_column}"].mean()
    above_average_intervals = results_df[results_df[f"diff_{reference_column}"] > average_diff].shape[0]
    below_average_intervals = results_df[results_df[f"diff_{reference_column}"] < average_diff].shape[0]
    above_average_intervals_df = results_df[results_df[f"diff_{reference_column}"] > average_diff]
    below_average_intervals_df = results_df[results_df[f"diff_{reference_column}"] < average_diff]
    total_pairs_informative = below_average_intervals_df["num_pairs"].sum()
    total_pairs = results_df["num_pairs"].sum()
    percentage = round((total_pairs_informative / total_pairs) * 100, 2)

    # Average of the "iqr_{reference_column}"
    average_iqr = results_df[f"iqr_{reference_column}"].mean()
    above_average_intervals_iqr = results_df[results_df[f"iqr_{reference_column}"] > average_iqr].shape[0]
    below_average_intervals_iqr = results_df[results_df[f"iqr_{reference_column}"] < average_iqr].shape[0]
    above_average_intervals_iqr_df = results_df[results_df[f"iqr_{reference_column}"] > average_iqr]
    below_average_intervals_iqr_df = results_df[results_df[f"iqr_{reference_column}"] < average_iqr]
    total_pairs_informative_iqr = below_average_intervals_iqr_df["num_pairs"].sum()
    percentage_iqr = round((total_pairs_informative_iqr / total_pairs) * 100, 2)

    # Keep only first and last rows of the dataframe
    above_average_intervals_df = above_average_intervals_df.iloc[[0, -1]]
    below_average_intervals_df = below_average_intervals_df.iloc[[0, -1]]
    above_average_intervals_iqr_df = above_average_intervals_iqr_df.iloc[[0, -1]]
    below_average_intervals_iqr_df = below_average_intervals_iqr_df.iloc[[0, -1]]

    # Format the text
    text = f"""

    MINIMUM AND MAXIMUM {reference_column_sequence}:

    - Average difference of {reference_column}: {average_diff:.2f}
    - Number of intervals above the average difference (non-informative): {above_average_intervals} out of {results_df.shape[0]} intervals
    - Number of intervals below the average difference (informative): {below_average_intervals} out of {results_df.shape[0]} intervals
    => Total number of pairs in the informative intervals: {total_pairs_informative} out of {total_pairs} ({percentage}%)

    Non-informative intervals:
    {above_average_intervals_df.iloc[:, 0:2]}

    Informative intervals:
    {below_average_intervals_df.iloc[:, 0:2]}

    INTERQUARTILE RANGE (0.05 and 0.95):
    - Average interquartile range of {reference_column}: {average_iqr:.2f}
    - Number of intervals above the average interquartile range (non-informative): {above_average_intervals_iqr} out of {results_df.shape[0]} intervals
    - Number of intervals below the average interquartile range (informative): {below_average_intervals_iqr} out of {results_df.shape[0]} intervals
    => Total number of pairs in the informative intervals (interquartile range): {total_pairs_informative_iqr} out of {total_pairs} ({percentage_iqr}%)

    Non-informative intervals:
    {above_average_intervals_iqr_df.iloc[:, 0:3]}

    Informative intervals:
    {below_average_intervals_iqr_df.iloc[:, 0:3]}
    """
    print(text)
    return text
#%%

# ------------------------
# Subset and process the data
# ------------------------

df = pd.read_csv(model_results_file, sep="\t", index_col=0)
set_id_pairs = pd.read_csv(set_id_pairs_file, sep="\t", index_col=0)

# Merge the set_id_pairs to df by gene_id_human and gene_id_mouse, keeping all rows from df
df_merged = df.merge(set_id_pairs, on=["gene_id_human", "gene_id_mouse"], how="left")


print(f"Number of pairs in the model results: {df.shape}")
print(f"Number of pairs in the set_id_pairs: {set_id_pairs.shape}")

# Subset the set_id_pairs into training, validation and test sets
train_pairs = set_id_pairs[set_id_pairs['model'] == 'training']
validation_pairs = set_id_pairs[set_id_pairs['model'] == 'validation']
test_pairs = set_id_pairs[set_id_pairs['model'] == 'test']

print(f"Number of training pairs {train_pairs.shape}")
print(f"Number of validation pairs {validation_pairs.shape}")
print(f"Number of test pairs {test_pairs.shape}")

# Extract a list of the pairs that are in the train set (gene_id_human, gene_id_mouse)
train_pairs_list = list(zip(train_pairs['gene_id_human'], train_pairs['gene_id_mouse']))

# Add a column to df that is True if the pair is in the train set and False otherwise
df['train'] = df.apply(lambda x: (x['gene_id_human'], x['gene_id_mouse']) in train_pairs_list, axis = 1)

# Divide the data frame into train and the rest
train_results = df[df['train']] 
print(train_results['label'].value_counts())
notrain_results = df[~df['train']]
print(notrain_results['label'].value_counts())
print(train_results.shape, notrain_results.shape)

# The test results, divide them into positive/negative and unknown (label column)
pos_neg = notrain_results[notrain_results['label'] != 'U'] # validation and test set
unknown = notrain_results[notrain_results['label'] == 'U']
print(pos_neg.shape, unknown.shape)

# %%

# ------------------------
# Plots
# ------------------------
exp_column = "variance"
#exp_column = "exp_sim_cosine"

# Plot a scatter plot of emb_sim_cosine against exp_sim_cosing for the train and test sets
plt.figure()
plt.scatter(unknown['emb_sim_cos'], - unknown[exp_column], s=2, alpha=0.5, c = 'blue')
plt.scatter(pos_neg['emb_sim_cos'], - pos_neg[exp_column], s=2, alpha=1, c = 'green')
plt.legend(['Unknown', 'Positive/Negative'])
plt.xlabel('Embedding similarity (cosine)')
plt.ylabel(f'Expression similarity {exp_column}')
plt.title(f'Embedding similarity vs Expression similarity ({exp_column})')
plt.savefig(f"./exploratory_analysis/plots/{folder_id}/{label_method}/emb_sim_vs_{exp_column}_{label_method}.png")


plt.figure()
plt.scatter(unknown['promoter_identity'], unknown[exp_column], s=2, alpha=0.5, c = 'blue')
plt.scatter(pos_neg['promoter_identity'], pos_neg[exp_column], s=2, alpha=1, c = 'green')
plt.legend(['Unknown', 'Positive/Negative'])
plt.xlabel('Promoter identity')
plt.ylabel(f'Expression similarity ({exp_column})')
plt.title(f'Promoter identity vs Expression similarity ({exp_column})')
plt.savefig(f"./exploratory_analysis/plots/{folder_id}/{label_method}/promoter_identity_vs_{exp_column}_{label_method}.png")


# # Same plot but log-transformed of the expression similarity
# plt.figure()
# plt.scatter(unknown['emb_sim_cos'], np.log(unknown[exp_column] + 0.1), s=2, alpha=0.5, c = 'blue')
# plt.scatter(pos_neg['emb_sim_cos'], np.log(pos_neg[exp_column] + 0.1), s=2, alpha=1, c = 'green')
# plt.legend(['Unknown', 'Positive/Negative'])
# plt.xlabel('Embedding similarity (cosine)')
# plt.ylabel('Log expression similarity (cosine)')
# plt.title('Embedding similarity vs Log expression similarity (cosine)')
# plt.savefig(f"./exploratory_analysis/plots/{folder_id}/{label_method}/emb_sim_vs_log_{exp_column}_{label_method}.png")
# %%

# -----------------------
# Compute Performance metrics
# -----------------------

# ------- PEARSON CORRELATION -------
corr = np.corrcoef(notrain_results['emb_sim_cos'], notrain_results[exp_column])
print(f"Correlation of embedding similarity and expression similarity (Cosine) for the non-train set: {corr[0][1]}")
# Compute correlation of emb_sim_cos and exp_sim_cosine for the train set
corr_train = np.corrcoef(train_results['emb_sim_cos'], train_results[exp_column])
print(f"Correlation of embedding similarity and expression similarity (Cosine) for the train set: {corr_train[0][1]}")


# Compute Pearson correlation of the log transformed data
corr_log = np.corrcoef(np.log(notrain_results['emb_sim_cos'] + 1), np.log(notrain_results[exp_column] + 1))
print(f"Correlation of log(embedding similarity) and log(expression similarity) for the non-train set: {corr_log[0][1]}")
# Compute correlation of log(emb_sim_cos) and log(exp_sim_cosine) for the train set
corr_log_train = np.corrcoef(np.log(train_results['emb_sim_cos'] + 1), np.log(train_results[exp_column] + 1))


# ------- SPEARMAN CORRELATION -------
# Compute Spearman correlation for the non-train set
spearman_corr, _ = spearmanr(notrain_results['emb_sim_cos'], notrain_results[exp_column])
print(f"Spearman correlation of embedding similarity and expression similarity (Cosine) for the non-train set: {spearman_corr}")

# Compute Spearman correlation for the train set
spearman_corr_train, _ = spearmanr(train_results['emb_sim_cos'], train_results[exp_column])
print(f"Spearman correlation of embedding similarity and expression similarity (Cosine) for the train set: {spearman_corr_train}")

# identity_stats = calculate_statistics(analyze_intervals(notrain_results, exp_column, 'promoter_identity'), exp_column, 'promoter_identity')
# promoter_stats = calculate_statistics(analyze_intervals(notrain_results, exp_column, 'emb_sim_cos'), exp_column, 'emb_sim_cos')



# -------- KENDALL CORRELATION --------
# Compute Kendall correlation for the non-train set
kendall_corr, _ = kendalltau(notrain_results['emb_sim_cos'], notrain_results[exp_column])
kendall_corr_train, _ = kendalltau(train_results['emb_sim_cos'], train_results[exp_column])
print(f"Kendall correlation of embedding similarity and expression similarity (Cosine) for the train set: {kendall_corr_train}")
print(f"Kendall correlation of embedding similarity and expression similarity (Cosine) for the non-train set: {kendall_corr}")

# %%
# -------- LINEAR REGRESSION --------
x = notrain_results['emb_sim_cos'].values
y = notrain_results[exp_column].values

# Reshape x for sklearn
x_reshaped = x.reshape(-1, 1)

# Fit linear regression
model = LinearRegression()
model.fit(x_reshaped, y)

slope = model.coef_[0]
intercept = model.intercept_
print(f"Slope: {slope}, Intercept: {intercept}")
# Compute R^2
r_squared = model.score(x_reshaped, y)
print(f"R^2: {r_squared}")


# Optional: Visualize
plt.figure()
plt.scatter(unknown['emb_sim_cos'], unknown[exp_column], s=2, alpha=0.5, c = 'blue' , label = 'Not labeled')
plt.scatter(pos_neg['emb_sim_cos'], pos_neg[exp_column], s=2, alpha=1, c = 'green', label = 'Negative/Positive')
# plt.scatter(x, y, alpha=0.5, label="Original data", s=2)
plt.plot(x, model.predict(x_reshaped), color='red', label="Linear fit")
plt.xlabel("Embedding similarity")
plt.ylabel("Expression similarity")
plt.legend()
# Add vertical line in -0.45
# plt.axvline(x=-0.45, color='black', linestyle='--', label="Threshold")
# plt.axvline(x=-0.6, color='black', linestyle='--', label="Threshold")
# # Add horizontal line in 0.5
# plt.axhline(y=0.5, color='black', linestyle='--', label="Threshold")
# plt.axhline(y=0.7, color='black', linestyle='--', label="Threshold")
plt.title("Linear regression on original data")
plt.savefig(f"./exploratory_analysis/plots/{folder_id}/{label_method}/emb_sim_vs_{exp_column}_{label_method}_linear_regression.png")


# --------- PAIRWISE AGREEMENT ----------
import random 

def pairwise_agreement_from_df(df, exp_col='exp_sim_cosine', sim_col='emb_sim_cos', num_samples=1000000):
    exp_sims = df[exp_col].values
    embed_sims = df[sim_col].values
    N = len(df)

    count_agree = 0
    total = 0

    for _ in range(num_samples):
        i1, i2 = random.sample(range(N), 2)

        exp_diff = exp_sims[i1] - exp_sims[i2]
        sim_diff = embed_sims[i1] - embed_sims[i2]

        # Compare the sign of differences
        if np.sign(exp_diff) == np.sign(sim_diff):
            count_agree += 1
        total += 1

    return count_agree / total

def find_twilight_threshold(df, exp_col='exp_sim_cosine', sim_col='emb_sim_cos', agreement_target=0.9):
    thresholds = np.linspace(df[sim_col].min(), df[sim_col].max(), 100)
    
    for t in thresholds:
        # print(t)
        subset = df[df[sim_col] >= t]
        if len(subset) < 10:
            continue  # not enough data for reliable agreement
        score = pairwise_agreement_from_df(subset, exp_col, sim_col, num_samples=100000)
        # print(score)
        if score >= agreement_target:
            return t
    return None

def get_agreement_curve(df, exp_col='exp_sim_cosine', sim_col='emb_sim_cos', min_points=10, num_samples=100000):
    import pandas as pd
    import numpy as np

    thresholds = np.linspace(df[sim_col].min(), df[sim_col].max(), 100)
    results = []

    for t in thresholds:
        subset = df[df[sim_col] >= t]
        if len(subset) < min_points:
            continue  # skip small subsets
        score = pairwise_agreement_from_df(subset, exp_col, sim_col, num_samples=num_samples)
        results.append({'threshold': t, 'agreement': score})

    return pd.DataFrame(results)

# %%

score = pairwise_agreement_from_df(notrain_results,exp_col=exp_column, sim_col='emb_sim_cos')
print(f"Agreement score: {score:.4f}")
agreement_target = 0.85
threshold = find_twilight_threshold(notrain_results, exp_col=exp_column, sim_col='emb_sim_cos', agreement_target = agreement_target)
print(f"Twilight zone threshold ({agreement_target} agreement): {threshold}")
results_df = get_agreement_curve(notrain_results)

# Plot the agreement curve
plt.figure()
plt.plot(results_df['threshold'], results_df['agreement'])
plt.xlabel('Threshold')
plt.ylabel('Agreement')
plt.title('Agreement Curve')
plt.axhline(y=agreement_target, color='r', linestyle='--', label=f'Agreement target ({agreement_target})')
plt.legend()
plt.savefig(f"./exploratory_analysis/plots/{folder_id}/{label_method}/agreement_curve_{label_method}.png")


results_promoter_identity = get_agreement_curve(notrain_results, exp_col='exp_sim_cosine', sim_col='promoter_identity', min_points=10, num_samples=100000)

plt.figure()
plt.plot(results_promoter_identity['threshold'], results_promoter_identity['agreement'])
plt.xlabel('Threshold')
plt.ylabel('Agreement')
plt.title('Agreement Curve (Promoter Identity)')
plt.axhline(y=agreement_target, color='r', linestyle='--', label=f'Agreement target ({agreement_target})')
plt.legend()
plt.savefig(f"./exploratory_analysis/plots/{folder_id}/{label_method}/agreement_curve_promoter_identity_{label_method}.png")

score_identity = pairwise_agreement_from_df(notrain_results, exp_col='exp_sim_cosine', sim_col='promoter_identity')
kendall_identity = kendalltau(notrain_results['promoter_identity'], notrain_results[exp_column])



# ----------- FALSE POSITIVE RANKING -----------
# Sort the dataframe by the embedding similarity
sorted_df = notrain_results.sort_values(by='emb_sim_cos', ascending=False)

# Get the index of the first pair labeled as Negative
first_negative_index = sorted_df[sorted_df['label'] == 'N'].index[0]
print(f"First negative index: {first_negative_index}")

# Get the embedding similarity of the first negative pair
first_negative_similarity = sorted_df.loc[first_negative_index, 'emb_sim_cos']
print(f"First negative similarity: {first_negative_similarity}")
#%%

# ------------------------
# WRITE REPORT
# ------------------------

# Write correlation values to a text file
with open(f"./results/downstream_reports/{folder_id}_{label_method}_report.txt", "w") as f:
    f.write(f"Expression similarity column: {exp_column}\n")
    f.write("NUMBER OF PAIRS: \n")
    f.write(f"Number of training pairs {train_pairs.shape}\n")
    f.write(f"Number of validation pairs {validation_pairs.shape}\n")
    f.write(f"Number of test pairs {test_pairs.shape}\n")

    f.write("CORRELATIONS: \n")
    f.write(f"Correlation of embedding similarity and expression similarity (Cosine) for the non-train set: {corr[0][1]}\n")
    f.write(f"Correlation of embedding similarity and expression similarity (Cosine) for the train set: {corr_train[0][1]}\n")

    
    f.write(f"Spearman correlation of embedding similarity and expression similarity (Cosine) for the non-train set: {spearman_corr}\n")
    f.write(f"Spearman correlation of embedding similarity and expression similarity (Cosine) for the train set: {spearman_corr_train}\n")

    f.write("LINEAR REGRESSION: \n")
    f.write(f"Linear regression slope: {slope}\n")
    f.write(f"Linear regression intercept: {intercept}\n")
    f.write(f"Linear regression R^2: {model.score(x_reshaped, y)}\n")

    f.write("KENDALL CORRELATION: \n")
    f.write(f"Kendall correlation of embedding similarity and expression similarity (Cosine) for the train set: {kendall_corr_train}\n")
    f.write(f"Kendall correlation of embedding similarity and expression similarity (Cosine) for the non-train set: {kendall_corr}\n")

    f.write("AGREEMENT SCORE: \n")
    f.write(f"Agreement score: {score:.4f}\n")
    f.write(f"Twilight zone threshold ({agreement_target} agreement): {threshold}\n")


# with open(f"./results/downstream_reports/{folder_id}_{label_method}_report.txt", "a") as text_file:
#     text_file.write(identity_stats)


# with open(f"./results/downstream_reports/{folder_id}_{label_method}_report.txt", "a") as text_file:
#     text_file.write(promoter_stats)




# ----- OTHER ANALYSIS -------

# Select which pairs have emb sim (-0.4, -0.6) and emb sim 0.4, 0.7
weird_pairs = notrain_results[(notrain_results['emb_sim_cos'] < -0.60) & (notrain_results['emb_sim_cos'] > -0.7) & (notrain_results['exp_sim_cosine'] > 0.25) & (notrain_results['exp_sim_cosine'] < 0.5)]
print(f"Weird pairs: {weird_pairs.shape}")

# See how many times each gene_id appears in the weird pairs
weird_pairs_counts = pd.concat([weird_pairs['gene_name_human'], weird_pairs['gene_name_mouse']]).value_counts()

# Optional: Visualize
plt.figure()
plt.scatter(unknown['emb_sim_cos'], unknown[exp_column], s=2, alpha=0.5, c = 'blue' , label = 'Not labeled')
plt.scatter(pos_neg['emb_sim_cos'], pos_neg[exp_column], s=2, alpha=1, c = 'green', label = 'Negative/Positive')
# plt.scatter(x, y, alpha=0.5, label="Original data", s=2)
plt.plot(x, model.predict(x_reshaped), color='red', label="Linear fit")
plt.xlabel("Embedding similarity")
plt.ylabel("Expression similarity")
plt.legend()
# Add vertical line in -0.45
plt.axvline(x=-0.6, color='black', linestyle='--', label="Threshold")
plt.axvline(x=-0.7, color='black', linestyle='--', label="Threshold")
# Add horizontal line in 0.5
plt.axhline(y=0.25, color='black', linestyle='--', label="Threshold")
plt.axhline(y=0.5, color='black', linestyle='--', label="Threshold")
plt.title("Linear regression on original data")
plt.savefig(f"./exploratory_analysis/plots/{folder_id}/{label_method}/emb_sim_vs_{exp_column}_{label_method}_linear_regression.png")




# Select which pairs have emb sim (-0.4, -0.6) and emb sim 0.4, 0.7
normal_pairs = notrain_results[(notrain_results['emb_sim_cos'] > 0) & (notrain_results['emb_sim_cos'] < 0.15) & (notrain_results['exp_sim_cosine'] > 0.25) & (notrain_results['exp_sim_cosine'] < 0.5)]
print(f"Normal pairs: {normal_pairs.shape}")

# See how many times each gene_id appears in the weird pairs
normal_pairs_counts = pd.concat([normal_pairs['gene_name_human'], normal_pairs['gene_name_mouse']]).value_counts()

# 
len(pd.concat([normal_pairs['gene_name_human'], normal_pairs['gene_name_mouse']]).unique())

# Optional: Visualize
plt.figure()
plt.scatter(unknown['emb_sim_cos'], unknown[exp_column], s=2, alpha=0.5, c = 'blue' , label = 'Not labeled')
plt.scatter(pos_neg['emb_sim_cos'], pos_neg[exp_column], s=2, alpha=1, c = 'green', label = 'Negative/Positive')
# plt.scatter(x, y, alpha=0.5, label="Original data", s=2)
plt.plot(x, model.predict(x_reshaped), color='red', label="Linear fit")
plt.xlabel("Embedding similarity")
plt.ylabel("Expression similarity")
plt.legend()
# Add vertical line in -0.45
plt.axvline(x=0, color='black', linestyle='--', label="Threshold")
plt.axvline(x=0.15, color='black', linestyle='--', label="Threshold")
# Add horizontal line in 0.5
plt.axhline(y=0.25, color='black', linestyle='--', label="Threshold")
plt.axhline(y=0.5, color='black', linestyle='--', label="Threshold")
plt.title("Linear regression on original data")
plt.savefig(f"./exploratory_analysis/plots/{folder_id}/{label_method}/emb_sim_vs_{exp_column}_{label_method}_linear_regression.png")


# orthologs_df wich Gene name is SBK3
orthologs_df[orthologs_df['Gene name'] == 'SBK3']
orthologs_df[orthologs_df['Gene name'] == 'SBK2']


# %%
# Select merged_df which gene_name_human or gene_name_mouse is SBK3 or SBK2
merged_df_sbk= df_merged[(df_merged['gene_name_human'] == 'SBK3') | (df_merged['gene_name_mouse'] == 'SBK3') | (df_merged['gene_name_human'] == 'SBK2') | (df_merged['gene_name_mouse'] == 'SBK2')]
print(merged_df_sbk.shape)

merged_df_sbk2 = df_merged[(df_merged['gene_name_human'] == 'SBK2') | (df_merged['gene_name_mouse'] == 'SBK2')]
print(merged_df_sbk2.shape)

merged_df_mapk13 = df_merged[(df_merged['gene_name_human'] == 'MAPK13') | (df_merged['gene_name_mouse'] == 'MAPK13')]
print(merged_df_mapk13.shape)

merged_df_kit = df_merged[(df_merged['gene_name_human'] == 'KIT') | (df_merged['gene_name_mouse'] == 'KIT')]
print(merged_df_kit.shape)

merged_df_prkcq = df_merged[(df_merged['gene_name_human'] == 'PRKCQ') | (df_merged['gene_name_mouse'] == 'PRKCQ')]
print(merged_df_prkcq.shape)


# Plot a scatter plot
plt.figure()
plt.scatter(merged_df_sbk['emb_sim_cos'], merged_df_sbk[exp_column], s=2, alpha=0.5, c = 'red')
plt.scatter(merged_df_mapk13['emb_sim_cos'], merged_df_mapk13[exp_column], s=2, alpha=1, c = 'green')
plt.legend(['SBK3/SBK2', 'MAPK13'])
plt.xlabel('Embedding similarity (cosine)')
plt.ylabel(f'Expression similarity (cosine)')
plt.title(f'Embedding similarity vs Expression similarity (cosine)')

# 
plt.figure()
plt.scatter(merged_df_prkcq['emb_sim_cos'], merged_df_prkcq[exp_column], s=2, alpha=0.5, c = 'red')
plt.scatter(merged_df_mapk13['emb_sim_cos'], merged_df_mapk13[exp_column], s=2, alpha=1, c = 'green')
plt.legend(['PRKCQ', 'MAPK13'])
plt.xlabel('Embedding similarity (cosine)')
plt.ylabel(f'Expression similarity (cosine)')
plt.title(f'Embedding similarity vs Expression similarity (cosine)')

# Same for domain identtity for prkcq
plt.figure()
plt.scatter(merged_df_prkcq['domain_identity'], merged_df_prkcq[exp_column], s=2, alpha=0.5, c = 'red')
plt.scatter(merged_df_mapk13['domain_identity'], merged_df_mapk13[exp_column], s=2, alpha=1, c = 'green')
plt.legend(['PRKCQ', 'MAPK13'])
plt.xlabel('Domain identity')
plt.ylabel(f'Expression similarity (cosine)')
plt.title(f'Domain identity vs Expression similarity (cosine)')




# Same for domain identtity
plt.figure()
plt.scatter(merged_df_sbk['domain_identity'], merged_df_sbk[exp_column], s=2, alpha=0.5, c = 'red')
plt.scatter(merged_df_mapk13['domain_identity'], merged_df_mapk13[exp_column], s=2, alpha=1, c = 'green')
plt.legend(['SBK3/SBK2', 'MAPK13'])
plt.xlabel('Domain identity')
plt.ylabel(f'Expression similarity (cosine)')
plt.title(f'Domain identity vs Expression similarity (cosine)')


# 
plt.figure()
plt.scatter(merged_df_sbk2['emb_sim_cos'], merged_df_sbk2[exp_column], s=2, alpha=0.5, c = 'red')
plt.scatter(merged_df_mapk13['emb_sim_cos'], merged_df_mapk13[exp_column], s=2, alpha=1, c = 'green')
plt.scatter(merged_df_kit['emb_sim_cos'], merged_df_kit[exp_column], s=2, alpha=1, c = 'blue')
plt.legend(['SBK2', 'MAPK13', 'KIT'])
plt.xlabel('Embedding similarity (cosine)')
plt.ylabel(f'Expression similarity (cosine)')
plt.title(f'Embedding similarity vs Expression similarity (cosine)')


# Same for domain identtity
plt.figure()
plt.scatter(merged_df_sbk2['domain_identity'], merged_df_sbk2[exp_column], s=2, alpha=0.5, c = 'red')
plt.scatter(merged_df_mapk13['domain_identity'], merged_df_mapk13[exp_column], s=2, alpha=1, c = 'green')
plt.scatter(merged_df_kit['domain_identity'], merged_df_kit[exp_column], s=2, alpha=1, c = 'blue')
plt.legend(['SBK2', 'MAPK13', 'KIT'])
plt.xlabel('Domain identity')
plt.ylabel(f'Expression similarity (cosine)')
plt.title(f'Domain identity vs Expression similarity (cosine)')




# Select which pairs have emb sim (-0.4, -0.6) and emb sim 0.4, 0.7
high_pairs = notrain_results[(notrain_results['emb_sim_cos'] > 0.6) & (notrain_results['emb_sim_cos'] < 0.75) & (notrain_results['exp_sim_cosine'] > 0.7) & (notrain_results['exp_sim_cosine'] < 0.9)]
print(f"Highl pairs: {high_pairs.shape}")

# See how many times each gene_id appears in the weird pairs
high_pairs_counts = pd.concat([high_pairs['gene_name_human'], high_pairs['gene_name_mouse']]).value_counts()

# 
len(pd.concat([normal_pairs['gene_name_human'], normal_pairs['gene_name_mouse']]).unique())



# Optional: Visualize
plt.figure()
plt.scatter(unknown['emb_sim_cos'], unknown[exp_column], s=2, alpha=0.5, c = 'blue' , label = 'Not labeled')
plt.scatter(pos_neg['emb_sim_cos'], pos_neg[exp_column], s=2, alpha=1, c = 'green', label = 'Negative/Positive')
# plt.scatter(x, y, alpha=0.5, label="Original data", s=2)
plt.plot(x, model.predict(x_reshaped), color='red', label="Linear fit")
plt.xlabel("Embedding similarity")
plt.ylabel("Expression similarity")
plt.legend()
# Add vertical line in -0.45
plt.axvline(x=0.6, color='black', linestyle='--', label="Threshold")
plt.axvline(x=0.75, color='black', linestyle='--', label="Threshold")
# Add horizontal line in 0.5
plt.axhline(y=0.7, color='black', linestyle='--', label="Threshold")
plt.axhline(y=0.9, color='black', linestyle='--', label="Threshold")
plt.title("Linear regression on original data")
plt.savefig(f"./exploratory_analysis/plots/{folder_id}/{label_method}/emb_sim_vs_{exp_column}_{label_method}_linear_regression.png")










# TRIANGLE UNEQUALITY
# %%
import numpy as np
import pandas as pd
import random
from tqdm import tqdm



# %%


def triangle_inequality_gene_pairs(df, sim_col='Exp Sim', num_triplets=10000, seed=42):
    """
    Check triangle inequality on scalar similarities across gene pairs.

    Arguments:
        df: DataFrame with expression similarity values for gene pairs
        sim_col: column with similarity values
        num_triplets: how many triplets to sample
        seed: random seed

    Returns:
        violation_rate: float
        total_violations: int
    """
    np.random.seed(seed)
    sims = df[sim_col].values
    N = len(sims)
    
    violations = 0
    for _ in tqdm(range(num_triplets)):
        i, j, k = np.random.choice(N, size=3, replace=False)
        #print(i, j, k)
        s1, s2, s3 = sims[i], sims[j], sims[k]
        #print(s1, s2, s3)

        d12 = abs(s1 - s2)
        d23 = abs(s2 - s3)
        d13 = abs(s1 - s3)
        #print(d12, d23, d13)

        if d12 + d23 < d13 or d12 + d13 < d23 or d13 + d23 < d12:
            violations += 1

    violation_rate = violations / num_triplets
    return violation_rate, violations

# %%
violation_rate, total = triangle_inequality_gene_pairs(df, sim_col='emb_sim_cos', num_triplets=10000)
print(f"Triangle inequality violated in {violation_rate*100:.2f}% of sampled triplets")


# Add a column with the squared values of emb_sim_cos
df['emb_sim_cos_squared'] = df['emb_sim_cos'] ** 2
# Add a column with the squared values of exp_sim_cosine
df['exp_sim_cosine_squared'] = df['exp_sim_cosine'] ** 2
# Add a column with the squared values of rho
rho_merged_df['Rho_squared'] = rho_merged_df['Rho'] ** 2
rho_merged_df['emb_sim_cos_squared'] = rho_merged_df['emb_sim_cos'] ** 2


# Plot the emb_sim_cos_squared against exp_sim_cosine_squared
plt.figure()
plt.scatter(df_clr_expanded['emb_distance'], df_clr_expanded['exp_sim_cosine'], s=2, alpha=0.5)
plt.xlabel('Embedding similarity squared')
plt.ylabel('Expression similarity')
plt.title('Embedding similarity squared vs Expression similarity')



# %%
violation_rate_rho, total = triangle_inequality_gene_pairs(rho_merged_df, sim_col='Rho', num_triplets=2)
print(f"Triangle inequality violated in {violation_rate*100:.2f}% of sampled triplets")
# %%

rho_merged_df = df.copy()
# Plot emb_sim_cos_squared against rho
plt.figure()
plt.scatter(rho_merged_df['emb_sim_cos'], rho_merged_df['Rho_squared'], s=2, alpha=0.5)
plt.xlabel('Embedding similarity')

# Correlation of emb_sim_cos_squared and rho
corr_rho = np.corrcoef(rho_merged_df['emb_sim_cos_squared'], rho_merged_df['Rho_squared'])
print(f"Correlation of embedding similarity squared and rho: {corr_rho[0][1]}")

# Correlation of emb_sim_cos_squared and rho
corr_rho_emb_sq = np.corrcoef(rho_merged_df['emb_sim_cos_squared'], rho_merged_df['Rho'])
print(f"Correlation of embedding similarity squared and rho: {corr_rho_emb_sq[0][1]}")

# Correlation of emb_sim_cos and rho
corr_rho_emb = np.corrcoef(rho_merged_df['emb_sim_cos'], rho_merged_df['Rho'])
print(f"Correlation of embedding similarity and rho: {corr_rho_emb[0][1]}")



plt.figure()
plt.scatter(rho_merged_df['emb_sim_cos'], rho_merged_df['Rho_distance'], s=2, alpha=0.5)
plt.xlabel('Embedding similarity')
plt.ylabel('Rho')

plt.figure()
plt.scatter(rho_merged_df['emb_sim_cos'], rho_merged_df['Rho'], s=2, alpha=0.5)
plt.xlabel('Embedding similarity')
plt.ylabel('Rho')

# Correlation of emb_sim_cos and rho_distance
corr_rho_distance = np.corrcoef(rho_merged_df['emb_sim_cos'], rho_merged_df['Rho_distance'])
print(f"Correlation of embedding similarity and rho distance: {corr_rho_distance[0][1]}")

# Correlation of emb_sim_cos and rho
corr_rho = np.corrcoef(rho_merged_df['emb_sim_cos'], rho_merged_df['neg_log_Rho'])
print(f"Correlation of embedding similarity and rho squared: {corr_rho[0][1]}")

print(f"Correlation of embedding similarity and rho: {corr_rho_distance[0][1]}")


df_clr_expanded = df.copy()


# For Rho


# Unique of gene_name_human
genes = rho_merged_df['gene_name_human'].unique()
# Choose 3 random genes
random_genes = random.sample(list(genes), 3)
print(random_genes)

gene_1 = random_genes[0]
gene_2 = random_genes[1]
gene_3 = random_genes[2]

# Select the rows with gene_name_human AND gene_name_mouse in the random genes --> 6 pairs
random_genes_df = rho_merged_df[(rho_merged_df['gene_name_human'].isin(random_genes)) & (rho_merged_df['gene_name_mouse'].isin(random_genes))]
print(random_genes_df.shape)

# Keep the even rows
random_genes_df = random_genes_df.iloc[::2, :]

# Check triangle inequality
metric_column = 'Rho'

s1, s2, s3 = random_genes_df[metric_column].iloc[0], random_genes_df[metric_column].iloc[1], random_genes_df[metric_column].iloc[2]
d12 = abs(s1 - s2)
d23 = abs(s2 - s3)
d13 = abs(s1 - s3)
print(f"Gene 1: {gene_1}, Gene 2: {gene_2}, Gene 3: {gene_3}")
print(f"Rho values: {s1}, {s2}, {s3}")
print(f"Distances: {d12}, {d23}, {d13}")
print(f"Triangle inequality violated: {d12 + d23 < d13 or d12 + d13 < d23 or d13 + d23 < d12}")

rho_merged_df_save = rho_merged_df.copy()
rho_merged_df = rho_merged_df.iloc[::2, :]

# Add a column with the squared values of Rho
rho_merged_df['Rho_squared'] = rho_merged_df['Rho'] ** 2
rho_merged_df['neg_log_Rho'] = -np.log(1 - rho_merged_df['Rho'])
rho_merged_df['Rho_distance'] = - (1 - rho_merged_df['Rho'])
rho_merged_df['exp_sim_cosine_distance'] = 1 - rho_merged_df['exp_sim_cosine']

num_triplets = 10000
violations = 0
# Keep even rows

metric_column = 'exp_sim_cosine_distance'

# Precompute unique genes to avoid recalculating in each iteration
unique_genes = rho_merged_df['gene_name_human'].unique()

for _ in tqdm(range(num_triplets)):
    # Randomly sample 3 unique genes
    random_genes = np.random.choice(unique_genes, size=3, replace=False)

    # Filter rows for the selected genes
    random_genes_df = rho_merged_df[
        (rho_merged_df['gene_name_human'].isin(random_genes)) & 
        (rho_merged_df['gene_name_mouse'].isin(random_genes))
    ]

    # Ensure we have exactly 3 rows to proceed
    if len(random_genes_df) != 3:
        print('hey')
        continue

    # Extract the metric values
    s1, s2, s3 = random_genes_df[metric_column].iloc[:3]
    d12 = abs(s1 - s2)
    d23 = abs(s2 - s3)
    d13 = abs(s1 - s3)

    # Check for triangle inequality violation
    if d12 + d23 < d13 or d12 + d13 < d23 or d13 + d23 < d12:
        violations += 1

# Calculate violation rate
violation_rate = violations / num_triplets * 100
print(f"For the metric {metric_column} the rate is {violation_rate}%")


#%%

# Load pairs, keep only the P/N and the columns of interest
pairs_df = pd.read_csv(pairs_file, sep = "\t")
pairs_df = pairs_df[pairs_df["label"] != "U"]
pairs_df = pairs_df[["gene_id_human", "gene_id_mouse", "label"]]


# Map sequences using the ENSEMBL ID
pairs_df["GeneA_Seq"] = pairs_df["gene_id_human"].map(human_sequences)
pairs_df["GeneB_Seq"] = pairs_df["gene_id_mouse"].map(mouse_sequences)


# 