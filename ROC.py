# Receiver operating characteristic plots

# Load WT and SW scores
WT_scores = pd.read_csv('peptide_scores_WT.csv', header = 0)
SW_scores = pd.read_csv('peptide_scores_SW.csv', header = 0)
WT_weighted_scores = pd.read_csv('peptide_scores_weighted_WT.csv', header = 0)
SW_weighted_scores = pd.read_csv('peptide_scores_weighted_SW.csv', header = 0)

# Concatenate data
pd_scores = pd.concat([WT_scores, SW_scores], axis = 1)
pd_weighted_scores = pd.concat([WT_weighted_scores, SW_weighted_scores], axis = 1)

# Add the difference between WT and SW scores in the dataframes
pd_scores['Diff'] = pd_scores['WT'] - pd_scores['SW']
pd_weighted_scores['Diff'] = pd_weighted_scores['WT'] - pd_weighted_scores['SW']