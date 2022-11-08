# Boxplots for the TCRpMHCmodel scores (leave one out scores)

from Bio.PDB import *
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats


# Load WT and SW scores
WT_scores = pd.read_csv('peptide_scores_WT_minusbias_usingmodels.csv', header = 0)
SW_scores = pd.read_csv('peptide_scores_SW_minusbias_usingmodels.csv', header = 0)
#WT_weighted_scores = pd.read_csv('peptide_scores_weighted_WT_minusbias.csv', header = 0)
#SW_weighted_scores = pd.read_csv('peptide_scores_weighted_SW_minusbias.csv', header = 0)

# Concatenate data
pd_scores = pd.concat([WT_scores, SW_scores], axis = 1)
#pd_weighted_scores = pd.concat([WT_weighted_scores, SW_weighted_scores], axis = 1)

# Add the difference between WT and SW scores in the dataframes
pd_scores['Diff'] = pd_scores['WT'] - pd_scores['SW']
#pd_weighted_scores['Diff'] = pd_weighted_scores['WT'] - pd_weighted_scores['SW']

# Seq ID percentages - highest match
seq_ID_WT = pd.read_csv('/mnt/c/Users/jeppe/Desktop/Speciale2022/Data/TCRpMHCmodels_new-main/Template_complex_seqidentity_highesthit_WT.csv',
                        header = 0)
seq_ID_SW = pd.read_csv('/mnt/c/Users/jeppe/Desktop/Speciale2022/Data/TCRpMHCmodels_new-main/Template_complex_seqidentity_highesthit_SW.csv',
                        header = 0)
#print(seq_ID_WT)
#print(WT_scores)

print(stats.mannwhitneyu(x=pd_scores['WT'], y=pd_scores['SW'], alternative = 'greater'))

# Make boxplots
# Non-weighted models
'''
plot1 = sns.boxplot(data=pd_scores[['WT','SW']], orient = 'v')
fig1 = plot1.get_figure()
plt.ylabel('Score')
plt.xlabel('Model')
plt.title('Scores for TCRpMHCmodels (WT = wildtype, SW = swapped) - minus bias, using TCRpMHCmodels for force field\np-value=0.755')
plt.gcf().set_size_inches(15, 8)
fig1.savefig("boxplot_scorings_minusbias_usingmodels.png")
'''
'''
# Weighted models
plot2 = sns.boxplot(data=pd_weighted_scores[['WT','SW']], orient = 'v')
fig2 = plot2.get_figure()
plt.ylabel('Score')
plt.xlabel('Model')
plt.title('Scores for weighted TCRpMHCmodels (WT = wildtype, SW = swapped) - minus bias')
plt.gcf().set_size_inches(15, 8)
fig2.savefig("boxplot_TCRpMHCmodels_weighted_minusbias.png")
'''
'''
# Make plot for the peptide scores
fig3 = plt.figure(figsize = (21, 16))
plt.bar(WT_scores['PDB_structure'], pd_scores['Diff'], color = 'Crimson')
plt.xlabel('PDB structure', fontsize = 15)
plt.ylabel('Score difference', fontsize = 15)
plt.xticks(rotation=90, fontsize = 8)
plt.title('Difference (WT - SW) in leave one out peptide scorings - minus bias, using TCRpMHCmodels for force field', fontsize = 20)
fig3.savefig("nonweighted_diff_scorings_minusbias_usingmodels.png")
'''
'''
# Make plot for the peptide scores
fig4 = plt.figure(figsize = (21, 16))
plt.bar(WT_scores['PDB_structure'], pd_weighted_scores['Diff'], color = 'mediumseagreen')
plt.xlabel('PDB structure', fontsize = 15)
plt.ylabel('Score difference', fontsize = 15)
plt.xticks(rotation=90, fontsize = 8)
plt.title('Difference (WT - SW) in weighted leave one out peptide scorings - minus bias', fontsize = 20)
fig4.savefig("weighted_diff_scorings_minusbias.png")
'''
'''
plot1 = sns.boxplot(data=pd_scores[['Diff']], orient = 'v')
fig1 = plot1.get_figure()
plt.ylabel('Score diff')
#plt.xlabel('')
plt.title('Boxplot - difference (WT - SW) in non-weighted leave one out peptide scorings')
plt.gcf().set_size_inches(15, 8)
fig1.savefig("boxplotdiff_TCRpMHCmodels_nonweighted.png")'''

'''plot1 = sns.boxplot(data=pd_weighted_scores[['Diff']], orient = 'v')
fig1 = plot1.get_figure()
plt.ylabel('Score diff')
#plt.xlabel('')
plt.title('Boxplot - difference (WT - SW) in weighted leave one out peptide scorings')
plt.gcf().set_size_inches(15, 8)
fig1.savefig("boxplotdiff_TCRpMHCmodels_weighted.png")'''

'''
# Make plot for the WT scores pr seq id
fig5 = plt.figure(figsize = (21, 16))
plt.scatter(seq_ID_WT['total_identity'], WT_scores['WT'])
plt.xlabel('PDB model seq ID %', fontsize = 15)
plt.ylabel('WT score', fontsize = 15)
plt.title('WT scorings pr seq ID %', fontsize = 20)
fig5.savefig("WTscoringsprseqID_minusbias.png")
'''

'''
# Make plot for the SW scores pr seq id
fig6 = plt.figure(figsize = (21, 16))
plt.scatter(seq_ID_SW['total_identity'], SW_scores['SW'])
plt.xlabel('PDB model seq ID %', fontsize = 15)
plt.ylabel('SW score', fontsize = 15)
plt.title('SW scorings pr seq ID %', fontsize = 20)
fig6.savefig("SWscoringsprseqID_minusbias.png")
'''
'''
# Make plot for the peptide delta scores pr delta seq id
fig7 = plt.figure(figsize = (21, 16))
plt.scatter(seq_ID_WT['total_identity'] - seq_ID_SW['total_identity'], pd_scores['Diff'])
plt.xlabel('PDB model seq ID % diff', fontsize = 15)
plt.ylabel('Delta score', fontsize = 15)
plt.title('Delta scorings pr delta seq ID %', fontsize = 20)
fig7.savefig("deltascoringsprseqID_minusbias.png")
'''