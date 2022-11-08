# # This code scores each position in WT TCRpMHC structures and structures where the peptide is swapped
# # out to another peptide which binds to the MHC and has a similarity of 6/9 AAs

import Bio
from Bio import PDB
from Bio.PDB import *
from Bio.PDB.Polypeptide import three_to_one
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import itertools

# Load energy matrix
energy_matrix = pd.read_csv('energy_matrix.csv', index_col = 0)

aa_list = ['ALA', 'CYS', 'ASP', 'GLU',
			   'PHE', 'GLY', 'HIS', 'ILE',
			   'LYS', 'LEU', 'MET', 'ASN',
			   'PRO', 'GLN', 'ARG', 'SER',
			   'THR', 'VAL', 'TRP', 'TYR']

# Empty matrix (position specific)
combined_interactions = np.zeros((20, 9))

def find_tcr_p_interactions(pdbfile):

    # Parse pdb file
    parser = PDBParser()
    structure = parser.get_structure('tcrpmhc', pdbfile)[0]

    # Get chains (alpha, beta and peptide)
    chain_a = structure['A']
    chain_b = structure['B']
    chain_p = structure['P']

    # Picking out only the CDR3 regions chain alpha and beta according to LYRA ranges for CDR3
    CDR3_a_resi = []
    CDR3_b_resi = []

    # Creating lists for axis labels in interaction matrix
    CDR3_a_header = []
    CDR3_b_header = []
    peptide_header = []

    # Fetching the CDR3 residues of alpha and beta chain
    for residue_a in chain_a:
        if residue_a.id[1] in range(107,129):
            CDR3_a_resi.append(residue_a)
            CDR3_a_header.append(residue_a.resname)

    for residue_b in chain_b:
        if residue_b.id[1] in range(104,129):
            CDR3_b_resi.append(residue_b)
            CDR3_b_header.append(residue_b.resname)

    # Make empty interaction matrices
    interaction_a = np.zeros((len(CDR3_a_resi),len(chain_p)))
    interaction_b = np.zeros((len(CDR3_b_resi),len(chain_p)))

    # For each of the atoms in the peptide see which atom interactions there are to chain alpha and beta
    for residue in chain_p:
        peptide_header.append(residue.resname)
        for atom in residue:

            if atom.name != 'CA':
                CDR3_a_count = 0
                CDR3_b_count = 0

                for CDR3_a_residue in CDR3_a_resi:
                    CDR3_a_atoms = Bio.PDB.Selection.unfold_entities(CDR3_a_residue, 'A')
                    ns_a = Bio.PDB.NeighborSearch(CDR3_a_atoms)

                    if ns_a.search(atom.coord, 4, 'R'):
                        interaction_a[CDR3_a_count][residue.id[1]-1] = 1

                    CDR3_a_count += 1

                for CDR3_b_residue in CDR3_b_resi:
                    CDR3_b_atoms = Bio.PDB.Selection.unfold_entities(CDR3_b_residue, 'A')
                    ns_b = Bio.PDB.NeighborSearch(CDR3_b_atoms)

                    if ns_b.search(atom.coord, 4, 'R'):
                        interaction_b[CDR3_b_count][residue.id[1]-1] = 1

                    CDR3_b_count += 1

    interaction_a = pd.DataFrame(interaction_a, columns = peptide_header, index = CDR3_a_header).astype(int)
    interaction_b = pd.DataFrame(interaction_b, columns = peptide_header, index = CDR3_b_header).astype(int)
    #interaction_ab = pd.DataFrame(combined_interactions, columns = peptide_header, index = aa_list).astype(float)

    #print(interaction_a)
    #print(interaction_b)
    #print(interaction_ab)

    '''# Adds the binary values / counts to the combined interaction dataframes
    for rowIndex, row in interaction_a.iterrows():
        for columnIndex, value in row.items():
            print(value)
            if value == 1 and interaction_ab.at[rowIndex, columnIndex] == 0:
                interaction_ab.at[rowIndex, columnIndex] += value
    for rowIndex, row in interaction_b.iterrows():
        for columnIndex, value in row.items():
            if value == 1 and interaction_ab.at[rowIndex, columnIndex] == 0:
                interaction_ab.at[rowIndex, columnIndex] += value'''

    return interaction_a, interaction_b

# File paths
path_of_dir_WT = '/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/Scoringstructures/PDB/WT'
path_of_dir_SW = '/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/Scoringstructures/PDB/SW'
ext = '.pdb'

# Make list of lists for position specific scores for WT and SW
pos_scores_WT = []
pos_scores_SW = []

for i in range(9):
    pos_scores_WT.append([])
    pos_scores_SW.append([])

# Loop over WT pdb files
# Score the leave WT TCR-p structure
for file in os.listdir(path_of_dir_WT):

    WTfile = os.path.join(path_of_dir_WT, file)

    WT_a = find_tcr_p_interactions(WTfile)[0]
    WT_b = find_tcr_p_interactions(WTfile)[1]

    # Storage for scores
    temp_pos_list_WT_a = [0,0,0,0,0,0,0,0,0]
    temp_pos_list_WT_b = [0,0,0,0,0,0,0,0,0]

    # Score the peptide positions
    for rowIndex, row in WT_a.iterrows():
        pos_count_WT = 0
        for columnIndex, value in row.items():
            if value == 1:
                temp_pos_list_WT_a[pos_count_WT] += energy_matrix.at[rowIndex, columnIndex]
            pos_count_WT += 1

    for rowIndex, row in WT_b.iterrows():
        pos_count_WT = 0
        for columnIndex, value in row.items():
            if value == 1:
                temp_pos_list_WT_b[pos_count_WT] += energy_matrix.at[rowIndex, columnIndex]
            pos_count_WT += 1

    # Add the scores to pos_scores_WT list of lists
    for i in range(len(temp_pos_list_WT_a)):
        pos_scores_WT[i].append(temp_pos_list_WT_a[i] + temp_pos_list_WT_b[i])

# Loop over SW pdb files
# Score the leave SW TCR-p structure
for file in os.listdir(path_of_dir_SW):

    SWfile = os.path.join(path_of_dir_SW, file)

    SW_a = find_tcr_p_interactions(SWfile)[0]
    SW_b = find_tcr_p_interactions(SWfile)[1]

    # Storage for scores
    temp_pos_list_SW_a = [0,0,0,0,0,0,0,0,0]
    temp_pos_list_SW_b = [0,0,0,0,0,0,0,0,0]

    # Score the peptide positions
    for rowIndex, row in SW_a.iterrows():
        pos_count_SW = 0
        for columnIndex, value in row.items():
            if value == 1:
                temp_pos_list_SW_a[pos_count_SW] += energy_matrix.at[rowIndex, columnIndex]
            pos_count_SW += 1

    for rowIndex, row in SW_b.iterrows():
        pos_count_SW = 0
        for columnIndex, value in row.items():
            if value == 1:
                temp_pos_list_SW_b[pos_count_SW] += energy_matrix.at[rowIndex, columnIndex]
            pos_count_SW += 1

    # Add the scores to pos_scores_SW list of lists
    for i in range(len(temp_pos_list_SW_a)):
        pos_scores_SW[i].append(temp_pos_list_SW_a[i] + temp_pos_list_SW_b[i])


pos_scores_WT_df = pd.DataFrame(pos_scores_WT, index = ['1', '2', '3', '4', '5', '6', '7', '8', '9']).transpose()
pos_scores_SW_df = pd.DataFrame(pos_scores_SW, index = ['1', '2', '3', '4', '5', '6', '7', '8', '9']).transpose()

print(pos_scores_WT_df)
print(pos_scores_SW_df)

plot1 = sns.boxplot(data=pos_scores_WT_df[['1', '2', '3', '4', '5', '6', '7', '8', '9']], orient = 'v')
fig1 = plot1.get_figure()
plt.ylabel('Score')
plt.xlabel('Peptide position')
plt.title('Peptide position scores for non-weighted wildtype TCRpMHCmodels')
plt.gcf().set_size_inches(15, 8)
fig1.savefig("posWT_boxplot_TCRpMHCmodels_nonweighted.png")

'''plot2 = sns.boxplot(data=pos_scores_SW_df[['1', '2', '3', '4', '5', '6', '7', '8', '9']], orient = 'v')
fig2 = plot2.get_figure()
plt.ylabel('Score')
plt.xlabel('Peptide position')
plt.title('Peptide position scores for weighted swapped TCRpMHCmodels')
plt.gcf().set_size_inches(15, 8)
fig2.savefig("posSW_boxplot_TCRpMHCmodels_weighted.png")'''