# This code scores the structures of WT TCRpMHC structures and structures where the peptide is swapped
# out to another peptide which binds to the MHC and has a similarity of 6/9 AAs

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

# Empty matrix
combined_interactions = np.zeros((20, 20))


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
    interaction_ab = pd.DataFrame(combined_interactions, columns = aa_list, index = aa_list).astype(float)

    # Adds the binary values / counts to the combined interaction dataframes
    for rowIndex, row in interaction_a.iterrows():
        for columnIndex, value in row.items():
            if value == 1 and interaction_ab.at[rowIndex, columnIndex] == 0:
                interaction_ab.at[rowIndex, columnIndex] += value
    for rowIndex, row in interaction_b.iterrows():
        for columnIndex, value in row.items():
            if value == 1 and interaction_ab.at[rowIndex, columnIndex] == 0:
                interaction_ab.at[rowIndex, columnIndex] += value

    return interaction_ab

peptide_scores_WT = []
peptide_scores_SW = []

file_list_WT = []
file_list_SW = []

path_of_dir_WT = '/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/Scoringstructures/PDB/WT'
path_of_dir_SW = '/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/Scoringstructures/PDB/SW'
ext = '.pdb'

# Loop over WT pdb files
# Score the leave WT TCR-p structure
for file in os.listdir(path_of_dir_WT):

    WTfile = os.path.join(path_of_dir_WT, file)

    peptide_score = 0

    WT_ab = find_tcr_p_interactions(WTfile)

    # Score the peptide
    for rowIndex, row in energy_matrix.iterrows():
        for columnIndex, value in row.items():
            if WT_ab.at[rowIndex, columnIndex] == 1:
                peptide_score += value

    # Append the score to peptide score list
    peptide_scores_WT.append(peptide_score)
    file_list_WT.append(file)

    print(peptide_scores_WT)

# Make plot for the peptide WT scores
fig = plt.figure(figsize = (21, 16))
plt.bar(file_list_WT, peptide_scores_WT, color = 'maroon')
plt.xlabel('PDB WT structure', fontsize = 15)
plt.ylabel('Score', fontsize = 15)
plt.xticks(rotation=90, fontsize = 8)
plt.title('Wildtype peptide scorings', fontsize = 20)
fig.savefig("WT_scorings.png")

# Loop over SW pdb files
# Score the leave SW TCR-p structure
for file in os.listdir(path_of_dir_SW):

    SWfile = os.path.join(path_of_dir_SW, file)

    peptide_score = 0

    SW_ab = find_tcr_p_interactions(SWfile)

    # Score the peptide
    for rowIndex, row in energy_matrix.iterrows():
        for columnIndex, value in row.items():
            if SW_ab.at[rowIndex, columnIndex] == 1:
                peptide_score += value

    # Append the score to peptide score list
    peptide_scores_SW.append(peptide_score)
    file_list_SW.append(file)

    print(peptide_scores_SW)

print(file_list_WT)

print(file_list_SW)

# Make plot for the peptide SW scores
fig = plt.figure(figsize=(21, 16))
plt.bar(file_list_SW, peptide_scores_SW, color='maroon')
plt.xlabel('PDB SW structure', fontsize = 15)
plt.ylabel('Score', fontsize = 15)
plt.xticks(rotation=90, fontsize = 8)
plt.title('Swapped peptide scorings', fontsize = 20)
fig.savefig("SW_scorings.png")

