# This script takes TCRpMHC structures and calculates the distance between all of them in a heatmap

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


# Initialize variables
combined_interactions = np.zeros((20,20))

aa_list = ['ALA', 'CYS', 'ASP', 'GLU',
			   'PHE', 'GLY', 'HIS', 'ILE',
			   'LYS', 'LEU', 'MET', 'ASN',
			   'PRO', 'GLN', 'ARG', 'SER',
			   'THR', 'VAL', 'TRP', 'TYR']

def find_tcr_p_interactions(pdbfile):

    # Parse pdb file
    parser = PDBParser()
    structure = parser.get_structure('tcrpmhc', pdbfile)[0]

    # Get chains (alpha, beta and peptide)
    chain_a = structure['A']
    chain_b = structure['B']
    chain_p = structure['P']

    # Picking out only the CDR3 regions chain alpha and beta according to LYRA range for CDR3
    CDR3_a_resi = []
    CDR3_b_resi = []

    # Creating lists for axis labels in interaction matrix
    CDR3_a_header = []
    CDR3_b_header = []
    peptide_header = []

    pd_a = pd.DataFrame(combined_interactions, columns=aa_list, index=aa_list).astype(float)
    pd_b = pd.DataFrame(combined_interactions, columns=aa_list, index=aa_list).astype(float)
    pd_ab = pd.DataFrame(combined_interactions, columns=aa_list, index=aa_list).astype(float)

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

    # Write into pd dataframes
    interaction_a = pd.DataFrame(interaction_a, columns = peptide_header, index = CDR3_a_header).astype(float)
    interaction_b = pd.DataFrame(interaction_b, columns = peptide_header, index = CDR3_b_header).astype(float)
    interaction_ab = pd.DataFrame(combined_interactions, columns = aa_list, index = aa_list).astype(float)

    # Adds the binary values / counts to the combined interaction dataframes (CDR3 alpha and beta)
    for rowIndex, row in interaction_a.iterrows():
        for columnIndex, value in row.items():
            if value == 1 and interaction_ab.at[rowIndex, columnIndex] == 0:
                interaction_ab.at[rowIndex, columnIndex] += value
    for rowIndex, row in interaction_b.iterrows():
        for columnIndex, value in row.items():
            if value == 1 and interaction_ab.at[rowIndex, columnIndex] == 0:
                interaction_ab.at[rowIndex, columnIndex] += value

    return interaction_ab

# Set directory path
path_of_dir = '/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/TCRpMHCinteractions/Data'
ext = '.pdb'

# Create lists
file_list = []
matrix_list = []

# Loop over all pdb files and find the binary interaction matrices for each of them
for file in os.listdir(path_of_dir):

    pdb_file = os.path.join(path_of_dir, file)

    # Append file name and interaction matrix for pdb
    file_list.append(file)
    matrix_list.append(find_tcr_p_interactions(pdb_file))


# Initialize distance matrix and matrix counter
distance_matrix = np.zeros((len(file_list),len(file_list)))
matrix_count = 0

for matrix in matrix_list:
    for i in range(len(file_list)):
        # Find the euclidean distance between the structures
        distance_matrix[matrix_count, i] = np.linalg.norm(matrix_list[i] - matrix)

    matrix_count += 1

# Make into df
distance_df = pd.DataFrame(distance_matrix, columns=file_list, index=file_list).astype(float)

print(distance_df)

# Alpha Beta combined - heatmap
plt.title('Distance between all TCRpMHC structures (AA interactions)', fontsize = 20)
sns_distance = sns.heatmap(distance_df, cmap = "vlag")
fig = sns_distance.get_figure()
plt.xlabel('PDB structure', fontsize = 15)
plt.ylabel('PDB structure', fontsize = 15)
plt.gcf().set_size_inches(21, 16)
fig.savefig("sns_structure_distance_heatmap.png")
