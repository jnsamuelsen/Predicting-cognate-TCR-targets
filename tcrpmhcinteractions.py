# This script calculates whether residues in the CDR3 region and the peptide interacts or not
# via atomic distances

import Bio
from Bio import PDB
from Bio.PDB import *
from Bio.PDB.Polypeptide import three_to_one
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# Function that takes as input a pdb file and outputs a matrix of binary residue interactions between CDR3 alpha and beta and the peptide
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

    return interaction_a, interaction_b



find_tcr_p_interactions('/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/CleanPDBfiles/data/cleaned_data/PDBs/1ao7_cleaned.pdb')


### Start of Program
# Merge all the interaction dataframes together and make one heatmap

path_of_dir = '/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/TCRpMHCinteractions/Data'
ext = '_cleaned.pdb'

aa_list = ['ALA', 'CYS', 'ASP', 'GLU',
			   'PHE', 'GLY', 'HIS', 'ILE',
			   'LYS', 'LEU', 'MET', 'ASN',
			   'PRO', 'GLN', 'ARG', 'SER',
			   'THR', 'VAL', 'TRP', 'TYR']

# Preparing data tables for storing of combined amino acid counts
combined_interactions = np.zeros((20,20))
combined_pd_a = pd.DataFrame(combined_interactions, columns = aa_list, index = aa_list).astype(int)
combined_pd_b = pd.DataFrame(combined_interactions, columns = aa_list, index = aa_list).astype(int)

# Takes 1 pdb file at a time and runs them through find_tcr_p_interactions
for files in os.listdir(path_of_dir):
    print(files)
    if files.endswith(ext):
        df_a = find_tcr_p_interactions(files)[0]
        df_b = find_tcr_p_interactions(files)[1]

        # Adds the binary values / counts to the combined interaction dataframes
        for rowIndex, row in df_a.iterrows():
            for columnIndex, value in row.items():
                combined_pd_a.at[rowIndex, columnIndex] += value
        for rowIndex, row in df_b.iterrows():
            for columnIndex, value in row.items():
                combined_pd_b.at[rowIndex, columnIndex] += value

print(combined_pd_a)
print(combined_pd_b)

# Make the combined heatmaps
# Alpha chains
plt.title('Amino acid interaction heatmap for peptide and CDR3 alpha')
sns_a = sns.heatmap(combined_pd_a, cmap = "Blues")
fig = sns_a.get_figure()
plt.xlabel('Peptide residue')
plt.ylabel('CDR3 alpha residue')
plt.gcf().set_size_inches(15, 8)
fig.savefig("sns_alpha.png")

# Beta chains
plt.title('Amino acid interaction heatmap for peptide and CDR3 beta')
sns_b = sns.heatmap(combined_pd_b, cmap = "Blues")
fig = sns_b.get_figure()
plt.xlabel('Peptide residue')
plt.ylabel('CDR3 beta residue')
plt.gcf().set_size_inches(15, 8)
fig.savefig("sns_beta.png")

# Alpha Beta combined
combined_pd_a_b = combined_pd_a.add(combined_pd_b)
plt.title('Amino acid interaction heatmap for peptide and CDR3 alpha/beta combined')
sns_ab = sns.heatmap(combined_pd_a_b, cmap = "Blues")
fig = sns_ab.get_figure()
plt.xlabel('Peptide residue')
plt.ylabel('CDR3 residue')
plt.gcf().set_size_inches(15, 8)
fig.savefig("sns_alpha_beta.png")