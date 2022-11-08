# Pseudocode
# For all CDR3 structures count the observed amount of all amino acids
# Repeat that for the peptides of all the structures
# For all the positions in the combined AA matrix:
# Add pseudocount
# calculate pos_ab = ln (pos_ab (obs) / (pos_a (exp) * pos_b (exp)))

# Repeat the matrix calculation but leave structure out and then score all the peptides one at a time

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

aa_list = ['ALA', 'CYS', 'ASP', 'GLU',
			   'PHE', 'GLY', 'HIS', 'ILE',
			   'LYS', 'LEU', 'MET', 'ASN',
			   'PRO', 'GLN', 'ARG', 'SER',
			   'THR', 'VAL', 'TRP', 'TYR']

combined_interactions = np.zeros((20,20))


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

def count_CDR3_peptide_residues(pdbfile):
    # Parse pdb file
    parser = PDBParser()
    structure = parser.get_structure('tcrpmhc', pdbfile)[0]

    # Dictionaries for counting residues in combined peptides and CDR3 alpha/beta residues
    peptide_dict = {'ALA': 0, 'CYS': 0, 'ASP': 0, 'GLU': 0,
                   'PHE': 0, 'GLY': 0, 'HIS': 0, 'ILE': 0,
                   'LYS': 0, 'LEU': 0, 'MET': 0, 'ASN': 0,
                   'PRO': 0, 'GLN': 0, 'ARG': 0, 'SER': 0,
                   'THR': 0, 'VAL': 0, 'TRP': 0, 'TYR': 0}
    CDR3_dict = {'ALA': 0, 'CYS': 0, 'ASP': 0, 'GLU': 0,
                   'PHE': 0, 'GLY': 0, 'HIS': 0, 'ILE': 0,
                   'LYS': 0, 'LEU': 0, 'MET': 0, 'ASN': 0,
                   'PRO': 0, 'GLN': 0, 'ARG': 0, 'SER': 0,
                   'THR': 0, 'VAL': 0, 'TRP': 0, 'TYR': 0}

    # Get chains (alpha, beta and peptide)
    chain_a = structure['A']
    chain_b = structure['B']
    chain_p = structure['P']

    # Counting residues in the CDR3 residues of alpha and beta chain
    for residue_a in chain_a:
        if residue_a.id[1] in range(107, 129):
            CDR3_dict[residue_a.resname] += 1

    for residue_b in chain_b:
        if residue_b.id[1] in range(104, 129):
            CDR3_dict[residue_a.resname] += 1

    for residue_p in chain_p:
        peptide_dict[residue_p.resname] += 1

    # To get frequencies divide by total number of residues
    # For P_obs divide by total numbers of contacts

    return peptide_dict, CDR3_dict


# Combine all the files together into two count matrices (peptide and CDR3)
# Make the final energy matrix and run leave one out

### Start of Program
# Merge all the interaction dataframes together and make one heatmap

path_of_dir = '/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/TCRpMHCinteractions/Data'
ext = '_cleaned.pdb'

# Preparing data tables for storing of combined amino acid counts
combined_peptide_dict = {'ALA': 0, 'CYS': 0, 'ASP': 0, 'GLU': 0,
                   'PHE': 0, 'GLY': 0, 'HIS': 0, 'ILE': 0,
                   'LYS': 0, 'LEU': 0, 'MET': 0, 'ASN': 0,
                   'PRO': 0, 'GLN': 0, 'ARG': 0, 'SER': 0,
                   'THR': 0, 'VAL': 0, 'TRP': 0, 'TYR': 0}

combined_CDR3_dict = {'ALA': 0, 'CYS': 0, 'ASP': 0, 'GLU': 0,
                   'PHE': 0, 'GLY': 0, 'HIS': 0, 'ILE': 0,
                   'LYS': 0, 'LEU': 0, 'MET': 0, 'ASN': 0,
                   'PRO': 0, 'GLN': 0, 'ARG': 0, 'SER': 0,
                   'THR': 0, 'VAL': 0, 'TRP': 0, 'TYR': 0}

combined_pd_a_b = pd.DataFrame(combined_interactions, columns = aa_list, index = aa_list).astype(float)

# Initialize count parameters
total_resp_count = 0
total_resCDR3_count = 0
total_contact_count = 0

# List for leave one out
file_list = []

# Takes 1 pdb file at a time and runs them through find_tcr_p_interactions
for files in os.listdir(path_of_dir):
    #print(files)
    if files.endswith(ext):
        # Add residue counts to dicts
        for key, val in count_CDR3_peptide_residues(files)[0].items():
            combined_peptide_dict[key] += val
            total_resp_count += val
        for key, val in count_CDR3_peptide_residues(files)[1].items():
            combined_CDR3_dict[key] += val
            total_resCDR3_count += val

        # Combine the contact counts into combined dataframes
        df_ab = find_tcr_p_interactions(files)
        #df_b = find_tcr_p_interactions(files)[1]

        #print(df_a)
        #print(df_b)

        # Adds the binary values / counts to the combined interaction dataframes
        for rowIndex, row in df_ab.iterrows():
            for columnIndex, value in row.items():
                combined_pd_a_b.at[rowIndex, columnIndex] += value
                total_contact_count += value


# Comined interaction matrix
# combined_pd_a_b = combined_pd_a.add(combined_pd_b)

print(combined_pd_a_b)

# Pseudocount
pseudocount = 1
total_contact_count = total_contact_count + pseudocount * 400

print(total_contact_count)

# Add pseudocount and normalize interaction matrix
for rowIndex, row in combined_pd_a_b.iterrows():
    for columnIndex, value in row.items():
        combined_pd_a_b.at[rowIndex, columnIndex] += pseudocount
        combined_pd_a_b.at[rowIndex, columnIndex] = combined_pd_a_b.at[rowIndex, columnIndex] / total_contact_count

print(combined_pd_a_b)

print(combined_pd_a_b.sum())

# Making empty matrix for the final calculation of the energy matrix
final_energy_pd = pd.DataFrame(combined_interactions, columns = aa_list, index = aa_list).astype(float)

# Making the energy matrix
for rowIndex, row in final_energy_pd.iterrows():
    for columnIndex, value in row.items():
        p_obs = (combined_pd_a_b.at[rowIndex, columnIndex])
        p_exp = (((combined_pd_a_b[rowIndex].sum())) * (
                (combined_pd_a_b[columnIndex].sum())))
        final_energy_pd.at[rowIndex, columnIndex] = np.log(p_obs / p_exp)

print(final_energy_pd)

# Save csv
final_energy_pd.to_csv('energy_matrix.csv')

# Alpha Beta combined - heatmap
plt.title('Energy heatmap for interactions between peptide and CDR3 alpha/beta residues')
sns_ab = sns.heatmap(final_energy_pd, cmap = "vlag")
fig = sns_ab.get_figure()
plt.xlabel('Peptide residue')
plt.ylabel('CDR3 residue')
plt.gcf().set_size_inches(15, 8)
fig.savefig("sns_energy_heatmap.png")
