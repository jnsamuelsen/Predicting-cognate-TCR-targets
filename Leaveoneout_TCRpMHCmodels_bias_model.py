# For all CDR3 structures count the observed amount of all amino acids
# Repeat that for the peptides of all the structures
# For all the positions in the combined AA matrix:
# Add pseudocount
# calculate pos_ab = ln (pos_ab (obs) / (pos_a (exp) * pos_b (exp)))

# Repeat the matrix calculation but leave structure out and then score all the peptides one at a time
import sys

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
import distance

aa_list = ['ALA', 'CYS', 'ASP', 'GLU',
			   'PHE', 'GLY', 'HIS', 'ILE',
			   'LYS', 'LEU', 'MET', 'ASN',
			   'PRO', 'GLN', 'ARG', 'SER',
			   'THR', 'VAL', 'TRP', 'TYR']

convert3_1  = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
			   'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
			   'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
			   'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
			   'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

# Empty matrix
combined_interactions = np.zeros((20, 20))

# Import alignment
#alignment_df = pd.read_csv('/mnt/c/Users/jeppe/Desktop/Speciale2022/Data/Alignments/Full_alignment_cleaned_pdb_files.csv',
#                           index_col = 0)

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

    # Getting the peptide sequence as a string
    peptide_seq = ''.join([convert3_1[x] for x in peptide_header])

    return interaction_ab, peptide_seq


# Combine all the files together into two count matrices (peptide and CDR3)
# Make the final energy matrix and run leave one out

### Start of Program
# Merge all the interaction dataframes together and make one heatmap

path_of_dir = '/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/Scoring_using_models/Data'
path_of_dir_scoring = '/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/Scoringstructures/PDB/SW'
ext = '_cleaned.pdb'

# List for leave one out
file_list = []

# Peptide scores
peptide_scores = []

# Pseudocount
pseudocount = 1

# Takes 1 pdb file at a time and runs them through find_tcr_p_interactions
for leaveoutfile in os.listdir(path_of_dir_scoring):
    print(leaveoutfile)

    WTfile = os.path.join(path_of_dir_scoring, leaveoutfile)

    # Initialize matrices
    combined_pd_a_b = pd.DataFrame(combined_interactions, columns=aa_list, index=aa_list).astype(float)

    peptide_score = 0

    # Initialize count parameters
    total_contact_count = 0

    # Leave out TCR-p structure and peptide sequence
    leave_out_ab = find_tcr_p_interactions(WTfile)[0]
    leave_out_pep_seq = find_tcr_p_interactions(WTfile)[1]

    for files in os.listdir(path_of_dir):

        # If the peptide is not a 9-mer
        # or if the peptide is a 9-mer and the peptide leaveoutfile is max 6/9 AA identical to the peptide from the file
        ## If the tcr A / B from the leaveoutfile is max 95 % identical to the tcr A / B from the file
        ## Proceed

        pep_seq = find_tcr_p_interactions(files)[1]

        if len(leave_out_pep_seq) != len(pep_seq):
            if files.endswith(ext):

                # Combine the contact counts into combined dataframes
                df_ab = find_tcr_p_interactions(files)[0]

                # Adds the binary values / counts to the combined interaction dataframes
                for rowIndex, row in df_ab.iterrows():
                    for columnIndex, value in row.items():
                        combined_pd_a_b.at[rowIndex, columnIndex] += value
                        total_contact_count += value

        elif distance.levenshtein(leave_out_pep_seq, pep_seq) >= 3:
            if files.endswith(ext):

                # Combine the contact counts into combined dataframes
                df_ab = find_tcr_p_interactions(files)[0]

                # Adds the binary values / counts to the combined interaction dataframes
                for rowIndex, row in df_ab.iterrows():
                    for columnIndex, value in row.items():
                        combined_pd_a_b.at[rowIndex, columnIndex] += value
                        total_contact_count += value

    total_contact_count = total_contact_count + pseudocount * 400

    # Add pseudocount to interaction matrix
    for rowIndex, row in combined_pd_a_b.iterrows():
        for columnIndex, value in row.items():
            combined_pd_a_b.at[rowIndex, columnIndex] += pseudocount
            combined_pd_a_b.at[rowIndex, columnIndex] = combined_pd_a_b.at[rowIndex, columnIndex] / total_contact_count

    # Making empty matrix for the final calculation of the energy matrix
    final_energy_pd = pd.DataFrame(combined_interactions, columns = aa_list, index = aa_list).astype(float)

    # Making the energy matrix
    for rowIndex, row in final_energy_pd.iterrows():
        for columnIndex, value in row.items():
            p_obs = (combined_pd_a_b.at[rowIndex, columnIndex])
            p_exp = (((combined_pd_a_b[rowIndex].sum())) * (
                    (combined_pd_a_b[columnIndex].sum())))
            final_energy_pd.at[rowIndex, columnIndex] = np.log(p_obs / p_exp)

    # Score the peptide
    for rowIndex, row in final_energy_pd.iterrows():
        for columnIndex, value in row.items():
            if leave_out_ab.at[rowIndex, columnIndex] == 1:
                peptide_score += value

    # Append leave out file to file list and append the score to peptide score list
    file_list.append(leaveoutfile)
    peptide_scores.append(peptide_score)

    print(file_list)
    print(peptide_scores)

# Make plot for the peptide scores
fig = plt.figure(figsize = (21, 16))
plt.bar(file_list, peptide_scores, color = 'crimson')
plt.xlabel('PDB SW structure', fontsize = 15)
plt.ylabel('Score', fontsize = 15)
plt.xticks(rotation=90, fontsize = 8)
plt.title('Leave one out swapped peptide scorings (using models created from TCRpMHCmodels)', fontsize = 20)
fig.savefig("SW_scorings_minusbias_usingmodels.png")

peptide_scores_df = pd.DataFrame(peptide_scores, index = file_list).astype(float)

peptide_scores_df.to_csv('peptide_scores_SW_minusbias_usingmodels.csv')