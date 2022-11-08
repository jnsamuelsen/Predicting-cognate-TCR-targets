# Alignment as csv

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

alignment_df = pd.read_csv('/mnt/c/Users/jeppe/Desktop/Speciale2022/Data/Alignments/Full_alignment_cleaned_pdb_files.csv',
                           index_col = 0)

pdb_list = pd.read_csv('/mnt/c/Users/jeppe/Desktop/Speciale2022/Scripts/CleanPDBfiles/pdb_list.csv', index_col = 0)

print(alignment_df.at['1oga', 'peptide_seq'])