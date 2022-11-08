# This script takes an energy matrix and calculates the distance between the amino acids

# Find distance between rows and columns in the asymmetric matrix

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
