"""
Scripts that creates the database of classifiers
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

# Input:
#  - one multiple sequence alignment (MSA) per marker gene. The MSA is obtained
#    from the function htc align, like:
#       >gene1\t0001010001010000100101000...
#       >gene2\t0000110001010100100101001...
#  - a taxonomy file that describes the taxonomy of the genes:
#       gene1\tBacteria\tFirmicutes\t...
#
# Output:
#  - a database file (hdf5) that can be used by htc classify

from sklearn.linear_model import LogisticRegression
import numpy as np
import sys
import random
import time
import pandas as pd

#===============================================================================
#                                 FUNCTIONS
#===============================================================================


# function to load an alignment produced by the "align" option =================
# Input:
#  - a file created by "align"
# Output:
#  - a panda object
# as a note, numpy.loadtxt is way slower than pandas read.csv
# It works also on .gz files
def load_alignment_from_file(file_name):
    alignment = pd.read_csv(file_name,delimiter='\t',index_col = 0, header=None)
    alignment = alignment.astype('bool') # apparently you cannot load directly
                                         # bool if the rownames are not bool
    return alignment

# function to load the taxonomy ================================================
