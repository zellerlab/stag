"""
Scripts that find the taxonomy of an aligned sequence
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import numpy as np
import sys
import time
import pandas as pd
import os
import h5py

#===============================================================================
#                            LOAD THE HDF5 DATABASE
#===============================================================================
def load_DB(hdf5_DB_path):
    f = h5py.File(hdf5_DB_path, 'r')
    f.close()





#===============================================================================
#                                      MAIN
#===============================================================================

def find_tax(database, fasta_input, verbose, output):
    # load the database
    load_DB(database)
