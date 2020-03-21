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
import tempfile

#===============================================================================
#                            LOAD THE HDF5 DATABASE
#===============================================================================
def load_DB(hdf5_DB_path):
    f = h5py.File(hdf5_DB_path, 'r')

    # first, we save a temporary file with the hmm file ------------------------
    hmm_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(hmm_file.name, 0o644)
    hmm_file.write(f['hmm_file'][0])
    hmm_file.flush()
    os.fsync(hmm_file.fileno())
    hmm_file.close()

    # second if we need to use cm_align ----------------------------------------
    use_cmalign = f['use_cmalign'] # bool

    # third: taxonomy ----------------------------------------------------------
    taxonomy = dict()
    for i in f['taxonomy']:
        taxonomy[i] = list(f['taxonomy/'+i])

    # fourth: tax_function -----------------------------------------------------
    tax_function = dict()
    for r in f['tax_function']:
        tax_function[int(r)] = dict()
        for e in f['tax_function/'+r]:
            tax_function[int(r)][int(e)] = list(f["tax_function/"+r+"/"+e])

    # fifth: the classifiers ---------------------------------------------------
    classifiers = dict()
    for c in f['classifiers']:
        if not (isinstance(f["classifiers/"+c][0], str)): # if it is a string (else statement), then it means it was not a classifier
            classifiers[c] = np.array(f["classifiers/"+c],dtype = np.float64)
        else:
            classifiers[c] = "no_negative_examples"

    f.close()

    return hmm_file.name, use_cmalign, taxonomy, tax_function, classifiers




#===============================================================================
#                                      MAIN
#===============================================================================

def classify(database, aligned_sequences, verbose, output):
    # load the database
    hmm_file_path, use_cmalign, taxonomy, tax_function, classifiers = load_DB(database)

    # analyse the aligned sequences

    # delete the hmm temp file that was created
    os.remove(hmm_file_path)
