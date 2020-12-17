"""
Scripts that saves the database to a directory
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

# Requirements:
# - numpy
# - h5py

import numpy as np
import sys
import time
import os
import h5py
import tempfile
import shutil

#===============================================================================
#                            LOAD THE HDF5 DATABASE
#===============================================================================
def load_DB(hdf5_DB_path):
    f = h5py.File(hdf5_DB_path, 'r')

    # tool version -------------------------------------------------------------
    db_tool_version = f['tool_version'][0]

    # use proteins -------------------------------------------------------------
    use_proteins = f['align_protein'][0] # bool

    # second if we need to use cm_align ----------------------------------------
    use_cmalign = f['use_cmalign'][0] # bool

    # we save a temporary file with the hmm file ------------------------
    hmm_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(hmm_file.name, 0o644)
    hmm_file.write(f['hmm_file'][0])
    hmm_file.flush()
    os.fsync(hmm_file.fileno())
    hmm_file.close()

    # third: taxonomy ----------------------------------------------------------
    taxonomy = dict()
    for i in f['taxonomy']:
        taxonomy[i] = list(f['taxonomy/'+i])

    # fourth: tax_function -----------------------------------------------------
    tax_function = dict()
    for c in f['tax_function']:
        tax_function[str(c)] = np.array(f["tax_function/"+c],dtype = np.float64)

    # fifth: the classifiers ---------------------------------------------------
    classifiers = dict()
    for c in f['classifiers']:
        if not (isinstance(f["classifiers/"+c][0], str)): # if it is a string (else statement), then it means it was not a classifier
            classifiers[c] = np.array(f["classifiers/"+c],dtype = np.float64)
        else:
            classifiers[c] = "no_negative_examples"

    f.close()

    return hmm_file.name, db_tool_version, use_proteins, use_cmalign, taxonomy, tax_function, classifiers

#===============================================================================
#                                      MAIN
#===============================================================================
def unzip_db(database, verbose, dir_output):
    # load the database
    hmm_file, db_tool_version, use_proteins, use_cmalign, taxonomy, tax_function, classifiers = load_DB(database)

    # check if the output dir exists already
    if os.path.isdir(dir_output):
        sys.stderr.write("Error, output dir exists already\n")
        sys.exit(1)
