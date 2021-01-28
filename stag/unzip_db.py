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

from stag.load_db import load_db
#===============================================================================
#                            LOAD THE HDF5 DATABASE
#===============================================================================
#def load_db(hdf5_DB_path, protein_fasta_input, aligned_sequences, dir_output=None):
#def load_DB(hdf5_DB_path, verbose, dir_output):
#    f = h5py.File(hdf5_DB_path, 'r')
#
#    o = open(dir_output+"/parameters.tsv","w")
#    # tool version -------------------------------------------------------------
#    db_tool_version = f['tool_version'][0]
#    o.write("Tool version: "+str(db_tool_version)+"\n")
#
#    # use proteins -------------------------------------------------------------
#    use_proteins = f['align_protein'][0] # bool
#    o.write("Use proteins for the alignment: "+str(use_proteins)+"\n")
#
#    # second if we need to use cm_align ----------------------------------------
#    use_cmalign = f['use_cmalign'][0] # bool
#    o.write("Use cmalign instead of hmmalign: "+str(use_cmalign)+"\n")
#    o.close()
#
#    # we save a temporary file with the hmm file ------------------------
#    hmm_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
#    os.chmod(hmm_file.name, 0o644)
#    hmm_file.write(f['hmm_file'][0])
#    hmm_file.flush()
#    os.fsync(hmm_file.fileno())
#    hmm_file.close()
#    shutil.move(hmm_file.name,dir_output+"/hmmfile.hmm")
#
#    # third: taxonomy ----------------------------------------------------------
#    o = open(dir_output+"/node_hierarchy.tsv","w")
#    o.write("Node\tChildren\n")
#    for i in f['taxonomy']:
#        o.write(i+"\t"+";".join(list(f['taxonomy/'+i]))+"\n")
#    o.close()
#
#    # fourth: tax_function -----------------------------------------------------
#    o = open(dir_output+"/taxonomy_function.tsv","w")
#    for c in f['tax_function']:
#        val_this = np.array(f["tax_function/"+c],dtype = np.float64)
#        o.write(str(c)+"\t"+ str(val_this)+"\n")
#    o.close()
#
#    # fifth: the classifiers ---------------------------------------------------
#    o = open(dir_output+"/classifiers_weights.tsv","w")
#    for c in f['classifiers']:
#        if not (isinstance(f["classifiers/"+c][0], str)): # if it is a string (else statement), then it means it was not a classifier
#            this_class = np.array(f["classifiers/"+c],dtype = np.float64)
#        else:
#            this_class = "no_negative_examples"
#        o.write(str(c)+"\t"+ str(this_class)+"\n")
#
#    o.close()
#    f.close()

#===============================================================================
#                                      MAIN
#===============================================================================
def unzip_db(database, verbose, dir_output):
    # check if the output dir exists already
    if os.path.isdir(dir_output):
        sys.stderr.write("Error, output dir exists already\n")
        sys.exit(1)
    else:
        try:
            os.mkdir(dir_output)
        except OSError:
            sys.stderr.write("Error, failed to create the output directory\n")
            sys.exit(1)

    # load the database and save it
    load_db(database, dir_output=dir_output)
