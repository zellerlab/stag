"""
Scripts that find the taxonomy of an aligned sequence
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
import contextlib


import stag.align as align

#===============================================================================
#                            LOAD THE HDF5 DATABASE
#===============================================================================
def load_db(hdf5_DB_path, protein_fasta_input=None, aligned_sequences=None, dir_output=None):
    with h5py.File(hdf5_DB_path, 'r') as db_in:
        params_out = open(os.path.join(dir_output, "parameters.tsv"), "w") if dir_output else contextlib.nullcontext()
        with params_out:
            # zero: tool version -------------------------------------------------------
            db_tool_version = db_in['tool_version'][0]
            use_proteins = db_in['align_protein'][0] # bool
            use_cmalign = db_in['use_cmalign'][0] # bool
            if dir_output:
                params_out.write("Tool version: "+str(db_tool_version)+"\n")
                params_out.write("Use proteins for the alignment: "+str(use_proteins)+"\n")
                params_out.write("Use cmalign instead of hmmalign: "+str(use_cmalign)+"\n")

        # check that it is the correct database, for 'classify', we need a single
        # gene
        if db_in['db_type'][0] != "single_gene":
            sys.stderr.write("[E::main] Error: this database is not designed to run with stag classify\n")
            sys.exit(1)
        # check if we used proteins
        if not aligned_sequences and not dir_output:
            if protein_fasta_input and not db_in['align_protein'][0]:
                # some proteins are provided in the classify but the db was constructed without using the proteins
                sys.stderr.write("Error: protein provided, but the database was constructed on genes.\n")
                sys.exit(1)
            elif not protein_fasta_input and db_in['align_protein'][0]:
                # the classify do not have proteins but the db was constructed WITH the proteins
                sys.stderr.write("Error: missing protein file (the database was constructed aligning proteins).\n")
                sys.exit(1)
    
        # first, we save a temporary file with the hmm file ------------------------
        hmm_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
        with hmm_file:
            os.chmod(hmm_file.name, 0o644)
            hmm_file.write(db_in['hmm_file'][0])
            hmm_file.flush()
            os.fsync(hmm_file.fileno())

        if dir_output:
            shutil.move(hmm_file.name, os.path.join(dir_output, "hmmfile.hmm"))

    
        # second if we need to use cm_align ----------------------------------------
        use_cmalign = db_in['use_cmalign'][0] # bool
    
        # third: taxonomy ----------------------------------------------------------
        taxonomy = {key: list(db_in['taxonomy/{}'.format(key)]) for key in db_in['taxonomy']}
        if dir_output:
            tax_out = open(os.path.join(dir_output, "node_hierarchy.tsv"), "w")
            with tax_out:
                print("Node", "Children", sep="\t", file=tax_out)
                for key, values in taxonomy.items():
                    print(key, *map(str, values), sep="\t", file=tax_out)
    
        # fourth: tax_function -----------------------------------------------------
        tax_function = {str(key): np.array(db_in['tax_function/{}'.format(key)], dtype=np.float64) 
                        for key in db_in['tax_function']}
        if dir_output:
            tax_func_out = open(os.path.join(dir_output, "taxonomy_function.tsv"), "w")
            with tax_func_out:
                for key, value in tax_function.items():
                    print(key, value, sep="\t", file=tax_func_out)
            
    
        # fifth: the classifiers ---------------------------------------------------
        classifiers = dict()
        class_out = open(os.path.join(dir_output, "classifiers_weights.tsv"), "w") if dir_output else contextlib.nullcontext()
        with class_out:
            for key in db_in['classifiers']:
                classifier = db_in['classifiers/{}'.format(key)]
                if not isinstance(classifier[0], str):
                    classifiers[key] = np.array(classifier, dtype=np.float64) 
                else:
                    classifiers[key] = "no_negative_examples"
                if dir_output:
                    print(key, *classifiers[key], sep="\t", file=class_out)
    

    return hmm_file.name, use_cmalign, taxonomy, tax_function, classifiers, db_tool_version
