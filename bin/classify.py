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

# load align routine -----------------------------------------------------------
path_this = os.path.realpath(__file__)
path_array = path_this.split("/")
relative_path = "/".join(path_array[0:-1])
# add /bin to the path ---------------------------------------------------------
try:
    if os.path.isdir(relative_path):
        sys.path.insert(0, relative_path)
    else:
        sys.stderr.write("[E::main] Error: "+relative_path+" directory is missing.\n")
        sys.exit(1)
except:
    sys.stderr.write("[E::main] Error: "+relative_path+" directory is missing.\n")
    sys.exit(1)

try:
    import align as align
except:
    sys.stderr.write("[E::main] Error: fail to load the script: "+relative_path+"/align.py\n")
    sys.exit(1)

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
    use_cmalign = f['use_cmalign'][0] # bool

    # third: taxonomy ----------------------------------------------------------
    taxonomy = dict()
    for i in f['taxonomy']:
        taxonomy[i] = list(f['taxonomy/'+i])

    # fourth: tax_function -----------------------------------------------------
    tax_function = dict()
    for c in f['tax_function']:
        tax_function[int(c)] = np.array(f["tax_function/"+c],dtype = np.float64)

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
#                       TAXONOMY ANNOTATE SEQUENCES
#===============================================================================
def run_lasso_prediction(seq, coeff):
    sm = coeff*seq
    np_sum = (sm).sum()
    score = 1/(1+np.exp(-np_sum))
    return score

# given many taxa (all sibilings) and a sequence, it finds taxa with the highest
# score. Returns the taxa name and the score
def find_best_score(test_seq, sibilings, classifiers):
    best_score = -1
    best_taxa = ""
    # check that sibilings is not empty:
    if len(sibilings) < 1:
        sys.stderr.write("Error. no sibilings")
    # if there is only one sibiling:
    if len(sibilings) == 1:
        best_score = 1
        best_taxa = sibilings[0]
    if len(sibilings) > 1:
        for s in sibilings:
            this_score = run_lasso_prediction(test_seq, classifiers[s])
            if this_score > best_score:
                best_score = this_score
                best_taxa = s
    return best_taxa, best_score

def predict_iter(test_seq, taxonomy, classifiers, tax, perc, arrived_so_far):
    # last iterative step
    if arrived_so_far in taxonomy:
        t,p = find_best_score(test_seq, taxonomy[arrived_so_far], classifiers)
        tax.append(t)
        perc.append(p)
        predict_iter(test_seq, taxonomy, classifiers, tax, perc, t)


#===============================================================================
#                      CALCULATE EMPIRICAL PROBABILITY
#===============================================================================
def calc_empirical_vals(perc, tax_function):
    prob_per_level = list()
    max_v = 0
    sel_lev = -1
    for l in sorted(list(tax_function)):
        seq = np.asarray(perc)
        prob_this_level = run_lasso_prediction(seq, tax_function[l])
        if prob_this_level > max_v:
            max_v = prob_this_level
            sel_lev = l
        prob_per_level.append(str(prob_this_level))

    # note that sel_lev represents which level was removed, so
    # 0 means that we cannot assign even at the lowest level
    # (there is one predictor more than the number of levels)
    return sel_lev-1, prob_per_level


#===============================================================================
#                             CLASSIFY ONE SEQUENCE
#===============================================================================

def classify_seq(al_seq, taxonomy, tax_function, classifiers, threads, verbose):
    # al_seq is a dictionary with one element, example:
    # {'gene1': array([False,  True, False,  True, False,  True, False,  True, False])}

    # name of the gene
    res_string = list(al_seq.keys())[0]
    # sequence in numpy format
    test_seq = al_seq[res_string]

    # now we evaluate across the taxonomy --------------------------------------
    tax = list()
    perc = list()
    # we arrived at the root, and now we classify from there
    predict_iter(test_seq, taxonomy, classifiers, tax, perc, "tree_root")

    # now we have the raw prediction, we compare to  ---------------------------
    # the empirical values to obtain a better result
    sel_lev, prob_per_level = calc_empirical_vals(perc, tax_function)

    # transform perc to string
    perc_text = list()
    for i in perc:
        perc_text.append(str(i))

    # return the result --------------------------------------------------------
    res_string = res_string + "\t" + "/".join(tax[0:sel_lev+1]) + "\t" + "/".join(tax) + "\t" + str(sel_lev) + "\t" + "/".join(perc_text) + "\t" + "/".join(prob_per_level)
    return res_string



#===============================================================================
#                                      MAIN
#===============================================================================

def classify(database, fasta_input, protein_fasta_input, verbose, threads, output):
    # load the database
    hmm_file_path, use_cmalign, taxonomy, tax_function, classifiers = load_DB(database)

    # align the sequences and classify them
    list_to_print = list()
    for al_seq in align.align_generator(fasta_input,protein_fasta_input,hmm_file_path, use_cmalign, threads, verbose, True):
        list_to_print.append(classify_seq(al_seq, taxonomy, tax_function, classifiers, threads, verbose))

    # delete the hmm temp file that was created --------------------------------
    os.remove(hmm_file_path)

    # print the sequences ------------------------------------------------------
    if not(output is None):
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(outfile.name, 0o644)
    else:
        outfile = sys.stdout

    outfile.write("sequence\ttaxonomy\tfull_taxonomy\tselected_level\tprob_from_classifiers\tprob_per_level\n")
    for i in list_to_print:
        outfile.write(i+"\n")

    # close
    if not(output is None):
        try:
            outfile.flush()
            os.fsync(outfile.fileno())
            outfile.close()
        except:
            sys.stderr.write("[E::main] Error: failed to save the result\n")
            sys.exit(1)
        try:
            #os.rename(outfile.name,output) # atomic operation
            shutil.move(outfile.name,output) #It is not atomic if the files are on different filsystems.
        except:
            sys.stderr.write("[E::main] Error: failed to save the profile\n")
            sys.stderr.write("[E::main] you can find the file here:\n"+outfile.name+"\n")
            sys.exit(1)

    sys.exit(0) # correct
