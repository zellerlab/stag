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

from stag import __version__ as tool_version
import stag.align as align
from stag.load_db import load_db

#===============================================================================
#                            LOAD THE HDF5 DATABASE
#===============================================================================
#def load_DB(hdf5_DB_path, protein_fasta_input, aligned_sequences):
#    f = h5py.File(hdf5_DB_path, 'r')
#
#    # zero: tool version -------------------------------------------------------
#    db_tool_version = f['tool_version'][0]
#    # check that it is the correct database, for 'classify', we need a single
#    # gene
#    if f['db_type'][0] != "single_gene":
#        sys.stderr.write("[E::main] Error: this database is not designed to run with stag classify\n")
#        sys.exit(1)
#    # check if we used proteins
#    if not aligned_sequences:
#        if protein_fasta_input and not f['align_protein'][0]:
#            # some proteins are provided in the classify but the db was constructed without using the proteins
#            sys.stderr.write("Error: protein provided, but the database was constructed on genes.\n")
#            sys.exit(1)
#        elif not protein_fasta_input and f['align_protein'][0]:
#            # the classify do not have proteins but the db was constructed WITH the proteins
#            sys.stderr.write("Error: missing protein file (the database was constructed aligning proteins).\n")
#            sys.exit(1)
#
#    # first, we save a temporary file with the hmm file ------------------------
#    with tempfile.NamedTemporaryFile(delete=False, mode="w") as hmm_file:
#        os.chmod(hmm_file.name, 0o644)
#        hmm_file.write(f['hmm_file'][0])
#        hmm_file.flush()
#        os.fsync(hmm_file.fileno())
#
#    # second if we need to use cm_align ----------------------------------------
#    use_cmalign = f['use_cmalign'][0] # bool
#
#    # third: taxonomy ----------------------------------------------------------
#    taxonomy = {key: list(f['taxonomy/{}'.format(key)]) for key in f['taxonomy']}
#
#    # fourth: tax_function -----------------------------------------------------
#    tax_function = {str(key): np.array(f['tax_function/{}'.format(key)], dtype=np.float64)
#                    for key in f['tax_function']}
#
#
#    # fifth: the classifiers ---------------------------------------------------
#    classifiers = dict()
#    for key in f['classifiers']:
#        classifier = f['classifiers/{}'.format(key)]
#        if not isinstance(classifier[0], str):
#            classifiers[key] = np.array(classifier, dtype=np.float64)
#        else:
#            classifiers[key] = "no_negative_examples"
#
#    f.close()
#
#    return hmm_file.name, use_cmalign, taxonomy, tax_function, classifiers, db_tool_version

#===============================================================================
#                    FUNCTION TO LOAD ALIGNED SEQUENCES
#===============================================================================
def alignment_reader(aligned_sequences):
    with open(aligned_sequences,"r") as align_in:
        for ali_line in align_in:
            gene_id = ali_line[:ali_line.find("\t")]
            aligned_seq = [int(col) for col in ali_line.split("\t")[1:].rstrip()]
            yield gene_id, np.array(aligned_seq, dtype=bool)

#===============================================================================
#                     TAXONOMICALLY ANNOTATE SEQUENCES
#===============================================================================
def run_logistic_prediction(seq, coeff_raw):
    # the first value of the coeff is the intercept
    coeff = coeff_raw[1:]
    intercept = coeff_raw[0]
    # calculate
    sm = coeff * seq
    np_sum = (sm).sum() + intercept
    score = 1 / (1 + np.exp(-np_sum))
    return score

# given many taxa (all siblings) and a sequence, it finds taxa with the highest
# score. Returns the taxa name and the score
def find_best_score(test_seq, siblings, classifiers):
    best_score = -1
    best_taxa = ""
    # check that siblings is not empty:
    if len(siblings) < 1:
        sys.stderr.write("Error. no siblings")
    # if there is only one sibiling:
    if len(siblings) == 1:
        best_score = 2 # if there are no siblings I put 2, it will be replaced after
        best_taxa = siblings[0]
    if len(siblings) > 1:
        for s in siblings:
            this_score = run_logistic_prediction(test_seq, classifiers[s])
            if this_score > best_score:
                best_score = this_score
                best_taxa = s
    return best_taxa, best_score

def predict_iter(test_seq, taxonomy, classifiers, tax, perc, arrived_so_far):
    # last iterative step
    if taxonomy.get(arrived_so_far) is not None:
        t, p = find_best_score(test_seq, taxonomy[arrived_so_far], classifiers)
        tax.append(t)
        perc.append(p)
        predict_iter(test_seq, taxonomy, classifiers, tax, perc, t)


#===============================================================================
#                    FIND TO WHICH TAXONOMIC LEVEL TO STOP
#===============================================================================
def find_correct_level(perc, tax_function):
    prob_per_level = list()
    max_probability = 0
    selected_level = -1
    for level in sorted(list(tax_function)):
        seq = np.asarray(perc)
        prob_this_level = run_logistic_prediction(seq, tax_function[level])
        if prob_this_level > max_probability:
            max_probability = prob_this_level
            selected_level = level
        prob_per_level.append(level+":"+str(prob_this_level))

    # note that selected_level represents which level to predict, so 2 means that we
    # predict 0,1,2. If it is "-1", then we cannot predict anything
    return selected_level, prob_per_level

#===============================================================================
#              FIND THE NUMBER OF MATCH FROM AN ALIGNED SEQUENCE
#===============================================================================
# test_seq is a numpy array, example:
# [False,  True, False,  True, False,  True, False,  True, False]
def find_n_aligned_characters(test_seq):
    # we go by multiple of five, for the one-hot encoding.
    # hence, 0\t0\t0\t0\t1 corresponds to "A". If a character doesnt match to an
    # internal state, then we have: "1\t0\t0\t0\t0".
    # Hence, it's enough to check positions 0,5,10,... and if there is a 0 (False)
    # it means that there was a match to an internal state
    pos_0_array = test_seq[0::5]
    n_False = np.size(pos_0_array) - np.sum(pos_0_array)
    return n_False


#===============================================================================
#                             CLASSIFY ONE SEQUENCE
#===============================================================================

def classify_seq(gene_id, test_seq, taxonomy, tax_function, classifiers, threads, verbose):
    # al_seq is a dictionary with one element, example:
    # {'gene1': array([False,  True, False,  True, False,  True, False,  True, False])}

    #gene_id, test_seq = list(al_seq.items())[0]
    # number of characters that map to the internal states of the HMM
    n_aligned_characters = find_n_aligned_characters(test_seq)

    # now we evaluate across the taxonomy --------------------------------------
    tax, perc = list(), list()
    # we arrived at the root, and now we classify from there
    predict_iter(test_seq, taxonomy, classifiers, tax, perc, "tree_root")

    # we change the predictions that came from having only one sibiling --------
    if perc[0] == 2:
        perc[0] = 1
    for i in range(len(perc)):
        if perc[i] == 2:
            perc[i] = perc[i-1]

    # now we have the raw prediction, we compare to  ---------------------------
    # the empirical values to obtain a better result
    selected_level, prob_per_level = find_correct_level(perc, tax_function)

    prob_per_level = "/".join(prob_per_level)
    perc_text = "\t".join([str(p) for p in perc])
    assigned_tax_text = ";".join(tax[0:(int(selected_level) + 1)])
    tax_text = "/".join(tax)

    result = [gene_id, assigned_tax_text, tax_text, selected_level,
              perc_text, prob_per_level, str(n_aligned_characters)]

    return result


#===============================================================================
#                                      MAIN
#===============================================================================

def classify(database, fasta_input=None, protein_fasta_input=None, verbose=3, threads=1, output=None,
             long_out=False, current_tool_version=tool_version,
             aligned_sequences=None, save_ali_to_file=None, min_perc_state=0, internal_call=False):
    t0 = time.time()
    # load the database
    db = load_db(database, protein_fasta_input=protein_fasta_input, aligned_sequences=aligned_sequences)
    hmm_file_path, use_cmalign, taxonomy, tax_function, classifiers, db_tool_version = db
    if verbose > 2:
        time_after_loading = time.time()
        sys.stderr.write("Load database: " + str("{0:.2f}".format(time_after_loading - t0))+" sec\n")

    alignment_length = None
    # align the sequences and classify them
    list_to_print = list()

    alignment_out, write_alignments = contextlib.nullcontext(), False
    if aligned_sequences:
        alignments = alignment_reader(aligned_sequences)
    else:
        alignments = align.align_generator(fasta_input, protein_fasta_input, hmm_file_path, use_cmalign,
                                           threads, verbose, True, min_perc_state)
        if save_ali_to_file:
            alignment_out, write_alignments = open(save_ali_to_file, "w"), True

    with alignment_out:
        for gene_id, ali in alignments:
            if not alignment_length:
                alignment_length = len(ali)
            list_to_print.append(classify_seq(gene_id, ali, taxonomy, tax_function, classifiers, threads, verbose))

            if write_alignments:
                ali_str = np.char.mod('%.0f', ali)
                print(gene_id, *ali_str, sep="\t", file=alignment_out)

    if verbose > 2:
        time_after_classification = time.time()
        sys.stderr.write("Classify sequences: " + str("{0:.2f}".format(time_after_classification - time_after_loading))+" sec\n")

    # delete the hmm temp file that was created --------------------------------
    os.remove(hmm_file_path)

    # print the sequences ------------------------------------------------------
    if output:
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(outfile.name, 0o644)
    else:
        outfile = sys.stdout


    out_header = ["sequence", "taxonomy", "full_taxonomy", "selected_level",
                  "prob_from_classifiers", "prob_per_level", "n_aligned_characters"]

    if not long_out or internal_call:
        out_header = out_header[:2]
        list_to_print = [item[:2] for item in list_to_print]

    if not internal_call:
        with outfile:
            print(*out_header, sep="\t", file=outfile)
            for line in list_to_print:
                print(*line, sep="\t", file=outfile)

            if output:
                try:
                    outfile.flush()
                    os.fsync(outfile.fileno())
                except:
                    sys.stderr.write("[E::main] Error: failed to save the result\n")
                    sys.exit(1)
                try:
                    #os.rename(outfile.name,output) # atomic operation
                    shutil.move(outfile.name, output) #It is not atomic if the files are on different filsystems.
                except:
                    sys.stderr.write("[E::main] Error: failed to save the profile\n")
                    sys.stderr.write("[E::main] you can find the file here:\n" + outfile.name + "\n")
                    sys.exit(1)

    return alignment_length, list_to_print
