import sys
import time
import os
import tempfile
import shutil
import contextlib

import numpy as np
import pandas as pd
import h5py

from . import __version__ as tool_version
from stag.taxonomy3 import Taxonomy
from stag.databases2 import load_db
import stag.align as align

def alignment_reader(aligned_sequences):
    with open(aligned_sequences,"r") as align_in:
        for ali_line in align_in:
            gene_id, *aligned_seq = ali_line.rstrip().split("\t")
            yield gene_id, np.array(list(map(int, aligned_seq)), dtype=bool)

def run_logistic_prediction(seq, coeff_raw):
    # the first value of the coeff is the intercept
    intercept, *coeff = coeff_raw
    sm = coeff * seq
    np_sum = (sm).sum() + intercept
    return 1 / (1 + np.exp(-np_sum))

# given many taxa (all siblings) and a sequence, it finds taxa with the highest
# score. Returns the taxa name and the score
def find_best_score(test_seq, siblings, classifiers):
    best_score, best_taxa = -1, ""
    if not siblings:
        pass
    elif len(siblings) == 1:
        # if there are no siblings I put 2, it will be replaced after
        best_score, best_taxa = 2, siblings[0]
    else:
        for sibling in siblings:
            this_score = run_logistic_prediction(test_seq, classifiers[sibling])
            if this_score > best_score:
                best_score, best_taxa = this_score, sibling
    return best_taxa, best_score

def predict_iter(test_seq, taxonomy, classifiers, tax, perc, arrived_so_far):
    # last iterative step
    if taxonomy.get(arrived_so_far) is not None:
        t, p = find_best_score(test_seq, taxonomy[arrived_so_far], classifiers)
        if t:
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
    # test_seq is a boolean numpy array corresponding to the encoded aligned sequence
    #print("TAX", taxonomy)
    # number of characters that map to the internal states of the HMM
    n_aligned_characters = find_n_aligned_characters(test_seq)

    # now we evaluate across the taxonomy --------------------------------------
    tax, perc = list(), list()
    # we arrived at the root, and now we classify from there
    predict_iter(test_seq, taxonomy, classifiers, tax, perc, Taxonomy.TREE_ROOT)

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
    perc_text = "/".join([str(p) for p in perc])
    assigned_tax_text = ";".join(tax[0:(int(selected_level) + 1)])
    tax_text = "/".join(tax)

    result = [gene_id, assigned_tax_text, tax_text, selected_level,
              perc_text, prob_per_level, str(n_aligned_characters)]

    return result


#===============================================================================
#                            DO NN CLASSIFICATION
#===============================================================================
# we change `list_to_print` directly inside the function adding the NN
# classification. Note that we have to change the position -1
def NN_classification(list_to_print,db,ali,base_save_features):
    annotation_this = "" # annotation based on the NN and the thresholds
    NN_this = "" # nearest neighbour sequence
    min_this = float('inf') # distance to the NN_this

    # --------------------------------------------------------------------------
    # check if there is a family annotation (note it can be different than
    # family the selected level)

    # Note that we check the full taxonomy, so we always have a taxonomy annotation
    current_annotation = list_to_print[-1][2].split("/") # changed from: list_to_print[-1][1]

    NN_start_level = db["NN_start_level"] # it's the `-L` option, default is `4` (family)
    if len(current_annotation) <= NN_start_level:
        # we do not have a family annotation, hence we don't do NN classification
        min_this = "-1"
    else:
        # there is a family annotation
        family_this = current_annotation[NN_start_level]
        lmnn_this = db["LMNN"][family_this]
        # the alignment was reduced, keeping only some columns (the one that do change)
        # hence we select the same columns
        col_to_keep = db["all_sel_positions"][family_this]
        ali_sel = ali[col_to_keep]
        # if we have to save the features, we do it now:
        if(base_save_features is not None):
            sequence_ID = list_to_print[-1][0]
            pd.DataFrame(np.array([ali_sel*1]), index=[sequence_ID]).to_csv(base_save_features+"/"+family_this+"_"+sequence_ID+'_original.gz', compression= 'gzip')
        # we need to check if there was enough data to do the training
        # if there is enough data, then we do LMNN transformation, otherwise no
        if isinstance(lmnn_this,(np.ndarray)):
        # same as: if lmnn_this != "NOT ENOUGH DATA":
            # we transform the vector. This is exactly the same as in the training
            # after transformation
            ali_sel = ali_sel.dot(lmnn_this)
            if(base_save_features is not None):
                sequence_ID = list_to_print[-1][0]
                pd.DataFrame(np.array([ali_sel]), index=[sequence_ID]).to_csv(base_save_features+"/"+family_this+"_"+sequence_ID+'_transformed.gz', compression= 'gzip')

        # calculate the distances to all centroid sequences and find the minimum
        for species in db["centroid_seq"][family_this].index.values:
            distance_this = dist_vectors(db["centroid_seq"][family_this],species,ali_sel)
            if distance_this < min_this:
                min_this = distance_this
                NN_this = species

        # now we found the closest sequence, and we can use the threshold to decide
        # the NN taxonomy
        sel_level_nn = -1
        for th in db["thresholds_NN"][family_this]:
            if th > sel_level_nn:
                if db["thresholds_NN"][family_this][th] > min_this:
                    sel_level_nn = th
        # prepare the new taxonomy
        if sel_level_nn == -1:
            annotation_this = ""
        else:
            annotation_this = db["species_to_tax"][NN_this][0:sel_level_nn+1]
    # --------------------------------------------------------------------------
    # change the `list_to_print`
    list_to_print[-1].append(";".join(annotation_this))
    list_to_print[-1].append(NN_this)
    list_to_print[-1].append(str(min_this))

# Calculate perc identity between two sequences
#  - ALI is the pandas DF with the centroids
#  - pos is the centroid we are cheking
#  - ali_sel is the sequence we are cheking
def dist_vectors(ALI,pos,ali_sel):
    seq1 = ALI.loc[pos,].to_numpy()
    # euclidean distance
    dist = np.linalg.norm(seq1-ali_sel)
    return dist


#===============================================================================
#                                      MAIN
#===============================================================================

def classify(database, fasta_input=None, protein_fasta_input=None, verbose=3, threads=1, output=None,
             long_out=False, current_tool_version=tool_version,
             aligned_sequences=None, save_ali_to_file=None, min_perc_state=0, base_save_features=None, internal_call=False):
    t0 = time.time()

    # ==========================================================================
    # load database
    db = load_db(database, protein_fasta_input=protein_fasta_input, aligned_sequences=aligned_sequences)
    # db is now a dictionary with all the data inside
    if verbose > 2:
        time_after_loading = time.time()
        sys.stderr.write("Load database: " + str("{0:.2f}".format(time_after_loading - t0))+" sec\n")

    # ==========================================================================
    # align sequences
    alignment_length = None
    # align the sequences and classify them
    list_to_print = list()

    alignment_out, write_alignments = contextlib.nullcontext(), False
    if aligned_sequences:
        alignments = alignment_reader(aligned_sequences)
    else:
        alignments = align.align_generator(fasta_input, protein_fasta_input, db["hmm_file_path"], db["use_cmalign"],
                                           threads, verbose, True, min_perc_state)
        if save_ali_to_file:
            alignment_out, write_alignments = open(save_ali_to_file, "w"), True

    # ==========================================================================
    # classify sequences along the tree
    with alignment_out:
        for gene_id, ali in alignments:
            if not alignment_length:
                alignment_length = len(ali)
            list_to_print.append(classify_seq(gene_id, ali, db["taxonomy"], db["tax_function"], db["classifiers"], threads, verbose))

            if write_alignments:
                ali_str = np.char.mod('%.0f', ali)
                print(gene_id, *ali_str, sep="\t", file=alignment_out)

            # we do NN classification per sample (not ideal, would be better to collect all the alignment from one family)
            NN_classification(list_to_print,db,ali,base_save_features)

    if verbose > 2:
        time_after_classification = time.time()
        sys.stderr.write("Classify sequences: " + str("{0:.2f}".format(time_after_classification - time_after_loading))+" sec\n")

    # delete the hmm temp file that was created --------------------------------
    os.remove(db["hmm_file_path"])

    # ==========================================================================
    # print sequences
    if output:
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(outfile.name, 0o644)
    else:
        outfile = sys.stdout

    out_header = ["sequence", "taxonomy", "full_taxonomy", "selected_level",
                  "prob_from_classifiers", "prob_per_level", "n_aligned_characters",
                  "annotation_NN","NN","distance"]

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
