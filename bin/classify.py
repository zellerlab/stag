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
def load_DB(hdf5_DB_path, protein_fasta_input):
    f = h5py.File(hdf5_DB_path, 'r')

    # zero: tool version -------------------------------------------------------
    db_tool_version = f['tool_version'][0]
    # check that it is the correct database, for 'classify', we need a single
    # gene
    if f['db_type'][0] != "single_gene":
        sys.stderr.write("[E::main] Error: this database is not designed to run with stag classify\n")
        sys.exit(1)
    # check if we used proteins
    if not(protein_fasta_input is None):
        # some proteins are provided in the classify
        if not f['align_protein'][0]:
            # but the db was constructed without using the proteins
            sys.stderr.write("Error: protein provided, but the database was constructed on genes.\n")
            sys.exit(1)
    else:
        # the classify do not have proteins
        if f['align_protein'][0]:
            # but the db was constructed WITH the proteins
            sys.stderr.write("Error: missing protein file (the database was constructed aligning proteins).\n")
            sys.exit(1)



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
        tax_function[str(c)] = np.array(f["tax_function/"+c],dtype = np.float64)

    # fifth: the classifiers ---------------------------------------------------
    classifiers = dict()
    for c in f['classifiers']:
        if not (isinstance(f["classifiers/"+c][0], str)): # if it is a string (else statement), then it means it was not a classifier
            classifiers[c] = np.array(f["classifiers/"+c],dtype = np.float64)
        else:
            classifiers[c] = "no_negative_examples"

    f.close()

    return hmm_file.name, use_cmalign, taxonomy, tax_function, classifiers, db_tool_version


#===============================================================================
#                    FUNCTION TO LOAD ALIGNED SEQUENCES
#===============================================================================
def file_2_generator(aligned_sequences):
    o = open(aligned_sequences,"r")
    for i in o:
        merged_fasta = i.rstrip()
        gene_id = merged_fasta.split("\t")[0]
        converted_ali = list()
        for u in merged_fasta.split("\t")[1:]:
            converted_ali.append(int(u))
        to_return = dict()
        to_return[gene_id] = np.array(converted_ali,dtype=bool)
        yield to_return
    o.close()



#===============================================================================
#                     TAXONOMICALLY ANNOTATE SEQUENCES
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
        best_score = 2 # if there are no sibilings I put 2, it will be replaced after
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
#                    FIND TO WHICH TAXONOMIC LEVEL WE STOP
#===============================================================================
def run_prediction_no_penalty(seq, coeff_raw):
    # the first value of the coeff is the intercept
    coeff = coeff_raw[1:]
    intercept = coeff_raw[0]
    # calculate
    sm = coeff*seq
    np_sum = (sm).sum() + intercept
    score = 1/(1+np.exp(-np_sum))
    return score

def find_correct_level(perc, tax_function):
    prob_per_level = list()
    max_v = 0
    sel_lev = -1
    for l in sorted(list(tax_function)):
        seq = np.asarray(perc)
        prob_this_level = run_prediction_no_penalty(seq, tax_function[l])
        if prob_this_level > max_v:
            max_v = prob_this_level
            sel_lev = l
        prob_per_level.append(l+":"+str(prob_this_level))

    # note that sel_lev represents which level to predict, so 2 means that we
    # predict 0,1,2. If it is "-1", then we cannot predict anything
    return sel_lev, prob_per_level


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

    # we change the predictions that came from having only one sibiling --------
    if perc[0] == 2:
        perc[0] = 1
    for i in range(len(perc)):
        if perc[i] == 2:
            perc[i] = perc[i-1]

    # now we have the raw prediction, we compare to  ---------------------------
    # the empirical values to obtain a better result
    sel_lev, prob_per_level = find_correct_level(perc, tax_function)

    # transform perc to string
    perc_text = list()
    for i in perc:
        perc_text.append(str(i))

    # return the result --------------------------------------------------------
    res_string = res_string + "\t" + ";".join(tax[0:(int(sel_lev)+1)]) + "\t" + "/".join(tax) + "\t" + sel_lev + "\t" + "/".join(perc_text) + "\t" + "/".join(prob_per_level)
    return res_string



#===============================================================================
#                                      MAIN
#===============================================================================

def classify(database, fasta_input, protein_fasta_input, verbose, threads, output, long_out, current_tool_version, aligned_sequences):
    t0 = time.time()
    # load the database
    hmm_file_path, use_cmalign, taxonomy, tax_function, classifiers, db_tool_version = load_DB(database, protein_fasta_input)
    if verbose>2:
        time_after_loading = time.time()
        sys.stderr.write("Load database: " + str("{0:.2f}".format(time_after_loading - t0))+" sec\n")

    # align the sequences and classify them
    list_to_print = list()
    if aligned_sequences is None:
        for al_seq in align.align_generator(fasta_input,protein_fasta_input,hmm_file_path, use_cmalign, threads, verbose, True):
            list_to_print.append(classify_seq(al_seq, taxonomy, tax_function, classifiers, threads, verbose))
    else:
        for al_seq in file_2_generator(aligned_sequences):
            list_to_print.append(classify_seq(al_seq, taxonomy, tax_function, classifiers, threads, verbose))


    if verbose>2:
        time_after_classification = time.time()
        sys.stderr.write("Classify sequences: " + str("{0:.2f}".format(time_after_classification - time_after_loading))+" sec\n")

    # delete the hmm temp file that was created --------------------------------
    os.remove(hmm_file_path)

    # print the sequences ------------------------------------------------------
    if not(output is None):
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(outfile.name, 0o644)
    else:
        outfile = sys.stdout

    if long_out:
        outfile.write("sequence\ttaxonomy\tfull_taxonomy\tselected_level\tprob_from_classifiers\tprob_per_level\n")
    else:
        outfile.write("sequence\ttaxonomy\n")

    for i in list_to_print:
        if long_out:
            outfile.write(i+"\n")
        else:
            outfile.write("\t".join(i.split("\t")[0:2])+"\n")

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
