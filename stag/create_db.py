"""
Scripts that creates the database of classifiers
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

# Input:
#  - one multiple sequence alignment (MSA) per marker gene. The MSA is obtained
#    from the function stag align, like:
#       >gene1\t0001010001010000100101000...
#       >gene2\t0000110001010100100101001...
#  - a taxonomy file that describes the taxonomy of the genes:
#       gene1\tBacteria;Firmicutes;...
#
# Output:
#  - a database file (hdf5) that can be used by stag classify

import numpy as np
import sys
import random
import time
import pandas as pd
import logging
import os
import math
from sklearn.linear_model import LogisticRegression
from sklearn.datasets import make_classification
import h5py
import tempfile
import shutil

from stag.taxonomy import Taxonomy

# Function to identify the rownames and number of columns in an alignment
def find_raw_names_ncol(file_name):
    gene_names = list()
    with open(file_name, "r") as f:
        for line in f.readlines():
            gene_names.append(line.split("\t")[0])
        n_col = len(line.split("\t"))
    return gene_names, n_col

# function to load an alignment produced by the "align" option =================
# Input:
#  - a file created by "align"
# Output:
#  - a panda object
# as a note, numpy.loadtxt is way slower than pandas read.csv
# It works also on .gz files
def load_alignment_from_file(file_name):
    # create empty pandas object of the correct size
    gene_names,ncol = find_raw_names_ncol(file_name)
    alignment = pd.DataFrame(False,index = gene_names,columns = range(ncol-1))
    # add correct values
    pos = 0
    with open(file_name, "r") as f:
        for line in f.readlines():
            vals = line.rstrip().split("\t")
            alignment.iloc[pos]= np.array([ False if x == "0" else True for x in vals[1:]])
            pos = pos + 1

    logging.info('   LOAD_AL: Number of genes: %s', str(len(list(alignment.index.values))))

    # we remove duplicates
    alignment = alignment.drop_duplicates()
    logging.info('   LOAD_AL: Number of genes, after removing duplicates: %s', str(len(list(alignment.index.values))))
    return alignment

# function to check that taxonomy and alignment are consistent =================
# 1. all genes in the alignment should be in the taxonomy,
# 2. the taxonomy can have more genes, than the one that are present in the
#    alignment, but we need to remove them, since the selection of the genes
#    for training and testing is done at the level of the taxonomy
def check_taxonomy_alignment_consistency(alignment, full_taxonomy):
    genes_in_alignment = list(alignment.index.values)
    genes_taxonomy = full_taxonomy.find_gene_ids(full_taxonomy.get_root())
    logging.info('   CHECK: genes in alignment: %s', str(len(genes_in_alignment)))
    logging.info('   CHECK: genes in taxonomy:  %s', str(len(genes_taxonomy)))

    # check that all genes in the alignment are in the taxonomy ----------------
    if not(set(genes_in_alignment).issubset(set(genes_taxonomy))):
        sys.stderr.write("Error: some genes in the alignment have no taxonomy.\n")
        sys.stderr.write("       Use the command 'check_input' to find more information.\n")
        logging.info(' Error: some genes in the alignment have no taxonomy.')
        for g in genes_in_alignment:
            if g not in genes_taxonomy:
                logging.info('    %s',g)
        sys.exit(1)
    else:
        logging.info('   CHECK: check all genes in the alignment have a taxonomy: correct')

    # check if we need to remove some genes from the taxonomy ------------------
    not_needed_gene_tax = set(genes_taxonomy).difference(set(genes_in_alignment))
    if len(not_needed_gene_tax) == 0:
        logging.info('   CHECK: check genes that we need to remove from the taxonomy: None')
    else:
        logging.info('   CHECK: check genes that we need to remove from the taxonomy: %s', str(len(not_needed_gene_tax)))
        full_taxonomy.remove_genes(list(not_needed_gene_tax))

    # double check that the number of genes is the same in the alignment and in
    # the taxonomy
    genes_taxonomy = full_taxonomy.find_gene_ids(full_taxonomy.get_root())
    if len(genes_taxonomy) != len(genes_in_alignment):
        sys.stderr.write("Error: even after correction, the genes in the taxonomy and the alignment do not agree\n")
        logging.info(' Error: even after correction, the genes in the taxonomy and the alignment do not agree.')
        sys.exit(1)








#===============================================================================
#                   FUNCTIONS TO TRAIN THE CLASSIFIERS
#===============================================================================

# function that finds positive and negative examples ===========================
def find_training_genes(node, sibilings, full_taxonomy, alignment):
    positive_examples = full_taxonomy.find_gene_ids(node)
    negative_examples = list()
    if len(sibilings) > 0:
        for s in sibilings:
            negative_examples = negative_examples + full_taxonomy.find_gene_ids(s)
    # "positive_examples" and "negative_examples" are list of gene ids

    if len(negative_examples) == 0:
        # it means that there was only one child, and we cannot do anything
        return positive_examples, negative_examples

    # From here, it means that there is at least one sibiling ==================
    # We make classes more balanced
    positive_examples_subsample = list(positive_examples)
    negative_examples_subsample = list(negative_examples)
    # 1. max 500 positive samples ----------------------------------------------
    if len(positive_examples_subsample) > 500:
        positive_examples_subsample = random.sample(positive_examples_subsample, 500)
    # 2. max 1000 negative samples ---------------------------------------------
    if len(negative_examples_subsample) > 1000:
        negative_examples_subsample = random.sample(negative_examples_subsample, 1000)
    # 3. max 20 times more negative than positive ------------------------------
    if len(negative_examples_subsample) > len(positive_examples_subsample)*20:
        negative_examples_subsample = random.sample(negative_examples_subsample, len(positive_examples_subsample)*20)
    # 4. we want to have at least 5 times more negative than positive ----------
    missing_neg = 0 # how many negative sequences we need to add
    if len(sibilings) == 1:
        # if there is only one other sibiling, we choose only 3 times more negative
        if len(negative_examples_subsample) > len(positive_examples_subsample)*3:
            negative_examples_subsample = random.sample(negative_examples_subsample, len(positive_examples_subsample)*3)
    if len(negative_examples_subsample) < len(positive_examples_subsample)*5:
        missing_neg = len(positive_examples_subsample)*5 - len(negative_examples_subsample)
    # add negative examples if necessary
    if missing_neg > 0:
        # positive examples
        X_clade = alignment.loc[positive_examples, : ].to_numpy()
        # always have 5 positive classes
        n_positive_class = len(X_clade)
        for i in range(n_positive_class,5):
            rr = random.choice(range(0,n_positive_class))
            X_clade = np.vstack((X_clade,X_clade[rr,]))

        # find possible genes to add additionaly to negarives
        possible_neg = list(set(alignment.index.values).difference(set(positive_examples + negative_examples)))
        if possible_neg: # if it is possible to add negatives
                         # note that at the highest level, it's not possible
            X_poss_na = alignment.loc[possible_neg, : ].to_numpy()
            len_poss_na = len(X_poss_na)

            # choose 5 random positive clades
            X_check_sim = X_clade[random.sample(range(len(X_clade)),5),]

            m_for_diff_0 = np.tile(X_check_sim[0,],(len_poss_na,1))
            m_for_diff_1 = np.tile(X_check_sim[1,],(len_poss_na,1))
            m_for_diff_2 = np.tile(X_check_sim[2,],(len_poss_na,1))
            m_for_diff_3 = np.tile(X_check_sim[3,],(len_poss_na,1))
            m_for_diff_4 = np.tile(X_check_sim[4,],(len_poss_na,1))

            differences_0 = np.sum(np.bitwise_xor(m_for_diff_0, X_poss_na), axis=1)
            differences_1 = np.sum(np.bitwise_xor(m_for_diff_1, X_poss_na), axis=1)
            differences_2 = np.sum(np.bitwise_xor(m_for_diff_2, X_poss_na), axis=1)
            differences_3 = np.sum(np.bitwise_xor(m_for_diff_3, X_poss_na), axis=1)
            differences_4 = np.sum(np.bitwise_xor(m_for_diff_4, X_poss_na), axis=1)

            non_zero_0 = np.sum(differences_0 != 0)
            differences_0 = np.where(differences_0 == 0,  np.nan, differences_0)
            corr_ord_0 = np.argsort(differences_0)[0:non_zero_0+1]

            non_zero_1 = np.sum(differences_1 != 0)
            differences_1 = np.where(differences_1 == 0,  np.nan, differences_1)
            corr_ord_1 = np.argsort(differences_1)[0:non_zero_1+1]

            non_zero_2 = np.sum(differences_2 != 0)
            differences_2 = np.where(differences_2 == 0,  np.nan, differences_2)
            corr_ord_2 = np.argsort(differences_2)[0:non_zero_2+1]

            non_zero_3 = np.sum(differences_3 != 0)
            differences_3 = np.where(differences_3 == 0,  np.nan, differences_3)
            corr_ord_3 = np.argsort(differences_3)[0:non_zero_3+1]

            non_zero_4 = np.sum(differences_4 != 0)
            differences_4 = np.where(differences_4 == 0,  np.nan, differences_4)
            corr_ord_4 = np.argsort(differences_4)[0:non_zero_4+1]

            to_add = list()
            for (a,b,c,d,e) in zip(list(corr_ord_0),list(corr_ord_1),list(corr_ord_2),list(corr_ord_3),list(corr_ord_4)):
                if not(a in to_add): to_add.append(a)
                if not(b in to_add): to_add.append(b)
                if not(c in to_add): to_add.append(c)
                if not(d in to_add): to_add.append(d)
                if not(e in to_add): to_add.append(e)
                if len(to_add) > missing_neg:
                    break # stop if we have enough similar genes

            # add list_genomes_to_add to the X_na
            for i in to_add:
                negative_examples_subsample.append(possible_neg[i])



    return positive_examples_subsample, negative_examples_subsample

# function that train the classifier for one node ==============================
def train_classifier(positive_examples,negative_examples,all_classifiers,alignment, node, penalty_v, solver_v):
    # check that we have at least 1 example for each class:
    if len(negative_examples) == 0:
        # when the node is the only child, then there are no negative examples
        logging.info('      Warning: no negative examples for "%s', node)
        return "no_negative_examples"
    if len(positive_examples) == 0:
        # There should be positive examples
        logging.info('      Error: no positive examples for "%s', node)
        return "ERROR_no_positive_examples"

    # select the genes from the pandas dataframe
    X = alignment.loc[ negative_examples + positive_examples , : ].to_numpy()
    train_labels = ["no"]*len(negative_examples)+["yes"]*len(positive_examples)
    # NOTE: we put first the negative class (0) because then the classifier will
    #       use this order. And when we will use only the coefficients, it will
    #       give the probability prediction of the secodn class

    y = np.asarray(train_labels)
    # train classifier
    clf = LogisticRegression(random_state=0, penalty = penalty_v, solver=solver_v)
    clf.fit(X, y)
    return clf


# train node and call the same function on all the children ====================
def train_node_iteratively(node, sibilings, all_classifiers, alignment, full_taxonomy, penalty_v, solver_v):
    # call the function on all the children
    # but only if they are not the last level
    if not(full_taxonomy.is_last_node(node)):
        children_of_node = full_taxonomy.find_children_node(node)
        for child in children_of_node:
            sibilings_child = list(children_of_node)
            sibilings_child.remove(child)
            train_node_iteratively(child, sibilings_child, all_classifiers, alignment, full_taxonomy, penalty_v, solver_v)

    # find genomes to use and to which class they belong to,
    # we need positive and negative examples
    logging.info('   TRAIN:"%s":Find genes', node)
    positive_examples, negative_examples = find_training_genes(node, sibilings, full_taxonomy, alignment)
    logging.info('      SEL_GENES:"%s": %s positive, %s negative', node,
                 str(len(positive_examples)),str(len(negative_examples)))

    # train the classifier
    logging.info('         TRAIN:"%s":Train classifier', node)
    all_classifiers[node] = train_classifier(positive_examples,negative_examples,
                                             all_classifiers, alignment, node, penalty_v, solver_v)


# function to train all classifiers ============================================
# this function will create a classifier for each node in the taxonomy
# Input:
#  - the aligned sequences as a pandas data frame
#  - the taxonomy (global variable)
# Output:
#  - a dictionary, where the keys are the node names and the values are a lasso
#                  classifier object
def train_all_classifiers(alignment, full_taxonomy, penalty_v, solver_v):
    all_classifiers = dict()
    children_of_root = full_taxonomy.find_children_node(full_taxonomy.get_root())
    for node in children_of_root:
        sibilings = list(children_of_root)
        sibilings.remove(node)
        train_node_iteratively(node, sibilings, all_classifiers, alignment, full_taxonomy, penalty_v, solver_v)
    return(all_classifiers)









#===============================================================================
#              FUNCTIONS TO LEARN THE FUNCTION FOR THE TAX. LEVEL
#===============================================================================
def predict_iter(test_seq, training_tax, classifiers_train, tax, perc, arrived_so_far):
    if training_tax.is_last_node(arrived_so_far):
        return
    max_perc = 0
    max_perc_taxa = ""
    # if there is only one child:
    if len(training_tax.find_children_node(arrived_so_far)) == 1:
        max_perc = 2 # if there are no sibilings I put 2, it will be replaced after
        max_perc_taxa = training_tax.find_children_node(arrived_so_far)[0]
    # if there are no children
    if len(training_tax.find_children_node(arrived_so_far)) < 1:
        sys.stderr.write("Error: no child\n")
    # if there is more than one child
    if len(training_tax.find_children_node(arrived_so_far)) > 1:
        for n in training_tax.find_children_node(arrived_so_far):
            clf = classifiers_train[n]
            res = str(clf.predict(test_seq)) # either "yes" or "no"
            predictions = clf.predict_proba(test_seq) # two predictions
            if res == "['no']":
                predicted_proba = np.amin(predictions) # if it predicts no, then the probability that we select is the smaller one
            else:
                predicted_proba = np.amax(predictions)
            # check if the prediction is higher
            if predicted_proba > max_perc:
                max_perc = predicted_proba
                max_perc_taxa = n

    tax.append(max_perc_taxa)
    perc.append(max_perc)
    predict_iter(test_seq, training_tax, classifiers_train, tax, perc, max_perc_taxa)

def predict_one_gene(test_seq, training_tax, classifiers_train):
    tax = list()
    perc = list()
    # we arrived at the root, and now we classify from there
    predict_iter(test_seq, training_tax, classifiers_train, tax, perc, training_tax.get_root())
    # we change the predictions that came from having only one sibiling --------
    if perc[0] == 2:
        perc[0] = 1
    for i in range(len(perc)):
        if perc[i] == 2:
            perc[i] = perc[i-1]

    return tax, perc

def predict(test_al, training_tax, classifiers_train):
    res = list()
    for i in test_al.index.values:
        r = list()
        r.append(i)
        predictions, percentages = predict_one_gene([test_al.loc[ i , : ].to_numpy()], training_tax, classifiers_train)
        r.append(predictions)
        r.append(percentages)
        res.append(r)
    return(res)

def learn_function_one_level(level_to_learn, alignment, full_taxonomy, penalty_v, solver_v):
    logging.info('  TEST:"%s" taxonomic level', str(level_to_learn))
    # 1. Identify which clades we want to remove (test set) and which to keep
    #    (training set)
    this_level_clades = full_taxonomy.find_node_level(level_to_learn)
    # we use 33% of the clades for testing
    perc_test_set = 0.33 # this cannot be higher than 0.5.
    test_set = set()
    training_set = set()
    for c in this_level_clades:
        # find how many to use for the test set:
        aval_clades = set(this_level_clades[c])
        n_test = round(len(aval_clades) * perc_test_set)
        if len(aval_clades) == 2:
            n_test = 0
        # add to test set
        for i in range(n_test):
            test_set.add(aval_clades.pop())
        # the remaining clades in aval_clades go to the trainig set
        training_set.update(aval_clades)
    logging.info('  TEST:"%s" level:test_set (%s):%s', str(level_to_learn),str(len(test_set)),str(test_set))
    logging.info('  TEST:"%s" level:trai_set (%s):%s', str(level_to_learn),str(len(training_set)),str(training_set))

    # 2. Create new taxonomy and alignment file & train the classifiers
    training_tax = full_taxonomy.copy()
    removed_genes = training_tax.remove_clades(list(test_set))
    training_al = alignment.loc[ training_tax.find_gene_ids(training_tax.get_root()) , : ]
    classifiers_train = train_all_classifiers(training_al, training_tax, penalty_v, solver_v)

    # 3. Classify the test set
    test_al = alignment.loc[ removed_genes , : ]
    pr = predict(test_al, training_tax, classifiers_train)
    for g in pr:
        # g is:
        # ["geneB",["A","B","D","species8"],[0.99,0.96,0.96,0.07]]
        correct_tax = full_taxonomy.extract_full_tax_from_gene(g[0])
        g.append(correct_tax)
        g.append(level_to_learn)

    return pr
    # return:
    #    GENE_ID         PREDICTED             PROB_PREDICTED        CORRECT        REMOVED_LEVEL
    # [["geneA",["A","B","C","species2"],[0.98,0.97,0.23,0.02],["A","B","Y","speciesX"],2]
    #  ["geneB",["A","B","D","species8"],[0.99,0.96,0.10,0.07],["A","B","U","speciesZ"],2]
    # .....                                                                               ]

def learn_function_genes_level(level_to_learn, alignment, full_taxonomy, penalty_v, solver_v):
    logging.info('  TEST:"%s" taxonomic level', str(level_to_learn))
    # 1. Identify which clades we want to remove (test set) and which to keep
    #    (training set)
    this_level_clades = full_taxonomy.get_last_level_to_genes() # now there are genes
    # we use 33% of the genes for testing
    perc_test_set = 0.33 # this cannot be higher than 0.5.
    test_set = set()
    training_set = set()
    for c in this_level_clades:
        # find how many to use for the test set:
        aval_clades = set(this_level_clades[c])
        n_test = round(len(aval_clades) * perc_test_set)
        # add to test set
        for i in range(n_test):
            test_set.add(aval_clades.pop())
        # the remaining clades in aval_clades go to the trainig set
        training_set.update(aval_clades)
    logging.info('  TEST:"%s" level:test_set (%s):%s', str(level_to_learn),str(len(test_set)),str(test_set))
    logging.info('  TEST:"%s" level:trai_set (%s):%s', str(level_to_learn),str(len(training_set)),str(training_set))

    # 2. Create new taxonomy and alignment file & train the classifiers
    training_tax = full_taxonomy.copy()
    training_tax.remove_genes(list(test_set))
    training_al = alignment.loc[ training_set , : ]
    classifiers_train = train_all_classifiers(training_al, training_tax, penalty_v, solver_v)

    # 3. Classify the test set
    test_al = alignment.loc[ test_set , : ]
    pr = predict(test_al, training_tax, classifiers_train)
    for g in pr:
        # g is:
        # ["geneB",["A","B","D","species8"],[0.99,0.96,0.96,0.07]]
        correct_tax = full_taxonomy.extract_full_tax_from_gene(g[0])
        g.append(correct_tax)
        g.append(level_to_learn)

    return pr

def estimate_function(all_calc_functions):
    # The all_calc_functions looks like:
    #    GENE_ID         PREDICTED             PROB_PREDICTED        CORRECT        REMOVED_LEVEL
    # [["geneA",["A","B","C","species2"],[0.98,0.97,0.23,0.02],["A","B","Y","speciesX"],2]
    #  ["geneB",["D","E","F","species8"],[0.99,0.96,0.95,0.07],["D","E","F","speciesZ"],3]
    #  ["geneC",["G","H","I","species9"],[0.99,0.96,0.95,0.94],["G","H","I","species9"],4]
    #   .....                                                                             ]
    # the REMOVED_LEVEL count starting from 0 (in the first example, 2 means that
    # "C" is removed)
    # when we have 4 (which is outside of the taxonomy levels, {0,1,2,3}), it
    # refers to the fact that we removed the genes

    # we remove duplicates with the same predicted probability -----------------
    all_uniq = dict()
    for line in all_calc_functions:
        v = ""
        for j in line[2]:
            v = v+str(j)
        all_uniq[v] = line
    logging.info('   LEARN_FUNCTION:Number of lines: %s (before removing duplicates: %s)',
                      str(len(all_uniq)),str(len(all_calc_functions)))
    # change all_calc_functions
    all_calc_functions = list()
    for j in all_uniq:
        all_calc_functions.append(all_uniq[j])

    # we find what is the correct value for the prediction level ---------------
    correct_level = list()
    for line in all_calc_functions:
        corr_level_this = -1
        cont = 0
        for p,c in zip(line[1],line[3]):
            cont = cont + 1
            if p == c:
                corr_level_this = cont-1 # we select to what level to predict
        correct_level.append(corr_level_this)
    # now in correct_level there is to which level to predict to. Example:
    # "A","B","C","species2"
    # with corr_level_this = 0, we should assign "A"
    # with corr_level_this = 2, we should assign "A","B","C"
    # with corr_level_this = -1, we should assign "" (no taxonomy)

    # check how many lines there are per correct level -------------------------
    for l in set(correct_level):
        cont = 0
        for j in correct_level:
            if j == l:
                cont = cont + 1
        logging.info('   LEARN_FUNCTION:Number of lines: level %s: %s',
                          str(l),str(cont))


    # we train the classifiers -------------------------------------------------
    all_classifiers = dict()
    for l in set(correct_level):
        # we create the feature matrix
        # NOTE: we always need the negative class to be first
        correct_order_lines = list()
        correct_order_labels = list()
        cont = 0
        for i in range(len(all_calc_functions)):
            if correct_level[i] != l:
                correct_order_lines.append(all_calc_functions[i][2])
                correct_order_labels.append(0)
            cont = cont + 1
        cont = 0
        for i in range(len(all_calc_functions)):
            if correct_level[i] == l:
                correct_order_lines.append(all_calc_functions[i][2])
                correct_order_labels.append(1)
            cont = cont + 1

        X = np.array([np.array(xi) for xi in correct_order_lines])
        y = np.asarray(correct_order_labels)
        # train classifier
        clf = LogisticRegression(random_state=0, penalty = "none", solver='saga',max_iter = 5000)
        clf.fit(X, y)
        all_classifiers[str(l)] = clf

    return all_classifiers

# create taxonomy selection function ===========================================
# This function define a function that is able to identify to which taxonomic
# level a new gene should be assigned to.
def learn_taxonomy_selection_function(alignment, full_taxonomy, save_cross_val_data, penalty_v, solver_v):
    # find number of levels
    n_levels = full_taxonomy.get_n_levels()

    # do the cross validation for each level
    all_calc_functions = list()
    for i in range(n_levels):
        all_calc_functions = all_calc_functions + learn_function_one_level(i, alignment, full_taxonomy, penalty_v, solver_v)
    # do the cross val. for the last level (using the genes)
    all_calc_functions = all_calc_functions + learn_function_genes_level(n_levels, alignment, full_taxonomy, penalty_v, solver_v)

    # save all_calc_functions if necessary -------------------------------------
    if not (save_cross_val_data is None):
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        outfile.write("gene\tpredicted\tprob\tground_truth\tremoved_level\n")
        os.chmod(outfile.name, 0o644)
        for vals in all_calc_functions:
            to_print = vals[0] + "\t" + "/".join(vals[1]) + "\t" # "geneB",["D","E","F","species8"]
            to_print = to_print + "/".join(str(x) for x in vals[2]) + "\t" # [0.99,0.96,0.95,0.07]
            to_print = to_print + "/".join(vals[3]) + "\t" # ["G","H","I","species9"]
            to_print = to_print + str(vals[4]) # removed level
            outfile.write(to_print+"\n")
        # save
        try:
            outfile.flush()
            os.fsync(outfile.fileno())
            outfile.close()
        except:
            sys.stderr.write("[E::main] Error: failed to save the cross validation results\n")
        try:
            #os.rename(outfile.name,output) # atomic operation
            shutil.move(outfile.name,save_cross_val_data) #It is not atomic if the files are on different filsystems.
        except:
            sys.stderr.write("[E::main] Error: failed to save the cross validation results\n")
            sys.stderr.write("[E::main] you can find the file here:\n"+outfile.name+"\n")
            sys.exit(1)

    # estimate the function ----------------------------------------------------
    f = estimate_function(all_calc_functions)
    return f



#===============================================================================
#                     FUNCTIONS TO SAVE TO A DATABASE
#===============================================================================
def save_to_file(classifiers, full_taxonomy, tax_function, use_cmalign, hmm_file_path, tool_version, output, protein_fasta_input):
    # where to save the file
    f = h5py.File(output, "w")
    string_dt = h5py.special_dtype(vlen=str)

    # zero: tool version -------------------------------------------------------
    f.create_dataset('tool_version',data=np.array([str(tool_version)],"S100"),dtype=string_dt)
    # and type of database
    f.create_dataset('db_type',data=np.array(["single_gene"],"S100"),dtype=string_dt)
    # was the alignment done at the protein level?
    if not(protein_fasta_input is None):
        f.create_dataset('align_protein',data=np.array([True]),dtype=bool)
    else:
        f.create_dataset('align_protein',data=np.array([False]),dtype=bool)

    # first we save the hmm file -----------------------------------------------
    line = ""
    o = open(hmm_file_path,"r")
    for i in o:
        line = line + i
    o.close()
    f.create_dataset('hmm_file',data=np.array([line],"S"+str(len(line)+100)),dtype=string_dt, compression="gzip")

    # second, save the use_cmalign info ----------------------------------------
    f.create_dataset('use_cmalign',data=np.array([use_cmalign]),dtype=bool)

    # third, we save the taxonomy ---------------------------------------------
    f.create_group("taxonomy")
    for i in full_taxonomy.child_nodes:
        f.create_dataset("taxonomy/"+i, data=np.array(list(full_taxonomy.child_nodes[i]),"S10000"),dtype=string_dt, compression="gzip")

    # fourth, the taxonomy function --------------------------------------------
    f.create_group("tax_function")
    for c in tax_function:
        # we append the intercept at the head (will have position 0)
        vals = np.append(tax_function[c].intercept_, tax_function[c].coef_)
        f.create_dataset("tax_function/"+str(c), data=vals, dtype=np.float64, compression="gzip")

    # fifth, save the classifiers ----------------------------------------------
    f.create_group("classifiers")
    for c in classifiers:
        if classifiers[c] != "no_negative_examples":
            vals = np.append(classifiers[c].intercept_, classifiers[c].coef_)
            f.create_dataset("classifiers/"+c, data=vals, dtype=np.float64, compression="gzip", compression_opts=8)
        else:
            # in this case, it always predict 1, we save it as an array of
            # with the string "no_negative_examples"
            f.create_dataset("classifiers/"+c,data=np.array(["no_negative_examples"],"S40"),dtype=string_dt, compression="gzip")

    # close hdm5 file ----------------------------------------------------------
    f.flush()
    f.close()








#===============================================================================
#                                      MAIN
#===============================================================================

def create_db(aligned_seq_file, tax_file, verbose, output, use_cmalign, hmm_file_path, save_cross_val_data, tool_version, protein_fasta_input, penalty_v, solver_v):
    # set log file
    filename_log = os.path.realpath(output)+'.log'
    logging.basicConfig(filename=filename_log,
                        filemode='w',
                        level=logging.INFO,
                        format='[%(asctime)s] %(message)s')
    logging.info('TIME:start')

    # 1. load the taxonomy into the tree (global variable)
    logging.info('MAIN:Load taxonomy')
    full_taxonomy = Taxonomy(tax_file)
    full_taxonomy.load_from_file()
    logging.info('TIME:Finish load taxonomy')

    # 2. load the alignment into a pandas dataframe
    logging.info('MAIN:Load alignment')
    alignment = load_alignment_from_file(aligned_seq_file)
    logging.info('TIME:Finish load alignment')

    # 3. check that the taxonomy and the alignment are consistent
    logging.info('MAIN:Check taxonomy and alignment')
    check_taxonomy_alignment_consistency(alignment, full_taxonomy)
    logging.info('TIME:Finish check-up')

    # 4. build a classifier for each node
    logging.info('MAIN:Train all classifiers')
    classifiers = train_all_classifiers(alignment, full_taxonomy, penalty_v, solver_v)
    logging.info('TIME:Finish train all classifiers')

    # 5. learn the function to identify the correct taxonomy level
    logging.info('MAIN:Learn taxonomy selection function')
    tax_function = learn_taxonomy_selection_function(alignment, full_taxonomy, save_cross_val_data, penalty_v, solver_v)
    logging.info('TIME:Finish learn taxonomy selection function')

    # 6. save the result
    logging.info('MAIN:Save to file')
    save_to_file(classifiers, full_taxonomy, tax_function, use_cmalign, hmm_file_path, tool_version, output, protein_fasta_input)
    logging.info('TIME:Finish save to file')

    logging.info('MAIN:Finished')
