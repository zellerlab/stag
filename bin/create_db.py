"""
Scripts that creates the database of classifiers
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

# Input:
#  - one multiple sequence alignment (MSA) per marker gene. The MSA is obtained
#    from the function htc align, like:
#       >gene1\t0001010001010000100101000...
#       >gene2\t0000110001010100100101001...
#  - a taxonomy file that describes the taxonomy of the genes:
#       gene1\tBacteria\tFirmicutes\t...
#
# Output:
#  - a database file (hdf5) that can be used by htc classify

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

#===============================================================================
#                          CLASS FOR THE TAXONOMY
#===============================================================================
class Taxonomy:
    # create class -------------------------------------------------------------
    def __init__(self, file_name):
        self.file_name = file_name
        self.child_nodes = dict()
        self.tree_root = "tree_root"
        self.child_nodes[self.tree_root] = set()
        self.last_level_to_genes = dict()
        self.all_gene_ids = list()
        self.number_of_taxonomic_levels = 0
        self.annotation_per_gene = dict()

    # load taxonomy from the defined file --------------------------------------
    def load_from_file(self):
        o = open(self.file_name,"r")

        self.number_of_taxonomic_levels = len(o.readline().rstrip().split("\t")) - 1
        o.seek(0)

        for line in o:
            # expected line: gene1\tBacteria\tFirmicutes\t...
            vals = line.rstrip().replace("/","-").split("\t")
            if self.number_of_taxonomic_levels != len(vals)-1 :
                sys.stderr.write("Error: taxonomy file does not have the same number of taxonomic levels in:\n")
                sys.stderr.write("  "+line+"\n")
                sys.exit(1)

            # we add the annotation_per_gene:
            self.annotation_per_gene[vals[0]] = list(vals[1:])

            # we enter the first level, to the root:
            self.child_nodes[self.tree_root].add(vals[1])
            if not(vals[1] in self.child_nodes):
                self.child_nodes[vals[1]] = set()
            # we enter all remaining levels
            for i in range(2,len(vals)-1):
                # first we enter that this is a child
                self.child_nodes[vals[i-1]].add(vals[i])
                # and second, we create a node if there is not already
                if not(vals[i] in self.child_nodes):
                    self.child_nodes[vals[i]] = set()
            # We add the last level
            self.child_nodes[vals[-2]].add(vals[-1])
            # Finally we add from the last level to the genes ids
            if not(vals[-1] in self.last_level_to_genes):
                self.last_level_to_genes[vals[-1]] = set()
            self.last_level_to_genes[vals[-1]].add(vals[0])
            # and we add it to the list of gene ids
            self.all_gene_ids.append(vals[0])

        self.all_gene_ids.sort() # we sort the list, so that search should be faster
        o.close()

    # make a copy of this taxonomy ---------------------------------------------
    def copy(self):
        temp = Taxonomy(self.file_name)
        temp.child_nodes = dict()
        temp.tree_root = "tree_root"
        for i in self.child_nodes:
            temp.child_nodes[i] = set(self.child_nodes[i])
        temp.last_level_to_genes = dict()
        for i in self.last_level_to_genes:
            temp.last_level_to_genes[i] = set(self.last_level_to_genes[i])
        temp.all_gene_ids = list(self.all_gene_ids)
        temp.number_of_taxonomic_levels = self.number_of_taxonomic_levels
        for i in self.annotation_per_gene:
            temp.annotation_per_gene[i] = list(self.annotation_per_gene[i])
        return temp

    # return the classification of one gene
    def extract_full_tax_from_gene(self, gene_id):
        return self.annotation_per_gene[gene_id]

    # return number of levels --------------------------------------------------
    def get_n_levels(self):
        return self.number_of_taxonomic_levels

    # return the root id -------------------------------------------------------
    def get_root(self):
        return self.tree_root

    # find children of a node --------------------------------------------------
    def find_children_node(self, node):
        if node in self.child_nodes:
            return list(self.child_nodes[node])
        else:
            return None
    # return the last level to genes -------------------------------------------
    def get_last_level_to_genes(self):
        last_level_to_genes_temp = dict()
        for i in self.last_level_to_genes:
            last_level_to_genes_temp[i] = set(self.last_level_to_genes[i])
        return last_level_to_genes_temp

    # check if it is the last node before the genes ----------------------------
    def is_last_node(self, node):
        if node in self.last_level_to_genes:
            return True
        else:
            return False
    # find all genes under a given node ----------------------------------------
    # return a list of all genes
    def find_gene_ids(self, node):
        all_leaves = list()
        self.find_leaves_recoursive(node, all_leaves)
        return all_leaves
    def find_leaves_recoursive(self, node, all_leaves):
        if node in self.last_level_to_genes:
            all_leaves.extend(self.last_level_to_genes[node])
        else:
            for c in self.child_nodes[node]:
                self.find_leaves_recoursive(c, all_leaves)

    # function to remove nodes (and genes underneath), given a list of nodes ---
    # it returns the gene ids that were removed
    def remove_clades(self, node_list):
        # remove all clades under
        list_removed_genes = list()
        for n in node_list:
            self.remove_clade_iter(n, list_removed_genes)
        # remove all clades on top
        for n in node_list:
            # now need to remove from the higher level in child_nodes
            for i in self.child_nodes:
                self.child_nodes[i].discard(n) # discard does not raise a KeyError.
        # if it was the only child, then we should remove also at higher level
        self.remove_unused_branches()
        return list(list_removed_genes)

    def remove_clade_iter(self, node, list_removed_genes):
        if node in self.last_level_to_genes:
            # we arrived at the end of the tree, we remove the genes, but first:
            # add to the set of removed genes
            list_removed_genes.extend(self.last_level_to_genes[node])
            # remove the genes from the gene list
            self.all_gene_ids = [e for e in self.all_gene_ids if e not in self.last_level_to_genes[node]]
            # and, finally, remove the node from the last_level_to_genes dict
            self.last_level_to_genes.pop(node,None)
        else:
            try:
                for n in self.child_nodes[node]:
                    self.remove_clade_iter(n, list_removed_genes)
                # remove from child_nodes
                self.child_nodes.pop(node,None)
            except KeyError as e:
                sys.stderr.write("WARNING: key not present when removing a clade ["+str(e)+"]\n\n")


    def remove_unused_branches(self):
        removed_any = False # this becomes True if we remove any node from
                            # child_nodes, in which case we re-run remove_unused_branches
        list_to_remove = list()
        for i in self.child_nodes:
            if len(self.child_nodes[i]) == 0:
                removed_any = True
                list_to_remove.append(i)
                # remove from taxonomy at higher level
                for j in self.child_nodes:
                    self.child_nodes[j].discard(i)
        # remove nodes that are empty from child_nodes
        for n in list_to_remove:
            self.child_nodes.pop(n,None)
        # call remove_unused_branches again if necessary
        if removed_any:
            self.remove_unused_branches()

    # function to remove genes from a list -------------------------------------
    def remove_genes(self, gene_list):
        # remove the genes from the gene list
        self.all_gene_ids = [e for e in self.all_gene_ids if e not in gene_list]
        # remove the genes from last_level_to_genes
        for g in gene_list:
            for node in self.last_level_to_genes:
                self.last_level_to_genes[node].discard(g)
        # Check if all the genes from one clade are removed, and hence we should
        # remove that clade
        list_to_remove = list()
        for node in self.last_level_to_genes:
            if len(self.last_level_to_genes[node]) == 0:
                list_to_remove.append(node)
        self.remove_clades(list_to_remove)

    # function that returns all nodes at one level, ----------------------------
    # as a dictionary of the parent nodes
    def find_tax_level_iter(self, current_node, current_level, result):
        if current_node in self.child_nodes:
            for n in self.child_nodes[current_node]:
                result[n] = current_level+1
                self.find_tax_level_iter(n, current_level+1, result)
    def find_node_level(self, tax_level_find):
        # find tax level for each node
        tax_level = dict()
        tax_level[self.tree_root] = 0
        self.find_tax_level_iter(self.tree_root,0,tax_level)
        # select only the one from the correct level
        res = dict()
        for n in self.child_nodes:
            if tax_level[n] == tax_level_find:
                res[n] = set(self.child_nodes[n])
        return res

    # print the values in the taxonomy class -----------------------------------
    def __str__(self):
        to_print = "NODES:\n"
        for i in self.child_nodes:
            to_print = to_print + "   (N):" + i + ": " + str(self.child_nodes[i]) + "\n"
        to_print = to_print + "\nGENES:\n"
        for i in self.last_level_to_genes:
            to_print = to_print + "   (G):" + i + ": " + str(self.last_level_to_genes[i]) + "\n"
        to_print = to_print + "\nLIST GENES:\n" + str(self.all_gene_ids) + "\n"
        to_print = to_print + "\nN LEVELS: " + str(self.number_of_taxonomic_levels) + "\n"
        return to_print




#===============================================================================
#                   FUNCTIONS TO LOAD AND CHECK THE ALIGNMENT
#===============================================================================

# function to load an alignment produced by the "align" option =================
# Input:
#  - a file created by "align"
# Output:
#  - a panda object
# as a note, numpy.loadtxt is way slower than pandas read.csv
# It works also on .gz files
def load_alignment_from_file(file_name):
    alignment = pd.read_csv(file_name,delimiter='\t',index_col = 0, header=None)
    logging.info('   LOAD_AL: Number of genes: %s', str(len(list(alignment.index.values))))
    alignment = alignment.astype('bool') # apparently you cannot load directly
                                         # bool if the rownames are not bool
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
        possible_neg = set(alignment.index.values).difference(set(positive_examples + negative_examples))
        if len(possible_neg) == 0: # if it is possible to add negatives
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
def train_classifier(positive_examples,negative_examples,all_classifiers,alignment, node):
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
    clf = LogisticRegression(random_state=0, penalty = "l1", solver='liblinear')
    clf.fit(X, y)
    return clf


# train node and call the same function on all the children ====================
def train_node_iteratively(node, sibilings, all_classifiers, alignment, full_taxonomy):
    # call the function on all the children
    # but only if they are not the last level
    if not(full_taxonomy.is_last_node(node)):
        children_of_node = full_taxonomy.find_children_node(node)
        for child in children_of_node:
            sibilings_child = list(children_of_node)
            sibilings_child.remove(child)
            train_node_iteratively(child, sibilings_child, all_classifiers, alignment, full_taxonomy)

    # find genomes to use and to which class they belong to,
    # we need positive and negative examples
    logging.info('   TRAIN:"%s":Find genes', node)
    positive_examples, negative_examples = find_training_genes(node, sibilings, full_taxonomy, alignment)
    logging.info('      SEL_GENES:"%s": %s positive, %s negative', node,
                 str(len(positive_examples)),str(len(negative_examples)))

    # train the classifier
    logging.info('         TRAIN:"%s":Train classifier', node)
    all_classifiers[node] = train_classifier(positive_examples,negative_examples,
                                             all_classifiers, alignment, node)


# function to train all classifiers ============================================
# this function will create a classifier for each node in the taxonomy
# Input:
#  - the aligned sequences as a pandas data frame
#  - the taxonomy (global variable)
# Output:
#  - a dictionary, where the keys are the node names and the values are a lasso
#                  classifier object
def train_all_classifiers(alignment, full_taxonomy):
    all_classifiers = dict()
    children_of_root = full_taxonomy.find_children_node(full_taxonomy.get_root())
    for node in children_of_root:
        sibilings = list(children_of_root)
        sibilings.remove(node)
        train_node_iteratively(node, sibilings, all_classifiers, alignment, full_taxonomy)
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
        max_perc = 1
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

def learn_function_one_level(level_to_learn, alignment, full_taxonomy):
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
    classifiers_train = train_all_classifiers(training_al, training_tax)

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

def learn_function_genes_level(level_to_learn, alignment, full_taxonomy):
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
    classifiers_train = train_all_classifiers(training_al, training_tax)

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

def estimate_function(all_calc_functions, n_levels):
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

    # some removed clades have more genes than others, we remove some genes from
    # the most abundant clades
    n_levels = len(all_calc_functions[0][1])
    # This will contain all the empirical distributions:
    all_values = dict()

    # parse each line
    line_pass_filter = dict()
    cont = -1
    for line in all_calc_functions:
        cont = cont + 1
        rem_lev = line[4]
        rem_clade = line[3][min( rem_lev,len(line[3])-1 )] # we use min, because if we removed a gene, then we put the species ID
        # inizialize all levels to zero
        if rem_lev not in all_values:
            all_values[rem_lev] = dict()

        # check that till the known clade we have the correct classification
        predicted = list(line[1])
        correct =  list(line[3])
        for i in range(rem_lev,n_levels):
            predicted[i] = "NA"
            correct[i] = "NA"
        s_predicted = "/".join(predicted)
        s_correct = "/".join(correct)
        if s_predicted != s_correct:
            s_perc = '/'.join(str(v) for v in line[2])
            logging.info('   LEARN_FUNCTION:WARNING:MISCLASSIFICATION: (gene):%s, (predicted): %s, (correct): %s, (perc):%s',
                              line[0],s_predicted, "/".join(line[3]), s_perc)
        # if they have the same taxonomy
        if s_predicted == s_correct:
            if rem_clade not in all_values[rem_lev]:
                all_values[rem_lev][rem_clade] = 0
            all_values[rem_lev][rem_clade] = all_values[rem_lev][rem_clade] + 1
            # we keep track of which lines passed the filter
            if rem_lev not in line_pass_filter:
                line_pass_filter[rem_lev] = list()
            line_pass_filter[rem_lev].append(cont)

    # now we have
    # all_values[removed_level][removed_clade] -> number of genes
    # all_values          [INT]         [INT]                [INT]

    # we find max number of genes per each level
    selected_lines = list()
    selected_lines_rem_level = list()
    for l in all_values:
        n_genes_per_clade = list()
        for i in all_values[l]:
            n_genes_per_clade.append(all_values[l][i])
        # we set the number of maximum genes per clade
        n_max = np.percentile(np.array(n_genes_per_clade), 75) # return 75th percentile
        logging.info('   LEARN_FUNCTION:NUMBER_PER_CLADE: level:%s, number:%s',
                          str(l),str(n_max))
        # we set to 0 the count of genes per clade
        for i in all_values[l]:
            all_values[l][i] = 0
        # we shuffle the list and start to add
        random.shuffle(line_pass_filter[l])
        for i in line_pass_filter[l]:
            line = all_calc_functions[i]
            rem_lev = line[4]
            rem_clade = line[3][min( rem_lev,len(line[3])-1 )] # we use min, because if we removed a gene, then we put the species ID
            if all_values[l][rem_clade] < n_max:
                all_values[l][rem_clade] = all_values[l][rem_clade] + 1
                selected_lines.append(line[2])
                selected_lines_rem_level.append(rem_lev)

    # now we use the selected lines to train a classifier for each level
    all_classifiers = dict()
    levels = list(line_pass_filter.keys())
    for l in levels:
        # we create the feature matrix
        # NOTE: we always need the negative class to be first
        correct_order_lines = list()
        correct_order_labels = list()
        cont = 0
        for i in selected_lines_rem_level:
            if i != l:
                correct_order_lines.append(selected_lines[cont])
                correct_order_labels.append(0)
            cont = cont + 1
        cont = 0
        for i in selected_lines_rem_level:
            if i == l:
                correct_order_lines.append(selected_lines[cont])
                correct_order_labels.append(1)
            cont = cont + 1

        X = np.array([np.array(xi) for xi in correct_order_lines])
        y = np.asarray(correct_order_labels)
        # train classifier
        clf = LogisticRegression(random_state=0, penalty = "l1", solver='liblinear',max_iter = 5000)
        clf.fit(X, y)
        all_classifiers[l] = clf

    return all_classifiers


# create taxonomy selection function ===========================================
# This function define a function that is able to identify to which taxonomic
# level a new gene should be assigned to.
def learn_taxonomy_selection_function(alignment, full_taxonomy):
    # find number of levels
    n_levels = full_taxonomy.get_n_levels()

    # do the cross validation for each level
    all_calc_functions = list()
    for i in range(n_levels):
        all_calc_functions = all_calc_functions + learn_function_one_level(i, alignment, full_taxonomy)
    # do the cross val. for the last level (using the genes)
    all_calc_functions = all_calc_functions + learn_function_genes_level(n_levels, alignment, full_taxonomy)

    # estimate the function
    f = estimate_function(all_calc_functions, n_levels)
    return f




#===============================================================================
#                     FUNCTIONS TO SAVE TO A DATABASE
#===============================================================================
def save_to_file(classifiers, full_taxonomy, tax_function, use_cmalign, hmm_file_path, output):
    # where to save the file
    f = h5py.File(output, "w")

    # first we save the hmm file -----------------------------------------------
    line = ""
    o = open(hmm_file_path,"r")
    for i in o:
        line = line + i
    o.close()
    string_dt = h5py.special_dtype(vlen=str)
    f.create_dataset('hmm_file',data=np.array([line],"S"+str(len(line)+100)),dtype=string_dt, compression="gzip")

    # second, save the use_cmalign info ----------------------------------------
    f.create_dataset('use_cmalign',data=np.array([use_cmalign]),dtype=bool)

    # third, we save the taxonomy ---------------------------------------------
    f.create_group("taxonomy")
    for i in full_taxonomy.child_nodes:
        f.create_dataset("taxonomy/"+i, data=np.array(list(full_taxonomy.child_nodes[i]),"S1000"),dtype=string_dt, compression="gzip")

    # fourth, the taxonomy function --------------------------------------------
    f.create_group("tax_function")
    for c in tax_function:
        f.create_dataset("tax_function/"+str(c), data=tax_function[c].coef_,dtype=np.float64, compression="gzip")

    # fifth, save the classifiers ----------------------------------------------
    f.create_group("classifiers")
    for c in classifiers:
        if classifiers[c] != "no_negative_examples":
            f.create_dataset("classifiers/"+c, data=classifiers[c].coef_,dtype=np.float64, compression="gzip", compression_opts=8)
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

def create_db(aligned_seq_file, tax_file, verbose, output, use_cmalign, hmm_file_path):
    # set log file
    filename_log = os.path.realpath(output)+'.log'
    logging.basicConfig(filename=filename_log,
                        filemode='w',
                        level=logging.INFO,
                        format='[%(asctime)s] %(message)s')
    logging.info('MAIN:start')

    # 1. load the taxonomy into the tree (global variable)
    logging.info('MAIN:Load taxonomy')
    full_taxonomy = Taxonomy(tax_file)
    full_taxonomy.load_from_file()
    logging.info('MAIN:Finish load taxonomy')

    # 2. load the alignment into a pandas dataframe
    logging.info('MAIN:Load alignment')
    alignment = load_alignment_from_file(aligned_seq_file)
    logging.info('MAIN:Finish load alignment')

    # 3. check that the taxonomy and the alignment are consistent
    logging.info('MAIN:Check taxonomy and alignment')
    check_taxonomy_alignment_consistency(alignment, full_taxonomy)
    logging.info('MAIN:Finish check-up')

    # 4. build a classifier for each node
    logging.info('MAIN:Train all classifiers')
    classifiers = train_all_classifiers(alignment, full_taxonomy)
    logging.info('MAIN:Finish train all classifiers')

    # 5. learn the function to identify the correct taxonomy level
    logging.info('MAIN:Learn taxonomy selection function')
    tax_function = learn_taxonomy_selection_function(alignment, full_taxonomy)
    logging.info('MAIN:Finish learn taxonomy selection function')

    # 6. save the result
    logging.info('MAIN:Save to file')
    save_to_file(classifiers, full_taxonomy, tax_function, use_cmalign, hmm_file_path, output)
    logging.info('MAIN:Finish save to file')
