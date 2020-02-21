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

from sklearn.linear_model import LogisticRegression
import numpy as np
import sys
import random
import time
import pandas as pd
import logging
import os
from sklearn.linear_model import LogisticRegression
from sklearn.datasets import make_classification

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

    # load taxonomy from the defined file --------------------------------------
    def load_from_file(self):
        o = open(self.file_name,"r")

        self.number_of_taxonomic_levels = len(o.readline().rstrip().split("\t"))
        o.seek(0)

        for line in o:
            # expected line: gene1\tBacteria\tFirmicutes\t...
            vals = line.rstrip().split("\t")
            if self.number_of_taxonomic_levels != len(vals):
                sys.stderr.write("Error: taxonomy file does not have the same number of taxonomic levels in:\n")
                sys.stderr.write("  "+line+"\n")
                sys.exit(1)

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
        return temp

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
        set_removed_genes = set()
        for n in node_list:
            remove_clade_iter(n, set_removed_genes)
            # now need to remove from the higher level in child_nodes
            for i in self.child_nodes:
                self.child_nodes[i].discard(n) # discard does not raise a KeyError.
        return list(set_removed_genes)

    def remove_clade_iter(self, node, set_removed_genes):
        if node in self.last_level_to_genes:
            # we arrived at the end of the tree, we remove the genes, but first:
            # add to the set of removed genes
            set_removed_genes = set_removed_genes.union(self.last_level_to_genes[node])
            # remove the genes from the gene list
            self.all_gene_ids = [e for e in self.all_gene_ids if e not in self.last_level_to_genes[node]]
            # and, finally, remove the node from the last_level_to_genes dict
            self.last_level_to_genes.pop(node,None)

        for n in self.child_nodes[node]:
            remove_clade_iter(n, set_removed_genes)

        # remove from child_nodes
        self.child_nodes.pop(node,None)




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
    alignment = alignment.astype('bool') # apparently you cannot load directly
                                         # bool if the rownames are not bool
    return alignment

# function to check that taxonomy and alignment are consistent =================
# 1. all genes in the alignment should be in the taxonomy
def check_taxonomy_alignment_consistency(alignment, full_taxonomy):
    ff = "sdf"






#===============================================================================
#                   FUNCTIONS TO TRAIN THE CLASSIFIERS
#===============================================================================

# function that finds positive and negative examples ===========================
def find_training_genes(node, sibilings, full_taxonomy):
    positive_examples = full_taxonomy.find_gene_ids(node)
    negative_examples = list()
    if len(sibilings) > 0:
        for s in sibilings:
            negative_examples = negative_examples + full_taxonomy.find_gene_ids(s)
    return positive_examples, negative_examples

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
    X = alignment.loc[ positive_examples + negative_examples , : ].to_numpy()
    train_labels = ["yes"]*len(positive_examples)+["no"]*len(negative_examples)
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
    positive_examples, negative_examples = find_training_genes(node, sibilings, full_taxonomy)
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

def learn_function_one_level(level_to_learn):
    return ["dummy"]
def learn_function_genes_level():
    return ["dummy"]

# create taxonomy selection function ===========================================
# This function define a function that is able to identify to which taxonomic
# level a new gene should be assigned to.
def learn_taxonomy_selection_function(alignment, full_taxonomy):
    # find number of levels
    n_levels = full_taxonomy.get_n_levels()

    # do the cross validation for each level
    all_calc_functions = list()
    for i in range(n_levels):
        all_calc_functions = all_calc_functions + learn_function_one_level(i)
    # do the cross val. for the last level (using the genes)
    all_calc_functions = all_calc_functions + learn_function_genes_level()
    # estimate the function









#===============================================================================
#                                      MAIN
#===============================================================================

def create_db(aligned_seq_file, tax_file, verbose, output):
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
    # TODO
    check_taxonomy_alignment_consistency(alignment, full_taxonomy)

    # 4. build a classifier for each node
    logging.info('MAIN:Train all classifiers')
    classifiers = train_all_classifiers(alignment, full_taxonomy)
    logging.info('MAIN:Finish train all classifiers')

    # 5. learn the function to identify the correct taxonomy level
    learn_taxonomy_selection_function(alignment, full_taxonomy)

    # 6. save the result
