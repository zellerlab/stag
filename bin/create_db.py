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
# Global variables for taxonomy tree
child_nodes = dict()
tree_root = "tree_root"
child_nodes[tree_root] = set()
last_level_to_genes = dict()
all_gene_ids = list()

#===============================================================================
#                                 FUNCTIONS
#===============================================================================


# function that find the leaves (hence gene ids) of a given node ===============
# Input:
#  - the node
# Output:
#  - a list of genes
def find_leaves_recoursive(node, all_leaves):
    if node in last_level_to_genes:
        all_leaves.extend(last_level_to_genes[node])
    else:
        for c in child_nodes[node]:
            find_leaves_recoursive(c, all_leaves)

def find_leaves(node):
    all_leaves = list()
    find_leaves_recoursive(node, all_leaves)
    return (all_leaves)

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

# function to load the taxonomy ================================================
def load_taxonomy(file_name):
    o = open(file_name,"r")

    number_of_taxonomic_levels = len(o.readline().rstrip().split("\t"))
    o.seek(0)

    for line in o:
        # expected line: gene1\tBacteria\tFirmicutes\t...
        vals = line.rstrip().split("\t")
        if number_of_taxonomic_levels != len(vals):
            sys.stderr.write("Error: taxonomy file does not have the same number of taxonomic levels in:\n")
            sys.stderr.write("  "+line+"\n")
            sys.exit(1)

        # we enter the first level, to the root:
        child_nodes[tree_root].add(vals[1])
        if not(vals[1] in child_nodes):
            child_nodes[vals[1]] = set()
        # we enter all remaining levels
        for i in range(2,len(vals)-1):
            # first we enter that this is a child
            child_nodes[vals[i-1]].add(vals[i])
            # and second, we create a node if there is not already
            if not(vals[i] in child_nodes):
                child_nodes[vals[i]] = set()
        # We add the last level
        child_nodes[vals[-2]].add(vals[-1])
        # Finally we add from the last level to the genes ids
        if not(vals[-1] in last_level_to_genes):
            last_level_to_genes[vals[-1]] = set()
        last_level_to_genes[vals[-1]].add(vals[0])
        # and we add it to the list of gene ids
        all_gene_ids.append(vals[0])

    all_gene_ids.sort() # we sort the list, so that search should be faster
    o.close()

# function that finds positive and negative examples ===========================
def find_training_genes(node, sibilings):
    positive_examples = find_leaves(node)
    negative_examples = list()
    if len(sibilings) > 0:
        for s in sibilings:
            negative_examples = negative_examples + find_leaves(s)
    return positive_examples, negative_examples

# function that train the classifier for one node ==============================
def train_classifier(positive_examples,negative_examples,all_classifiers,alignment):
    # select the genes from the pandas dataframe
    
    # train classifier
    clf = LogisticRegression(random_state=0, penalty = "l1", solver='liblinear')
    clf.fit(X, y)
    return clf


# train node and call the same function on all the children ====================
def train_node_iteratively(node, sibilings, all_classifiers, alignment):
    # call the function on all the children
    # but only if they are not the last level
    if not(node in last_level_to_genes):
        for child in child_nodes[node]:
            sibilings_child = list(child_nodes[node])
            sibilings_child.remove(child)
            train_node_iteratively(child, sibilings_child, all_classifiers, alignment)

    # find genomes to use and to which class they belong to,
    # we need positive and negative examples
    logging.info('   TRAIN:"%s":Find genes', node)
    positive_examples, negative_examples = find_training_genes(node, sibilings)
    logging.info('      SEL_GENES:"%s": %s positive, %s negative', node,
                 str(len(positive_examples)),str(len(negative_examples)))

    # train the classifier
    logging.info('   TRAIN:"%s":Train classifier', node)
    all_classifiers[node] = train_classifier(positive_examples,negative_examples,
                                             all_classifiers, alignment)


# function to train all classifiers ============================================
# this function will create a classifier for each node in the taxonomy
# Input:
#  - the aligned sequences as a pandas data frame
#  - the taxonomy (global variable)
# Output:
#  - a dictionary, where the keys are the node names and the values are a lasso
#                  classifier object
def train_all_classifiers(alignment):
    all_classifiers = dict()
    for node in child_nodes[tree_root]:
        sibilings = list(child_nodes[tree_root])
        sibilings.remove(node)
        train_node_iteratively(node, sibilings, all_classifiers, alignment)
    return(all_classifiers)

# main function ================================================================
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
    load_taxonomy(tax_file)
    logging.info('MAIN:Finish load taxonomy')

    # 2. load the alignment into a pandas dataframe
    logging.info('MAIN:Load alignment')
    alignment = load_alignment_from_file(aligned_seq_file)
    logging.info('MAIN:Finish load alignment')

    # 3. build a classifier for each node
    logging.info('MAIN:Train all classifiers')
    classifiers = train_all_classifiers(alignment)
    logging.info('MAIN:Finish train all classifiers')
