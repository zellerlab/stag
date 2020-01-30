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

#===============================================================================
# Global variables for taxonomy tree
child_nodes = dict()
tree_root = "tree_root"
child_nodes[tree_root] = set()
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
    if node in all_gene_ids:
        all_leaves.append(node)
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
        for i in range(2,len(vals)):
            # first we enter that this is a child
            child_nodes[vals[i-1]].add(vals[i])
            # and second, we create a node if there is not already
            if not(vals[i] in child_nodes):
                child_nodes[vals[i]] = set()
        # Finally, we add the gene id
        child_nodes[vals[-1]].add(vals[0])
        all_gene_ids.append(vals[0])

    all_gene_ids.sort() # we sort the list, so that search should be faster
    o.close()

# main function ================================================================
def create_db(aligned_seq_file, tax_file, verbose, output):
    # 1. load the taxonomy into the tree (global variable)
    load_taxonomy(tax_file)
    # 2. load the alignment into a pandas dataframe
    alignment = load_alignment_from_file(aligned_seq_file)
    # 3. build a classifier for each node
    classifiers = train_all_classifiers(alignment)
