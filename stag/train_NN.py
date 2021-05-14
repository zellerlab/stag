"""
Functions to train the nearest neighnour classifiers
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import numpy as np
import pandas as pd
import sys
import metric_learn


# ------------------------------------------------------------------------------
# load taxonomy.
# we need a different way to load the taxonomy here:
def load_tax_line(tax_file, ALI):
    selected_seq = list(ALI.index.values)
    res = dict()
    o = open(tax_file,"r")
    for line in o:
        vals = line.rstrip().split("\t")
        if vals[0] in selected_seq:
            res[vals[0]] = vals[1].split(";")
    o.close()
    return res

# ------------------------------------------------------------------------------
# function to estimate the weights
def estimate_weights(ALI, tax, sel_level):
    # find all families (or other level), we test with sel_level = 4
    all_clades = dict()
    for seq in tax:
        if not tax[seq][sel_level] in all_clades:
            all_clades[tax[seq][sel_level]] = list()
        all_clades[tax[seq][sel_level]].append(seq)
    # now we do the analysis by clade
    all_LMNN = dict()
    all_transformed = dict()
    for clade in all_clades:
        print(clade)
        # for the X --------------------------------------------
        # we subselect the training
        this_ALI = ALI.loc[all_clades[clade],:]
        # transform from pandas to numpy (for the features)
        X = 1*this_ALI.to_numpy()
        #X = X + 0.1

        # for the y  -------------------------------------------
        # we need to create the species ground thruth
        rownames = this_ALI.index.values
        list_species = list()
        for i in rownames:
            list_species.append(tax[i][-1])
        # we need numbers
        species_2_num = dict()
        cont = 1
        for i in list(set(list_species)):
            species_2_num[i] = cont
            cont = cont + 1
        y = list()
        for i in list_species:
            y.append(species_2_num[i])
        y = np.array(y)

        # we learn the transformation --------------------------
        lmnn = metric_learn.LMNN(k=1, learn_rate=1e-2,regularization = 0.4)
        # fit the data
        lmnn.fit(X, y)
        #TODO: check that it converges, you have to parse the output printed
        #      with verbose

        # transform our input space ----------------------------
        X_lmnn = lmnn.transform(X)
        # create a panda object with the transformed space and the correct
        # rownames
        X_lmnn_PD = pd.DataFrame(X_lmnn, index=rownames)

        # add to dict ------------------------------------------
        all_LMNN[clade] = lmnn
        all_transformed[clade] = X_lmnn_PD

    return all_LMNN, all_transformed


#===============================================================================
#                                      MAIN
#===============================================================================
def train_NN_classifiers(alignment, tax_file, NN_start_level):
    # 0. load the taxonomy
    tax = load_tax_line(tax_file, alignment)
    # 1. we find the transformations and we transform the original space
    all_LMNN, all_transformed = estimate_weights(alignment, tax, NN_start_level)
    return "dummy"
