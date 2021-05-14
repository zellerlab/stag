"""
Functions to train the nearest neighnour classifiers
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import numpy as np
import pandas as pd
import sys
import metric_learn

#===============================================================================
#                                      UTIL
#===============================================================================
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

#===============================================================================
#                          TRANSFORM THE SPACE
#===============================================================================
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
        all_transformed[clade] = X_lmnn_PD # it's a pandas object

    return all_LMNN, all_transformed



#===============================================================================
#                     FIND CENTROIDS AND EVAL. THRESHOLDS
#===============================================================================
# we want to select only one sequence per species (centroid) and then calculate
# centorids vs all to identify the threshold for genus and species (if for
# example the NN_start_level is family)

# ------------------------------------------------------------------------------
# Calculate perc identity between two sequences
#  - ALI is the pandas DF
#  - pos1 and pos2 are part of the rownames of ALI
def dist_vectors(ALI,pos1,pos2):
    seq1 = ALI.loc[pos1,].to_numpy()
    seq2 = ALI.loc[pos2,].to_numpy()
    # euclidean distance
    dist = np.linalg.norm(seq1-seq2)
    return dist



# ------------------------------------------------------------------------------
# FIND CENTROIDS
# It will return a dictionary where the keys are the centroid sequences and the
# value is the species
def find_centroids(ALI,tax):
    # first find all species and their map to the seq id
    all_species = dict()
    for seq in ALI.index.values:
        if not tax[seq][-1] in all_species:
            all_species[tax[seq][-1]] = list()
        all_species[tax[seq][-1]].append(seq)
    # now we go by species and find the centroid
    result = dict()
    for species in all_species:
        # if there is only one
        if len(all_species[species]) == 1:
            result[all_species[species][0]] = species
            continue
        # if there are 2, we choose randomly
        if len(all_species[species]) == 2:
            result[all_species[species][0]] = species
            continue
        # if there are at least 3, then we calculate the centroid
        all_id_dist = dict()
        for id in all_species[species]:
            all_id_dist[id] = 0
        for i in range(len(all_species[species])-1):
            for j in range(i+1,len(all_species[species])):
                id_i = all_species[species][i]
                id_j = all_species[species][j]
                # find vector and distance
                dist_this = dist_vectors(ALI,id_i,id_j)
                # add distance
                all_id_dist[id_i] = all_id_dist[id_i] + dist_this
                all_id_dist[id_j] = all_id_dist[id_j] + dist_this
        # now we need to identify the gene with the minimum distance (i.e. centroid)
        sel_min = ""
        val_min = float('inf')
        for i in all_id_dist:
            if all_id_dist[i] < val_min:
                val_min = all_id_dist[i]
                sel_min = i
        # we add it
        result[species] = sel_min

    # the result is a ditionary where the values are:
    # "E.coli" -> "gene1"
    # "P.copri" -> "gene236 "

    # return the result
    return result

def find_thresholds(all_transformed, tax):
    # find centroids per species
    for clade in all_transformed:
        # find centroids per species
        centroids_this = find_centroids(all_transformed[clade],tax)

    return "dummy1", "dummy2"

#===============================================================================
#                                      MAIN
#===============================================================================
def train_NN_classifiers(alignment, tax_file, NN_start_level):
    # 0. load the taxonomy
    tax = load_tax_line(tax_file, alignment)
    # 1. we calculate the transformations and we transform the original space
    all_LMNN, all_transformed = estimate_weights(alignment, tax, NN_start_level)

    # 2. find centroids and find the threshold distances
    thresholds_NN, centroid_seq = find_thresholds(all_transformed, tax)

    return all_LMNN, thresholds_NN, centroid_seq
