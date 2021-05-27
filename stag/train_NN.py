"""
Functions to train the nearest neighnour classifiers
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import numpy as np
import pandas as pd
import sys
import metric_learn
import statistics
from sklearn.metrics import precision_recall_curve

logging = "global_logging"

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
    all_ali = dict()
    for clade in all_clades:
        logging.info('    TRAIN_NN_3: Clade: %s', clade)
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

        if len(y) > 5:
            # we learn the transformation --------------------------
            lmnn = metric_learn.LMNN(k=1, learn_rate=1e-2,regularization = 0.4)
            # fit the data
            lmnn.fit(X, y)
            #TODO: check that it converges, you have to parse the output printed
            #      with verbose

            # transform our input space ----------------------------
            #X_lmnn = lmnn.transform(X)
            # create a panda object with the transformed space and the correct
            # rownames
            #X_lmnn_PD = pd.DataFrame(X_lmnn, index=rownames)

            # add to dict ------------------------------------------
            all_LMNN[clade] = lmnn
        else:
            all_LMNN[clade] = "NOT ENOUGH DATA"

        # we need the original alignment for calculating distances
        all_ali[clade] = pd.DataFrame(X, index=rownames)

    return all_LMNN, all_ali



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
def dist_vectors(ALI,pos1,pos2,LMNN_model):
    xi = ALI.loc[pos1,].to_numpy()
    xj = ALI.loc[pos2,].to_numpy()
    if LMNN_model != "NOT ENOUGH DATA":
        # Get mahalanobis matrix
        m = LMNN_model.get_mahalanobis_matrix()
        # distance
        dist = np.sqrt( (( xi-xj ) .dot(m) ).dot(xi-xj))
    else:
        # if there was not enough data to calculate the LMNN, then we calculate the
        # euclidean distance on the untransformed alignments
        seq1 = ALI.loc[pos1,].to_numpy()
        seq2 = ALI.loc[pos2,].to_numpy()
        # euclidean distance
        dist = np.linalg.norm(seq1-seq2)
    return dist



# ------------------------------------------------------------------------------
# FIND CENTROIDS
# It will return a dictionary where the keys are the centroid sequences and the
# value is the species
def find_centroids(ALI,tax,LMNN_model_this):
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
            result[species] = all_species[species][0]
            continue
        # if there are 2, we choose randomly
        if len(all_species[species]) == 2:
            result[species] = all_species[species][0]
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
                dist_this = dist_vectors(ALI,id_i,id_j,LMNN_model_this)
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
    # "P.copri" -> "gene236"

    # return the result
    return result


# ------------------------------------------------------------------------------
# calculate the distance of all vs centroids
def calc_all_dist_to_centroids(centroids_this,ALI,tax,LMNN_this):
    # here we are already inside one clade
    n_tax_level = len(list(tax.values())[0])-1
    # find all possible genes
    all_genes = list(ALI.index.values)

    # create dict for the result
    all_dist = dict()

    for species in centroids_this:
        sel_centroid = centroids_this[species]
        sel_centroid_tax = tax[sel_centroid]
        # now we calculate all possible distances from this centroid
        all_genes_this = all_genes
        all_genes_this.remove(sel_centroid)
        for gene in all_genes_this:
            # distance
            d = dist_vectors(ALI,sel_centroid,gene,LMNN_this)
            # check taxonomy
            for i in range(n_tax_level,-1,-1):
                if sel_centroid_tax[i] == tax[gene][i]:
                    sel_level = i
                    break
            # now in sel_level there is the level to predict
            # i.e. if sel_level == 1 it means they have the same phylum
            if not sel_level in all_dist:
                all_dist[sel_level] = list()
            all_dist[sel_level].append(d)

    # we check for log
    to_print = ""
    for i in all_dist:
        to_print = to_print + str(i) + "["+str(len(all_dist[i]))+"] "
    logging.info('     TRAIN_NN_4: Training vals: %s', to_print)

    return all_dist



# ------------------------------------------------------------------------------
# find the thresholds given the distances
# distances is like (if we use `-L 4`):
# {4:[6.3,9.5,11.4,22.4,3.6],
#  5:[2.3,2.4,2.0,2.1],
#  6:[0.3,0.6,0.5,1.2,1.0]}
# where 5 means genus
def find_thresholds_from_dist(distances):
    res = dict()
    # we check for each level:
    for i in distances:
        if i != min(list(distances.keys())): # the lowest (example 4) was already assigned
            if i in distances and i-1 in distances:
                negative_vals = distances[i-1]
                positive_vals = distances[i]
                y_true = np.array(([0]*len(negative_vals)) + ([1]*len(positive_vals)))
                y_scores = np.array(negative_vals + positive_vals)
                precision, recall, thresholds = precision_recall_curve(y_true, y_scores)

                # TODO does precision_recall_curve check enough values?!
                logging.info(f'     TRAIN_NN_4: Number of thresholds: {i}: {len(thresholds)}')

                # we calcualte the F1 score and take the maximum
                n_valid_thresholds = 0
                maxF1 = 0
                sel_threshold = 0
                for pr, re, th in zip(precision, recall, thresholds):
                    if pr !=0 and re != 0:
                        n_valid_thresholds = n_valid_thresholds + 1
                        F1_this = 2*( (pr*re)/(pr+re) )
                        if F1_this > maxF1:
                            maxF1 = F1_this
                            sel_threshold = th
                logging.info(f'     TRAIN_NN_4: Number of real thresholds: {i}: {n_valid_thresholds}')

                # save to res
                res[i] = sel_threshold
    return res


# ------------------------------------------------------------------------------
def find_thresholds(all_ali, tax, all_LMNN):
    # prepare the result
    threshold_clades = dict()
    list_centroids_all = dict() # will contain pandas array

    for clade in all_ali:
        logging.info('   TRAIN_NN_2: Clade: %s', clade)
        logging.info('    TRAIN_NN_3: Find centroids')
        # find centroids per species
        gene_centroids = list()
        species_centroids = list()
        # run find_centroids
        centroids_this = find_centroids(all_ali[clade],tax,all_LMNN[clade])
        # add them to the result
        for c in centroids_this:
            gene_centroids.append(centroids_this[c])
            species_centroids.append(c)
        # create a panda array
        centroids_this_pd = all_ali[clade].loc[gene_centroids,].to_numpy()
        # where the rownames are the species
        centroids_this_pd = pd.DataFrame(centroids_this_pd, index=species_centroids)
        list_centroids_all[clade] = centroids_this_pd

        # calc all distances for this clade --------------
        logging.info('    TRAIN_NN_3: Calculate distances')
        dist_all_vs_centroids = calc_all_dist_to_centroids(centroids_this,all_ali[clade],tax,all_LMNN[clade])

        # find the thresholds -------------------
        logging.info('    TRAIN_NN_3: Find thresholds')
        threshold_clades[clade] = find_thresholds_from_dist(dist_all_vs_centroids)

    return threshold_clades, list_centroids_all


#===============================================================================
#                                      MAIN
#===============================================================================
def train_NN_classifiers(alignment, tax_file, NN_start_level,logging_):
    # set logging
    global logging
    logging = logging_

    # 0. load the taxonomy
    logging.info('  TRAIN_NN_1: load tax')
    tax = load_tax_line(tax_file, alignment)

    # 1. we calculate the transformations and divide the alignment
    logging.info('  TRAIN_NN_1: calculate LMNN')
    all_LMNN, all_ali = estimate_weights(alignment, tax, NN_start_level)

    # 2. find centroids and find the threshold distances
    logging.info('  TRAIN_NN_1: find centroids and thresholds')
    thresholds_NN, centroid_seq = find_thresholds(all_ali, tax, all_LMNN)

    return all_LMNN, thresholds_NN, centroid_seq
