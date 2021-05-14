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
    # "P.copri" -> "gene236"

    # return the result
    return result


# ------------------------------------------------------------------------------
# calculate the distance of all vs centroids
def calc_all_dist_to_centroids(centroids_this,ALI,tax):
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
            d = dist_vectors(ALI,sel_centroid,gene)
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
        res[i] = 0
    return res


# ------------------------------------------------------------------------------
def find_thresholds(all_transformed, tax):
    # find centroids per species
    threshold_clades = dict()

    list_centroids_all = list()

    for clade in all_transformed:
        # find centroids per species
        gene_centroids = list()
        species_centroids = list()
        # run find_centroids
        centroids_this = find_centroids(all_transformed[clade],tax)
        # add them to the result
        for c in centroids_this:
            gene_centroids.append(centroids_this[c])
            species_centroids.append(c)
        # create a panda array
        centroids_this_pd = all_transformed[clade].loc[gene_centroids,].to_numpy()
        # where the rownames are the species
        centroids_this_pd = pd.DataFrame(centroids_this_pd, index=species_centroids)
        list_centroids_all.append(centroids_this_pd)

        # calc all distances for this clade --------------
        dist_all_vs_centroids = calc_all_dist_to_centroids(centroids_this,all_transformed[clade],tax)

        # find the thresholds -------------------
        threshold_clades[clade] = find_thresholds_from_dist(dist_all_vs_centroids)

    centroids_all = pd.concat(list_centroids_all)
    return threshold_clades, centroids_all


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
