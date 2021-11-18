"""
Functions to train the nearest neighnour classifiers
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import numpy as np
import pandas as pd
import sys
import metric_learn
import statistics
import multiprocessing as mp

from sklearn.metrics import precision_recall_curve

logging = "global_logging"
verbose = "global_verbose"

#===============================================================================
#                                      UTIL
#===============================================================================
# ------------------------------------------------------------------------------
# load taxonomy.
# we need a different way to load the taxonomy here:
def load_tax_line(tax_file, ALI):
    species_to_tax = dict()
    selected_seq = list(ALI.index.values)
    res = dict()
    o = open(tax_file,"r")
    for line in o:
        vals = line.rstrip().split("\t")
        if vals[0] in selected_seq:
            res[vals[0]] = vals[1].split(";")
            species_to_tax[vals[1].split(";")[-1]] = vals[1].split(";")
    o.close()
    return res, species_to_tax

#===============================================================================
#                          TRANSFORM THE SPACE
#===============================================================================
# given a 1-hot ancoding numpy array, it removes columns that are all the same
def remove_invariant_columns(X_full):
    num_rows, num_cols = X_full.shape
    # if there is only one sequence
    if num_rows == 1:
        col_to_keep = list(range(num_cols))
        logging.info('    keeping %s columns (%s percent)', str(len(col_to_keep)),str(len(col_to_keep)/num_cols))
        return X_full, col_to_keep

    # normal case: more than one seqeunce
    # colsum
    colsum = X_full.sum(axis=0)
    # which are not all zeros
    pos_all_non_0 = np.where(colsum != 0)[0].tolist()
    # which are not all ones
    pos_all_non_1 = np.where(colsum != num_rows)[0].tolist()
    # we intersect the two sets
    col_to_keep = list(set(pos_all_non_0) & set(pos_all_non_1))
    col_to_keep.sort()

    logging.info('    keeping %s columns (%s percent)', str(len(col_to_keep)),str(len(col_to_keep)/num_cols))

    return X_full[:,col_to_keep], col_to_keep



# MAIN function to estimate the weights
def estimate_weights(ALI, tax, sel_level, procs=1):

    def get_x_columns(ALI, clade):
        # we subselect the training
        this_ALI = ALI.loc[all_clades[clade],:]
        # transform from pandas to numpy (for the features)
        X_full = 1 * this_ALI.to_numpy()
        # find rownames
        rownames = this_ALI.index.values

        return X_full,rownames


    # find all families (or other level), we test with sel_level = 4
    all_clades = dict()
    for seq in tax:
        all_clades.setdefault(tax[seq][sel_level], list()).append(seq)
    # now we do the analysis by clade
    all_LMNN = dict()
    all_transformed = dict()
    all_sel_positions = dict()

    clades_to_compute = set()
    for clade in all_clades:
        logging.info('    TRAIN_NN_3: Clade: %s', clade)
        X_full,rownames = get_x_columns(ALI, clade)
        X, sel_positions = remove_invariant_columns(X_full)
        # which columns were selected for this clade
        all_sel_positions[clade] = sel_positions


        # for the y  -------------------------------------------
        # we need to create the species ground thruth
        list_species = [tax[row][-1] for row in rownames]
        # we need numbers
        species_2_num = {
            species: i for i, species in enumerate(set(list_species), start=1)
        }
        y = np.array([species_2_num[species] for species in list_species])

        if len(set(y)) == 1 or len(y) <= 5:
            if verbose > 5: sys.stderr.write("------------------- "+clade+": NOT ENOUGH DATA\n")
            all_LMNN[clade] = "NOT ENOUGH DATA"
            message = f'Not enough data ({len(y)})' if len(y) <= 5 else 'Only one species'
            logging.info(f'     TRAIN_NN_4: {message}')
        else:
            logging.info('     TRAIN_NN_4: Fit the data')
            if verbose > 4: sys.stderr.write("----- "+clade+"("+str(len(y))+"): will train\n")
            if procs == 1:
                all_LMNN[clade], all_transformed[clade] = estimate_weights_for_clade(X, y,rownames)
            else:
                clades_to_compute.add((clade, tuple(y)))

    if clades_to_compute:
        with mp.Pool(processes=procs) as pool:
            results = [
                pool.apply_async(
                    estimate_weights_for_clade,
                    args=(get_x_columns(ALI, clade)[0][:, all_sel_positions[clade]], np.array(y),get_x_columns(ALI, clade)[0])
                )
                for clade, y in clades_to_compute
            ]

            for res in results:
                all_LMNN[clade], all_transformed[clade] = res.get()


    return all_LMNN, all_transformed, all_sel_positions

def estimate_weights_for_clade(X, y, rownames):
    # we learn the transformation --------------------------
    if verbose > 4: sys.stderr.write("---------- ("+str(len(y))+"): start fit\n")
    lmnn = metric_learn.LMNN(k=1, learn_rate=1e-2, regularization=0.4)
    lmnn.fit(X, y)
    if verbose > 4: sys.stderr.write("---------- ("+str(len(y))+"): finish fit\n")

    #TODO: check that it converges, you have to parse the output printed
    #      with verbose

    logging.info('     TRAIN_NN_4: Transform the data')
    X_lmnn = lmnn.transform(X)

    # create a panda object with the transformed space and the correct
    # rownames
    X_lmnn_PD = pd.DataFrame(X_lmnn, index=rownames)

    return lmnn, X_lmnn_PD




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
def find_thresholds(all_transformed, tax):
    # prepare the result
    threshold_clades = dict()
    list_centroids_all = dict() # will contain pandas array

    for clade in all_transformed:
        logging.info('   TRAIN_NN_2: Clade: %s', clade)
        logging.info('    TRAIN_NN_3: Find centroids')
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
        list_centroids_all[clade] = centroids_this_pd

        # calc all distances for this clade --------------
        logging.info('    TRAIN_NN_3: Calculate distances')
        dist_all_vs_centroids = calc_all_dist_to_centroids(centroids_this,all_transformed[clade],tax)

        # find the thresholds -------------------
        logging.info('    TRAIN_NN_3: Find thresholds')
        threshold_clades[clade] = find_thresholds_from_dist(dist_all_vs_centroids)

    return threshold_clades, list_centroids_all


#===============================================================================
#                                      MAIN
#===============================================================================
def train_NN_classifiers(alignment, tax_file, NN_start_level,logging_, verbose_, procs=1):
    # set logging
    global logging
    logging = logging_
    # set verbose
    global verbose
    verbose = verbose_

    # 0. load the taxonomy
    logging.info('  TRAIN_NN_1: load tax')
    if verbose > 4: sys.stderr.write("-- Load taxonomy\n")
    tax, species_to_tax = load_tax_line(tax_file, alignment)

    # 1. we calculate the transformations and we transform the original space
    logging.info('  TRAIN_NN_1: calculate LMNN')
    if verbose > 4: sys.stderr.write("-- Calculate LMNN\n")
    all_LMNN, all_transformed, all_sel_positions = estimate_weights(alignment, tax, NN_start_level, procs=procs)

    # 2. find centroids and find the threshold distances
    logging.info('  TRAIN_NN_1: find centroids and thresholds')
    if verbose > 4: sys.stderr.write("-- Find centroids and thresholds\n")
    thresholds_NN, centroid_seq = find_thresholds(all_transformed, tax)

    return all_LMNN, thresholds_NN, centroid_seq, species_to_tax, all_sel_positions