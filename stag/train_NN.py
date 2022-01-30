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

from sklearn.metrics import f1_score

from stag import UTIL_log

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

# ------------------------------------------------------------------------------
# save distances to file
# dist_all_vs_centroids is like (if we use `-L 4`):
# {4:[6.3,9.5,11.4,22.4,3.6],
#  5:[2.3,2.4,2.0,2.1],
#  6:[0.3,0.6,0.5,1.2,1.0]}
# where 5 means genus
def save_distances_to_file(distances, clade, intermediate_dist_for_NN):
    with open(intermediate_dist_for_NN, "a") as o:
        for dis in distances:
            o.write(clade+"\t"+str(dis))
            for i in distances[dis]:
                o.write("\t"+str(i))
            o.write("\n")


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
def estimate_weights(ALI, tax, sel_level, min_training_data_lmnn, base_save_features, procs=1):

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

        if len(set(y)) == 1 or len(y) <= min_training_data_lmnn:
            if verbose > 5: sys.stderr.write("------------------- "+clade+": NOT ENOUGH DATA\n")
            all_LMNN[clade] = "NOT ENOUGH DATA"
            # we add the untransformed data
            all_transformed[clade] = pd.DataFrame(X, index=rownames)
            message = f'Not enough data ({len(y)})' if len(y) <= min_training_data_lmnn else 'Only one species'
            logging.info(f'     TRAIN_NN_4: {message}')
        else:
            logging.info('     TRAIN_NN_4: Fit the data')
            if verbose > 4: sys.stderr.write("----- "+clade+"("+str(len(y))+"): will train\n")
            if procs == 1:
                all_LMNN[clade], all_transformed[clade],unused = estimate_weights_for_clade(X, y, rownames, clade, verbose, base_save_features)
            else:
                clades_to_compute.add((clade, tuple(y)))

    if verbose > 4:
        UTIL_log.print_log("  Finished first pass through all clades\n")
        logging.info('       TRAIN_NN_5: Finished first pass')
    if clades_to_compute:
        if verbose > 4:
            UTIL_log.print_log("  Enter multiprocessing\n")
            logging.info('       TRAIN_NN_5: Enter multiprocessing')
        with mp.get_context("spawn").Pool(processes=procs) as pool:
            results = [
                pool.apply_async(
                    estimate_weights_for_clade,
                    args=(get_x_columns(ALI, clade)[0][:, all_sel_positions[clade]], np.array(y), get_x_columns(ALI, clade)[1], clade, verbose, base_save_features)
                )
                for clade, y in clades_to_compute
            ]

            for res in results:
                lmnn, transformed, cla = res.get()
                all_LMNN[cla] = lmnn
                all_transformed[cla] = transformed


    if verbose > 4:
        UTIL_log.print_log("  Finished train of all NN\n")
        logging.info(' Finished train of all NN')
    return all_LMNN, all_transformed, all_sel_positions

def estimate_weights_for_clade(X, y, rownames, clade, verb, base_save_features=None):
    # we learn the transformation --------------------------
    if verb > 4:
        UTIL_log.print_log("---------- ("+str(len(y))+"): start fit\n")
    lmnn = metric_learn.LMNN(k=1, learn_rate=1e-2, regularization=0.4)
    lmnn.fit(X, y)
    if verb > 4:
        UTIL_log.print_log("---------- ("+str(len(y))+"): finish fit\n")

    #TODO: check that it converges, you have to parse the output printed
    #      with verbose

    X_lmnn = lmnn.transform(X)

    # create a panda object with the transformed space and the correct
    # rownames
    if verb > 5: sys.stderr.write("--------------- ("+str(len(y))+"): dimX_lmnn="+str(X_lmnn.shape[0])+"x"+str(X_lmnn.shape[1])+"; dim_rownames="+str(len(rownames))+"\n")
    X_lmnn_PD = pd.DataFrame(X_lmnn, index=rownames)

    # if we have to save the features, we do it now:
    if(base_save_features is not None):
        X_lmnn_PD.to_csv(base_save_features+"/"+clade+'_transformed.gz', compression= 'gzip')
        pd.DataFrame(X, index=rownames).to_csv(base_save_features+"/"+clade+'_original.gz', compression= 'gzip')

    return lmnn, X_lmnn_PD, clade




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
# help function for the find_thresholds_from_dist
def find_thresholds_this(negative_vals,positive_vals):
    # find all possible thresholds, so like given these distances:
    # [1,3,4,8,10,20] the intermediate values are:
    # [2,3.5,6,9,15] which are the possible thresholds
    all_vals = negative_vals+positive_vals
    all_vals.sort()
    #
    res = list()
    for i in range(len(all_vals)-1):
        res.append((all_vals[i+1]+all_vals[i])/2)
    return res

def measureF1(negative_vals,positive_vals,thresholds):
    res_f1 = list()
    for t in thresholds:
        y_true = np.array(([0]*len(negative_vals)) + ([1]*len(positive_vals)))
        y_pred = list()
        for i in negative_vals:
            if i < t:
                # we measure distances, hence this would be a FP
                y_pred.append(1)
            else:
                y_pred.append(0)
        for i in positive_vals:
            if i < t:
                # we measure distances, hence this would be a TP
                y_pred.append(1)
            else:
                y_pred.append(0)
        y_pred = np.array(y_pred)

        if len(set(y_true) - set(y_pred)) > 0:
            # it means y_pred contains only zeros, we set the F1 score to zero
            # directly, so that I don't get a warning from np.f1_score
            res_f1.append(0)
        else:
            res_f1.append(f1_score(y_true, y_pred,average='weighted'))
    return res_f1

def find_max(thresholds_in,thresholdsF1_in,negative_vals,positive_vals):
    max_val = 0
    max_pos = -1
    thresh_1s = list()
    for i in range(len(thresholds_in)):
        if thresholdsF1_in[i] > max_val:
            max_val = thresholdsF1_in[i]
            max_pos = i
            if thresholdsF1_in[i] == 1:
                thresh_1s.append(i)
    #
    # if there is no value equal to 1 (maximum), or only 1
    if len(thresh_1s)==0 or len(thresh_1s)==1:
        return thresholds_in[max_pos], thresholdsF1_in[max_pos]
    # if we arrive here there are more than one perfect F1
    new_best = (max(thresh_1s) + min(thresh_1s)) / 2
    # we need also to check that the value in between has still a F1 of 1
    if measureF1(negative_vals,positive_vals,new_best)[0] == 1:
        return new_best, 1
    else:
        sys.stderr.write("Error. Cannot find best F1, will choose randomly between the best.")
        return min(thresh_1s), 1






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
                #
                thresholds_start = find_thresholds_this(negative_vals,positive_vals)
                thresholds_start_F1 = measureF1(negative_vals,positive_vals,thresholds_start)
                best_threshold,best_f1 = find_max(thresholds_start,thresholds_start_F1,negative_vals,positive_vals)
                res[i] = best_threshold
    return res


# ------------------------------------------------------------------------------
def find_thresholds(all_transformed, tax, intermediate_dist_for_NN):
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

        # check if saving the intermediate distances -------
        if(intermediate_dist_for_NN is not None):
            save_distances_to_file(dist_all_vs_centroids, clade, intermediate_dist_for_NN)

        # find the thresholds -------------------
        logging.info('    TRAIN_NN_3: Find thresholds')
        threshold_clades[clade] = find_thresholds_from_dist(dist_all_vs_centroids)

    return threshold_clades, list_centroids_all


#===============================================================================
#                                      MAIN
#===============================================================================
def train_NN_classifiers(alignment, tax_file, NN_start_level,logging_, verbose_,intermediate_dist_for_NN, min_training_data_lmnn, base_save_features, procs=1):
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
    all_LMNN, all_transformed, all_sel_positions = estimate_weights(alignment, tax, NN_start_level, min_training_data_lmnn, base_save_features, procs=procs)

    # 2. find centroids and find the threshold distances
    logging.info('  TRAIN_NN_1: find centroids and thresholds')
    if verbose > 4: sys.stderr.write("-- Find centroids and thresholds\n")
    thresholds_NN, centroid_seq = find_thresholds(all_transformed, tax, intermediate_dist_for_NN)

    return all_LMNN, thresholds_NN, centroid_seq, species_to_tax, all_sel_positions
