# Input:
#  - one multiple sequence alignment (MSA) per marker gene. The MSA is obtained
#    from the function stag align, like:
#       >gene1\t0001010001010000100101000...
#       >gene2\t0000110001010100100101001...
#  - a taxonomy file that describes the taxonomy of the genes:
#       gene1\tBacteria;Firmicutes;...
#
# Output:
#  - a database file (hdf5) that can be used by stag classify

import sys
import random
import logging
import os
import time
import tempfile
import shutil
from collections import Counter

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
import h5py

from stag.taxonomy3 import Taxonomy
from stag.databases2 import save_to_file
from stag.alignment import load_alignment_from_file
from stag.train_NN import train_NN_classifiers


# function to check the NN level -----------------------------------------------
def check_level_NN_start_level(full_taxonomy, NN_start_level):
    n_level_tax = full_taxonomy.get_n_levels()
    # NN_start_level = 0 kingdom
    if NN_start_level >= n_level_tax-1:
        logging.info('Error the selected level for the NN classification is too high')
        print("Error: -L should be smaller than "+str(n_level_tax-1))
        sys.exit(1)

# function that finds positive and negative examples ===========================
def find_training_genes(node, siblings, full_taxonomy, alignment):
    t00 = time.time()
    # "positive_examples" and "negative_examples" are list of gene ids
    positive_examples = full_taxonomy.find_gene_ids(node)
    t_pos = time.time() - t00
    t0 = time.time()
    negative_examples = list()
    for s in siblings:
        negative_examples.extend(full_taxonomy.find_gene_ids(s))
    t_neg = time.time() - t0

    if not negative_examples:
        # it means that there was only one child, and we cannot do anything
        return positive_examples, negative_examples

    # From here, it means that there is at least one sibling ==================
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
    # but if there is only one other sibling, we choose only 3 times more negative
    max_negative_samples = len(positive_examples_subsample) * (20 if len(siblings) > 1 else 3)
    if len(negative_examples_subsample) > max_negative_samples:
        negative_examples_subsample = random.sample(negative_examples_subsample, max_negative_samples)
    # 4. we want to have at least 5 times more negative than positive ----------
    missing_neg = 0 # how many negative sequences we need to add
    min_negative_samples = len(positive_examples_subsample) * 5
    if len(negative_examples_subsample) < min_negative_samples:
        missing_neg = min_negative_samples - len(negative_examples_subsample)
    # add negative examples if necessary
    if missing_neg > 0:
        # positive examples
        X_clade = alignment.loc[positive_examples, : ].to_numpy()
        # always have 5 positive classes
        n_positive_class = len(X_clade)
        for i in range(n_positive_class, 5):
            rr = random.choice(range(n_positive_class))
            X_clade = np.vstack((X_clade, X_clade[rr,]))

        # find possible genes to add additionaly to negatives
        possible_neg = list(set(alignment.index.values).difference(set(positive_examples + negative_examples)))
        if possible_neg: # if it is possible to add negatives
                         # note that at the highest level, it's not possible
            X_poss_na = alignment.loc[possible_neg, : ].to_numpy()

            # choose 5 random positive clades
            X_check_sim = X_clade[random.sample(range(len(X_clade)), 5), ]

            random_clades = list()
            for clade_i in range(5):
                m_for_diff = np.tile(X_check_sim[clade_i,], (len(X_poss_na), 1))
                differences = np.sum(np.bitwise_xor(m_for_diff, X_poss_na), axis=1)
                non_zero = np.sum(differences != 0)
                differences = np.where(differences == 0, np.nan, differences)
                corr_ord = np.argsort(differences)[:non_zero + 1]
                random_clades.append(list(corr_ord))
            clade_indices = set()
            for indices in zip(*random_clades):
                clade_indices.update(indices)
                if len(clade_indices) > missing_neg:
                    break
            negative_examples_subsample.extend(possible_neg[i] for i in clade_indices)

    t_total = time.time() - t00
    logging.info(f"find_training_genes\t{node}\t{len(positive_examples)}\t{len(negative_examples)}\t{t_pos:.3f}\t{t_neg:.3f}\t{t_total:.3f}\t{os.getpid()}")

    return positive_examples_subsample, negative_examples_subsample


def get_classification_input(taxonomy, alignment):
    for node, siblings in taxonomy.get_all_nodes(get_root=True):
        logging.info(f'   TRAIN:"{node}":Find genes')
        positive_examples, negative_examples = find_training_genes(node, siblings, taxonomy, alignment)
        logging.info(f'      SEL_GENES:"{node}": {len(positive_examples)} positive, {len(negative_examples)} negative')

        # check that we have at least 1 example for each class:
        if not negative_examples:
            # when the node is the only child, then there are no negative examples
            logging.info('      Warning: no negative examples for "%s"', node)
            X, y = "no_negative_examples", None
        elif not positive_examples:
            # There should be positive examples
            logging.info('      Error: no positive examples for "%s"', node)
            X, y = "ERROR_no_positive_examples", None
        else:
            X = alignment.loc[ negative_examples + positive_examples , : ].to_numpy()
            y = np.asarray(["no"] * len(negative_examples) + ["yes"] * len(positive_examples))
        yield node, X, y

def train_all_classifiers_nonmp(alignment, full_taxonomy, penalty_v, solver_v, procs=None):
    #print("train_all_classifiers_nonmp - single-proc")
    all_classifiers = dict()
    for node, X, y in get_classification_input(full_taxonomy, alignment):
        if y is not None:
            clf = LogisticRegression(random_state=0, penalty = penalty_v, solver=solver_v)
            clf.fit(X, y)
            all_classifiers[node] = clf
        else:
            all_classifiers[node] = X
    return all_classifiers

def perform_training(X, y, penalty_v, solver_v, node):
    if y is None:
        return node, X
    # logging.info('         TRAIN:"%s":Train classifier', node)
    clf = LogisticRegression(random_state=0, penalty=penalty_v, solver=solver_v)
    clf.fit(X, y)
    return node, clf

def get_classification_input_mp(node, siblings, taxonomy, alignment, penalty_v, solver_v):
    logging.info(f'   TRAIN:"{node}":Find genes (proc={os.getpid})')
    positive_examples, negative_examples = find_training_genes(node, siblings, taxonomy, alignment)
    logging.info(f'      SEL_GENES:"{node}": {len(positive_examples)} positive, {len(negative_examples)} negative')

    # check that we have at least 1 example for each class:
    if not negative_examples:
        # when the node is the only child, then there are no negative examples
        logging.info('      Warning: no negative examples for "%s', node)
        X, y = "no_negative_examples", None
    elif not positive_examples:
        # There should be positive examples
        logging.info('      Error: no positive examples for "%s', node)
        X, y = "ERROR_no_positive_examples", None
    else:
        X = alignment.loc[ negative_examples + positive_examples , : ].to_numpy()
        y = np.asarray(["no"] * len(negative_examples) + ["yes"] * len(positive_examples))

    return perform_training(X, y, penalty_v, solver_v, node)

def get_classification_input_mp2(nodes, taxonomy, alignment, penalty_v, solver_v):
    results = list()
    for node, siblings in nodes:
        #logging.info(f'   TRAIN:"{node}":Find genes')
        t00 = time.time()
        positive_examples, negative_examples = find_training_genes(node, siblings, taxonomy, alignment)
        t_select = time.time() - t00
        #logging.info(f'      SEL_GENES:"{node}": {len(positive_examples)} positive, {len(negative_examples)} negative')

        # check that we have at least 1 example for each class:
        if not negative_examples:
            # when the node is the only child, then there are no negative examples
            logging.info('      Warning: no negative examples for "%s', node)
            X, y = "no_negative_examples", None
        elif not positive_examples:
            # There should be positive examples
            logging.info('      Error: no positive examples for "%s', node)
            X, y = "ERROR_no_positive_examples", None
        else:
            X = alignment.loc[ negative_examples + positive_examples , : ].to_numpy()
            y = np.asarray(["no"] * len(negative_examples) + ["yes"] * len(positive_examples))
        t0 = time.time()
        results.append(perform_training(X, y, penalty_v, solver_v, node))
        t1 = time.time()
        t_train, t_total = t1 - t0, t1 - t00

        # logging.info(f'   "{node}": {len(positive_examples)} positive, {len(negative_examples)} negative\tselection: {t_select:.3f}s, training: {t_train:.3f}s\tpid={os.getpid()}')
        logging.info("\t".join(map(str, [node, len(positive_examples), len(negative_examples), f"{t_select:.3f}s", f"{t_train:.3f}s", f"{t_total:.3f}s", os.getpid()])))
    return results

def train_all_classifiers_mp(alignment, full_taxonomy, penalty_v, solver_v, procs=2):
    import multiprocessing as mp
    print(f"train_all_classifiers_mp with {procs} processes.")
    logging.info("\t".join(["                  node", "positive", "negative", "t_select", "t_train", "t_total", "pid"]))
    with mp.Pool(processes=procs) as pool:
        nodes = list(full_taxonomy.get_all_nodes(get_root=False))
        step = len(nodes) // procs
        results = [
            pool.apply_async(get_classification_input_mp2, args=(nodes[i:i+step], full_taxonomy, alignment, penalty_v, solver_v))
            for i in range(0, len(nodes), step)
        ]

        res_d = dict()
        for res in results:
            res_d.update(res.get())
        return res_d

        #results = [
        #    pool.apply_async(get_classification_input_mp, args=(node, siblings, full_taxonomy, alignment, penalty_v, solver_v))
        #    for node, siblings in full_taxonomy.get_all_nodes(get_root=True)
        #]

        #results = [
        #    pool.apply_async(perform_training, args=(X, y, penalty_v, solver_v, node))
        #    for node, X, y in list(get_classification_input(full_taxonomy, alignment))
        #]

        # return dict(p.get() for p in results)

def train_all_classifiers(*args, procs=None):
    train_f = train_all_classifiers_mp if (procs and procs > 1) else train_all_classifiers_nonmp
    return train_f(*args, procs=procs)


#===============================================================================
#              FUNCTIONS TO LEARN THE FUNCTION FOR THE TAX. LEVEL
#===============================================================================
def predict_iter(test_seq, training_tax, classifiers_train, tax, perc, arrived_so_far):
    if training_tax.is_last_node(arrived_so_far):
        return
    max_perc, max_perc_taxa  = 0, ""
    children = training_tax.find_children_node(arrived_so_far)
    if not children:
        print("Error: no child", file=sys.stderr)
    elif len(children) == 1:
        # if there are no siblings I put 2, it will be replaced after
        max_perc, max_perc_taxa = 2, children[0]
    else:
        for child in children:
            clf = classifiers_train[child]
            predictions = clf.predict_proba(test_seq)
            predicted_proba = np.amin(predictions) if clf.predict(test_seq)[0] == "no" else np.amax(predictions)
            if predicted_proba > max_perc:
                max_perc, max_perc_taxa = predicted_proba, child

    tax.append(max_perc_taxa)
    perc.append(max_perc)

    predict_iter(test_seq, training_tax, classifiers_train, tax, perc, max_perc_taxa)


def predict_one_gene(test_seq, training_tax, classifiers_train):
    tax = list()
    perc = list()
    # we arrived at the root, and now we classify from there
    predict_iter(test_seq, training_tax, classifiers_train, tax, perc, training_tax.get_root())
    # we change the predictions that came from having only one sibling --------
    if perc[0] == 2:
        perc[0] = 1
    for i in range(len(perc)):
        if perc[i] == 2:
            perc[i] = perc[i-1]

    return tax, perc

def predict(test_al, training_tax, classifiers_train):
    return [
        [gene, *predict_one_gene([test_al.loc[ gene , : ].to_numpy()], training_tax, classifiers_train)]
        for gene in test_al.index.values
    ]

def learn_function(level_to_learn, alignment, full_taxonomy, penalty_v, solver_v, perc_test_set=0.33, gene_level=False, procs=None):
    # perc_test_set <= 0.5 !
    logging.info(f'  TEST:"{level_to_learn}" taxonomic level')
    # 1. Identify which clades we want to remove (test set) and which to keep (training set)
    test_set, training_set = set(), set()
    clades = full_taxonomy.get_last_level_to_genes() if gene_level else full_taxonomy.find_node_level(level_to_learn)
    for node, children in clades.items():
        aval_clades = set(children)
        # find how many to use for the test set:
        n_test = 0 if (not gene_level and len(aval_clades) == 2) else round(len(aval_clades) * perc_test_set)
        test_set.update(aval_clades.pop() for _ in range(n_test))
        training_set.update(aval_clades)
    logging.info(f'  TEST:"{level_to_learn}" level:test_set ({len(test_set)}):{test_set}')
    logging.info(f'  TEST:"{level_to_learn}" level:trai_set ({len(training_set)}):{training_set}')

    # 2. Create new taxonomy and alignment file & train the classifiers
    training_tax = full_taxonomy.copy()
    if gene_level:
        training_tax.remove_genes(list(test_set))
        training_filter = training_set
        test_filter = test_set
    else:
        test_filter = training_tax.remove_clades(list(test_set))
        training_filter = training_tax.find_gene_ids(training_tax.get_root())

    training_al = alignment.loc[ training_filter, : ]
    classifiers_train = train_all_classifiers(training_al, training_tax, penalty_v, solver_v, procs=procs)

    # 3. Classify the test set
    test_al = alignment.loc[ test_filter , : ]
    pr = predict(test_al, training_tax, classifiers_train)
    for g in pr:
        # g is:
        # ["geneB",["A","B","D","species8"],[0.99,0.96,0.96,0.07]]
        correct_tax = full_taxonomy.extract_full_tax_from_gene(g[0])
        g.extend([correct_tax, level_to_learn])

    return pr
    # return:
    #    GENE_ID         PREDICTED             PROB_PREDICTED        CORRECT        REMOVED_LEVEL
    # [["geneA",["A","B","C","species2"],[0.98,0.97,0.23,0.02],["A","B","Y","speciesX"],2]
    #  ["geneB",["A","B","D","species8"],[0.99,0.96,0.10,0.07],["A","B","U","speciesZ"],2]
    # .....                                                                               ]

def estimate_function(all_calc_functions):
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

    # we remove duplicates with the same predicted probability -----------------
    all_uniq = {tuple(round(v, 2) for v in item[2]): item for item in all_calc_functions}
    logging.info(f'   LEARN_FUNCTION:Number of lines: {len(all_uniq)}/{len(all_calc_functions)}')

    correct_level = list()
    for _, predicted, _, ground_truth, _ in all_uniq.values():
        corr_level_this = -1
        for cont, (p, c) in enumerate(zip(predicted, ground_truth)):
            if p == c:
                corr_level_this = cont # we select to what level to predict
        correct_level.append(corr_level_this)

    # now in correct_level there is to which level to predict to. Example:
    # "A","B","C","species2"
    # with corr_level_this = 0, we should assign "A"
    # with corr_level_this = 2, we should assign "A","B","C"
    # with corr_level_this = -1, we should assign "" (no taxonomy)

    level_counter = Counter(correct_level)
    for level, count in sorted(level_counter.items()):
        logging.info(f'   LEARN_FUNCTION:Number of lines: level {level}: {count}')

    all_classifiers = dict()
    for uniq_level in sorted(level_counter):
        # NOTE: we always need the negative class to be first
        correct_order = [[], []]
        for level, (_, _, prob, *_) in zip(correct_level, all_uniq.values()):
            correct_order[int(uniq_level == level)].append(prob)

        X = np.array([np.array(xi) for xi in correct_order[0] + correct_order[1]])
        y = np.asarray([0] * len(correct_order[0]) + [1] * len(correct_order[1]))
        clf = LogisticRegression(random_state=0, penalty = "none", solver='saga', max_iter = 5000)
        clf.fit(X, y)
        all_classifiers[str(uniq_level)] = clf

    return all_classifiers


# create taxonomy selection function ===========================================
# This function define a function that is able to identify to which taxonomic
# level a new gene should be assigned to.
def learn_taxonomy_selection_function(alignment, full_taxonomy, save_cross_val_data, penalty_v, solver_v, procs=None):
    # find number of levels
    n_levels = full_taxonomy.get_n_levels()

    # do the cross validation for each level
    all_calc_functions = list()
    for level in range(n_levels):
        all_calc_functions.extend(learn_function(level, alignment, full_taxonomy, penalty_v, solver_v, procs=procs))
        #all_calc_functions.extend(learn_function_one_level(level, alignment, full_taxonomy, penalty_v, solver_v, procs=procs))
    # do the cross val. for the last level (using the genes)
    # all_calc_functions.extend(learn_function_genes_level(n_levels, alignment, full_taxonomy, penalty_v, solver_v, procs=procs))
    all_calc_functions.extend(learn_function(n_levels, alignment, full_taxonomy, penalty_v, solver_v, gene_level=True, procs=procs))

    # save all_calc_functions if necessary -------------------------------------
    if save_cross_val_data:
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        with outfile:
            os.chmod(outfile.name, 0o644)
            print("gene", "predicted", "prob", "ground_truth", "removed_level", sep="\t", file=outfile)
            for gene, predicted, prob, ground_truth, removed_level in all_calc_functions:
                predicted, prob, ground_truth = ("/".join(s) for s in (predicted, ["{:.2f}".format(pr) for pr in prob], ground_truth))
                print(gene, predicted, prob, ground_truth, removed_level, sep="\t", file=outfile)
            try:
                outfile.flush()
                os.fsync(outfile.fileno())
            except:
                print("[E::main] Error: failed to save the cross validation results", file=sys.stderr)
        try:
            shutil.move(outfile.name, save_cross_val_data)
        except:
            print("[E::main] Error: failed to save the cross validation results\n" + \
                  f"[E::main] you can find the file here: \n{outfile.name}\n", file=sys.stderr)

    return estimate_function(all_calc_functions)


def create_db(aligned_seq_file, tax_file, verbose, output, use_cmalign, hmm_file_path, save_cross_val_data, protein_fasta_input, penalty_v, solver_v, NN_start_level, procs=None):
    filename_log = os.path.realpath(output)+'.log'
    logging.basicConfig(filename=filename_log,
                        filemode='w',
                        level=logging.INFO,
                        format='[%(asctime)s] %(message)s')
    logging.info('TIME:start')

    # 1. load the taxonomy into the tree (global variable)
    logging.info('MAIN:Load taxonomy')
    full_taxonomy = Taxonomy(tax_file)
    full_taxonomy.load_from_file()
    logging.info('TIME:Finish load taxonomy')

    # 2. load the alignment into a pandas dataframe
    logging.info('MAIN:Load alignment')
    alignment = load_alignment_from_file(aligned_seq_file)
    logging.info('TIME:Finish load alignment')

    # 3. check that the taxonomy and the alignment are consistent
    logging.info('MAIN:Check taxonomy and alignment')
    full_taxonomy.ensure_geneset_consistency(list(alignment.index.values))
    check_level_NN_start_level(full_taxonomy, NN_start_level)
    logging.info('TIME:Finish check-up')

    # 4. build a classifier for each node
    logging.info('MAIN:Train all classifiers')
    classifiers = train_all_classifiers(alignment, full_taxonomy, penalty_v, solver_v, procs=procs)
    logging.info('TIME:Finish train all classifiers')

    # 5. learn the function to identify the correct taxonomy level
    logging.info('MAIN:Learn taxonomy selection function')
    tax_function = learn_taxonomy_selection_function(alignment, full_taxonomy, save_cross_val_data, penalty_v, solver_v, procs=procs)
    logging.info('TIME:Finish learn taxonomy selection function')

    # 6. train classifiers for the nearest neighbour
    logging.info('MAIN:Train classifiers for nearest neighbour')
    all_LMNN, thresholds_NN, centroid_seq, species_to_tax, all_sel_positions = train_NN_classifiers(alignment, tax_file, NN_start_level, logging)
    logging.info('TIME:Finish train classifiers for nearest neighbour')

    # 7. save the result
    logging.info('MAIN:Save to file')
    save_to_file(classifiers, full_taxonomy, tax_function, use_cmalign, output, all_LMNN, thresholds_NN, centroid_seq, species_to_tax,all_sel_positions,NN_start_level, hmm_file_path=hmm_file_path, protein_fasta_input=protein_fasta_input)
    logging.info('TIME:Finish save to file')

    logging.info('MAIN:Finished')
