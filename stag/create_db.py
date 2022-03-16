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

import logging
import multiprocessing as mp
import os
import pickle
import random
import shutil
import sys
import tempfile
import time

from collections import Counter

import numpy as np
from sklearn.linear_model import LogisticRegression

from stag.taxonomy3 import Taxonomy
from stag.databases import save_to_file
from stag.alignment import EncodedAlignment

import warnings
warnings.filterwarnings("error")


class BalancingParameters:
    # 1. max 500 positive samples ----------------------------------------------
    MAX_POSITIVE_SAMPLES = 500
    # 2. max 1000 negative samples ---------------------------------------------
    MAX_NEGATIVE_SAMPLES = 1000
    # 3. max 20 times more negative than positive ------------------------------
    # but if there is only one other sibling, we choose only 3 times more negative
    NEG_POS_FACTOR_1_SIBLING = 3
    NEG_POS_FACTOR_N_SIBLINGS = 20
    # 4. we want to have at least 5 times more negative than positive ----------
    MIN_NEGATIVE_SAMPLE_FACTOR = 5

    @staticmethod
    def get_neg_pos_sample_factor(n_siblings):
        if n_siblings > 1:
            return BalancingParameters.NEG_POS_FACTOR_N_SIBLINGS

        return BalancingParameters.NEG_POS_FACTOR_1_SIBLING

    @staticmethod
    def get_missing_negatives(n_positives, n_negatives):
        return max(0, n_positives * BalancingParameters.MIN_NEGATIVE_SAMPLE_FACTOR - n_negatives)


def balance_neg_pos_ratio(positive_examples, negative_examples, alignment, n_siblings):

    if len(positive_examples) > BalancingParameters.MAX_POSITIVE_SAMPLES:
        positive_examples = random.sample(positive_examples, BalancingParameters.MAX_POSITIVE_SAMPLES)
    n_positives = len(positive_examples)

    max_negatives = min(
        BalancingParameters.MAX_NEGATIVE_SAMPLES,
        BalancingParameters.get_neg_pos_sample_factor(n_siblings) * n_positives
    )

    if len(negative_examples) > max_negatives:
        negative_examples = random.sample(negative_examples, max_negatives)

    # i still think there should be an else here,
    # otherwise we downsample, then upsample again in some (hypothetical) cases
    missing_negatives = BalancingParameters.get_missing_negatives(
        n_positives, len(negative_examples)
    )
    if missing_negatives:
        negative_examples += alignment.get_auxiliary_negative_examples(
            positive_examples, negative_examples, missing_negatives
        )

    return positive_examples, negative_examples


# function that finds positive and negative examples ===========================
def find_training_genes(node, siblings, full_taxonomy, alignment):
    t00 = time.time()
    # "positive_examples" and "negative_examples" are list of gene ids
    positive_examples = full_taxonomy.find_gene_ids(node)
    t_pos = time.time() - t00
    t0 = time.time()
    negative_examples = list()
    for s in siblings:
        negative_examples += full_taxonomy.find_gene_ids(s)
    t_neg = time.time() - t0

    # if there are negative examples, then balance negative/positive ratio
    if negative_examples:
        positive_examples, negative_examples = balance_neg_pos_ratio(
            positive_examples, negative_examples, alignment, len(siblings)
        )

    t_total = time.time() - t00
    logging.info(
        f"find_training_genes\t{node}\t{len(positive_examples)}\t{len(negative_examples)}"
        f"\t{t_pos:.3f}\t{t_neg:.3f}\t{t_total:.3f}\t{os.getpid()}"
    )

    return positive_examples, negative_examples


def train_all_classifiers_nonmp(alignment, taxonomy, penalty_v, solver_v, max_iter=5000, procs=None):
    return list(
        do_training(node, siblings, taxonomy, alignment, penalty_v, solver_v, max_iter)
        for node, siblings in taxonomy.get_all_nodes(get_root=False)
    )


def do_training(node, siblings, taxonomy, alignment, penalty_v, solver_v, max_iter=5000):
    t00 = time.time()
    positive_examples, negative_examples = find_training_genes(node, siblings, taxonomy, alignment)
    t_select, t_train = time.time() - t00, 0

    if positive_examples and negative_examples:
        t0 = time.time()
        clf = LogisticRegression(random_state=0, penalty=penalty_v, solver=solver_v, max_iter=max_iter)
        clf.fit(
            alignment.get_rows(node, negative_examples + positive_examples),
            np.asarray(["no"] * len(negative_examples) + ["yes"] * len(positive_examples))
        )
        t_train = time.time() - t0
    elif positive_examples:
        msg, clf = f'      Warning: no negative examples for "{node}"', None
        logging.info(msg)
    else:
        raise ValueError(f'This literally cannot happen: no positive examples for "{node}"!')

    t_total = time.time() - t00
    logging.info(
        f"{node}\t{len(positive_examples)}\t{len(negative_examples)}"
        f"{t_select:.3f}s\t{t_train:.3f}s\t{t_total:.3f}s\t{os.getpid()}"
    )

    return node, clf


def process_chunk(nodes, taxonomy, alignment, penalty_v, solver_v):
    return [
        do_training(node, siblings, taxonomy, alignment, penalty_v, solver_v)
        for node, siblings in nodes
    ]


def train_all_classifiers_mp(alignment, full_taxonomy, penalty_v, solver_v, max_iter=5000, procs=2):
    print(f"train_all_classifiers_mp with {procs} processes.")
    logging.info("\t".join([
        "                  node", "positive", "negative", "t_select", "t_train", "t_total", "pid"
    ]))
    with mp.Pool(processes=procs) as pool:
        nodes = list(full_taxonomy.get_all_nodes(get_root=False))
        # try to distribute the load a little (the higher nodes have more data attached)
        # but still keep it reproducible, hence reset the seed every time
        random.seed(313)
        random.shuffle(nodes)
        step = len(nodes) // procs + 1

        results = pool.starmap_async(
            do_training,
            [(node, siblings, full_taxonomy, alignment, penalty_v, solver_v, max_iter) for node, siblings in nodes],
            step)
        return list(results.get())


def train_all_classifiers(*args, max_iter=5000, procs=None):
    train_f = train_all_classifiers_mp if (procs and procs > 1) else train_all_classifiers_nonmp

    return train_f(*args, max_iter=5000, procs=procs)


# ===============================================================================
#              FUNCTIONS TO LEARN THE FUNCTION FOR THE TAX. LEVEL
# ===============================================================================
def predict_iter(test_seq, training_tax, classifiers_train, tax, perc, arrived_so_far):
    if training_tax.is_last_node(arrived_so_far):
        return
    max_perc, max_perc_taxa = 0, ""
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
        [gene, *predict_one_gene(test_al.get_rows(gene, [gene]), training_tax, classifiers_train)]
        for gene in test_al.alignment.index.values
    ]


def learn_function(
    level_to_learn, alignment, full_taxonomy, penalty_v, solver_v,
    perc_test_set=0.33, gene_level=False, max_iter=5000, procs=None
):
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
        training_filter = training_tax.find_gene_ids()  # training_tax.get_root())

    logging.info(f'  TEST: training_taxonomy has {len(training_tax)} nodes.')

    classifiers_train = dict(
        train_all_classifiers(
            alignment.filter_alignment(training_filter),
            training_tax,
            penalty_v, solver_v, max_iter=max_iter, procs=procs
        )
    )

    # 3. Classify the test set
    pr = predict(
        alignment.filter_alignment(test_filter),
        training_tax,
        classifiers_train
    )

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


def estimate_function(all_calc_functions, max_iter=5000):
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
                corr_level_this = cont  # we select to what level to predict
        correct_level.append(corr_level_this)

    # now in correct_level there is to which level to predict to. Example:
    # "A","B","C","species2"
    # with corr_level_this = 0, we should assign "A"
    # with corr_level_this = 2, we should assign "A","B","C"
    # with corr_level_this = -1, we should assign "" (no taxonomy)
    level_counter = Counter(correct_level)
    for level, count in sorted(level_counter.items()):
        logging.info(f'   LEARN_FUNCTION:Number of lines: level {level}: {count}')

    all_classifiers = list()
    for uniq_level in sorted(level_counter):
        # NOTE: we always need the negative class to be first
        correct_order = [[], []]
        for level, (_, _, prob, *_) in zip(correct_level, all_uniq.values()):
            correct_order[int(uniq_level == level)].append(prob)

        if correct_order[0] and correct_order[1]:
            X = np.array([np.array(xi) for xi in correct_order[0] + correct_order[1]])
            y = np.asarray([0] * len(correct_order[0]) + [1] * len(correct_order[1]))
            clf = LogisticRegression(random_state=0, penalty="none", solver='saga', max_iter=max_iter)
            clf.fit(X, y)
            all_classifiers.append((str(uniq_level), clf))
        else:
            logging.info(
                f'Could not train classifier {uniq_level}: '
                f'neg={len(correct_order[0])} pos={len(correct_order[1])}'
            )

    return all_classifiers


# create taxonomy selection function ===========================================
# This function define a function that is able to identify to which taxonomic
# level a new gene should be assigned to.
def learn_taxonomy_selection_function(
    alignment, full_taxonomy, save_cross_val_data, penalty_v, solver_v, max_iter=5000, procs=None
):
    # find number of levels
    n_levels = full_taxonomy.get_n_levels()

    # do the cross validation for each level
    all_calc_functions = list()
    for level in range(n_levels):
        all_calc_functions += learn_function(
            level, alignment, full_taxonomy, penalty_v, solver_v, procs=procs
        )
    # do the cross val. for the last level (using the genes)
    all_calc_functions += learn_function(
        n_levels, alignment, full_taxonomy, penalty_v, solver_v, gene_level=True, procs=procs
    )

    # save all_calc_functions if necessary -------------------------------------
    if save_cross_val_data:
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        with outfile:
            os.chmod(outfile.name, 0o644)
            print(
                "gene", "predicted", "prob", "ground_truth", "removed_level",
                sep="\t", file=outfile
            )
            for gene, predicted, prob, ground_truth, removed_level in all_calc_functions:
                predicted, prob, ground_truth = (
                    "/".join(s)
                    for s in (predicted, [f"{pr:.2f}" for pr in prob], ground_truth)
                )
                print(gene, predicted, prob, ground_truth, removed_level, sep="\t", file=outfile)
            try:
                outfile.flush()
                os.fsync(outfile.fileno())
            except Exception:
                print("[E::main] Error: failed to save the cross validation results", file=sys.stderr)
        try:
            shutil.move(outfile.name, save_cross_val_data)
        except Exception:
            print(
                "[E::main] Error: failed to save the cross validation results\n"
                f"[E::main] you can find the file here: \n{outfile.name}\n",
                file=sys.stderr
            )

    return estimate_function(all_calc_functions, max_iter=max_iter)


def create_db(
    aligned_seq_file, tax_file, verbose, output, use_cmalign,
    hmm_file_path, save_cross_val_data, protein_fasta_input,
    penalty_v, solver_v, max_iter=5000, procs=None
):
    filename_log = os.path.realpath(output)+'.log'
    logging.basicConfig(filename=filename_log,
                        filemode='w',
                        level=logging.INFO,
                        format='[%(asctime)s] %(message)s')
    logging.info('TIME:start')

    # 1. load the taxonomy into the tree (global variable)
    logging.info('MAIN:Loading taxonomy')
    full_taxonomy = Taxonomy(tax_file)
    logging.info(f'TIME:Finished loading taxonomy - {len(full_taxonomy)} nodes in taxonomy')

    # 2. load the alignment into a pandas dataframe
    logging.info('MAIN:Loading alignment')
    alignment = EncodedAlignment(aligned_seq_file)
    logging.info('TIME:Finished loading alignment')

    # 3. check that the taxonomy and the alignment are consistent
    logging.info('MAIN:Checking taxonomy and alignment')
    full_taxonomy.ensure_geneset_consistency(alignment.get_index())
    logging.info(f'TIME:Finished check-up - {len(full_taxonomy)} nodes in taxonomy')

    # 4. build a classifier for each node
    logging.info('MAIN:Training all classifiers')
    classifiers_file = output + ".classifiers.dat"
    if all((
        os.path.exists(f) for f in (classifiers_file, classifiers_file + ".ok")
    )):
        with open(classifiers_file, "rb") as classifiers_in:
            classifiers = pickle.load(classifiers_in)
    else:
        classifiers = train_all_classifiers(
            alignment, full_taxonomy, penalty_v, solver_v,
            max_iter=max_iter, procs=procs
        )
        classifiers = [
            (node, np.append(clf.intercept_, clf.coef_) if clf else None)
            for node, clf in classifiers
        ]
        with open(classifiers_file, "wb") as clf_out:
            pickle.dump(classifiers, clf_out)
        open(classifiers_file + ".ok", "w").close()
    logging.info('TIME:Finished training all classifiers')

    # 5. learn the function to identify the correct taxonomy level
    logging.info('MAIN:Learning taxonomy selection function')
    taxfunc_file = output + ".taxfunc.dat"
    if all((os.path.exists(f) for f in (taxfunc_file, taxfunc_file + ".ok"))):
        tax_function = pickle.load(open(taxfunc_file, "rb"))
    else:
        learned_function = learn_taxonomy_selection_function(
            alignment, full_taxonomy, save_cross_val_data, penalty_v,
            solver_v, max_iter=max_iter, procs=procs
        )
        tax_function = [
            (node, np.append(clf.intercept_, clf.coef_) if clf else None)
            for node, clf in learned_function
        ]
        with open(taxfunc_file, "wb") as clf_out:
            pickle.dump(tax_function, clf_out)
        open(taxfunc_file + ".ok", "w").close()

    logging.info('TIME:Finished learning taxonomy selection function')

    # 6. save the result
    logging.info('MAIN:Saving database to file')
    save_to_file(
        classifiers, full_taxonomy, tax_function, use_cmalign, output,
        hmm_file_path=hmm_file_path, protein_fasta_input=protein_fasta_input
    )
    logging.info('TIME:Finished saving database')

    logging.info('MAIN:Finished')
