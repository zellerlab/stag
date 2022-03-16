import argparse
import logging
import os
import pickle

import numpy as np

from stag.taxonomy3 import Taxonomy
from stag.databases import save_to_file
from stag.alignment import EncodedAlignment
from stag.create_db import estimate_function


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("taxonomy", type=str)
    ap.add_argument("alignment", type=str)
    ap.add_argument("classifiers", type=str)
    ap.add_argument("hmm_file_path", type=str)
    ap.add_argument("protein_fasta_input", type=str)
    ap.add_argument("levels", nargs="*")
    ap.add_argument("-o", "--output", type=str, default="stagdb")
    ap.add_argument("--use_cmalign", action="store_true")
    args = ap.parse_args()

    output, tax_file, aligned_seq_file = args.output, args.taxonomy, args.alignment

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

    logging.info('MAIN:Loading trained data and estimating learn function')
    tax_function = list()
    for f in sorted(args.levels):
        tax_function += pickle.load(open(f, "rb"))

    with open(f"{output}.cross_val", "w") as outfile:
        print("gene", "predicted", "prob", "ground_truth", "removed_level", sep="\t", file=outfile)
        for gene, predicted, prob, ground_truth, removed_level in tax_function:
            predicted, prob, ground_truth = (
                "/".join(s)
                for s in (predicted, ["{:.2f}".format(pr) for pr in prob], ground_truth)
            )
            print(gene, predicted, prob, ground_truth, removed_level, sep="\t", file=outfile)

    tax_function = [
        (node, np.append(clf.intercept_, clf.coef_) if clf else None)
        for node, clf in estimate_function(tax_function)
    ]

    classifiers = pickle.load(open(args.classifiers, "rb"))
    logging.info('TIME:Finished loading/estimating')

    # 4. save the result
    logging.info('MAIN:Saving database to file')
    save_to_file(
        classifiers, full_taxonomy, tax_function, args.use_cmalign,
        f"{output}.stagDB", hmm_file_path=args.hmm_file_path, protein_fasta_input=args.protein_fasta_input
    )
    logging.info('TIME:Finished saving database')

    logging.info('MAIN:Finished')


if __name__ == "__main__":
    main()
