import argparse
import sys
import random
import logging
import os
import time
import tempfile
import shutil
from collections import Counter
import random
import pickle
import multiprocessing as mp

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from stag.taxonomy3 import Taxonomy
from stag.databases import save_to_file
from stag.alignment import load_alignment_from_file, EncodedAlignment
from stag.create_db import train_all_classifiers


def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("taxonomy", type=str)
	ap.add_argument("alignment", type=str)
	ap.add_argument("-o", "--output", type=str, default="stagdb")
	ap.add_argument('-e', action="store", default="l1", dest='penalty_logistic', help='penalty for the logistic regression',choices=['l1','l2','none'])
	ap.add_argument('-E', action="store", default="liblinear", dest='solver_logistic', help='solver for the logistic regression',choices=['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'])
	ap.add_argument("-t", "--threads", type=int, default=1)
	args = ap.parse_args()

	output, tax_file, aligned_seq_file, procs = args.output, args.taxonomy, args.alignment, args.threads
	solver_v, penalty_v = args.solver_logistic, args.penalty_logistic

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
	classifiers = [
		(node, np.append(clf.intercept_, clf.coef_) if clf else None)
		for node, clf in train_all_classifiers(alignment, full_taxonomy, penalty_v, solver_v, procs=procs)
	]
	with open(classifiers_file, "wb") as clf_out:
		pickle.dump(classifiers, clf_out)
	logging.info('TIME:Finished training all classifiers')

	logging.info('MAIN:Finished')


if __name__ == "__main__":
	main()
