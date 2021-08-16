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
from stag.create_db import estimate_function

"""
def create_db(aligned_seq_file, tax_file, verbose, output, use_cmalign, hmm_file_path, save_cross_val_data, protein_fasta_input, penalty_v, solver_v, procs=None):
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
	if all((os.path.exists(f) for f in (classifiers_file, classifiers_file + ".ok"))):
		classifiers = pickle.load(open(classifiers_file, "rb"))
	else:
		classifiers = [
			(node, np.append(clf.intercept_, clf.coef_) if clf else None)
			for node, clf in train_all_classifiers(alignment, full_taxonomy, penalty_v, solver_v, procs=procs)
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
		tax_function = [
			(node, np.append(clf.intercept_, clf.coef_) if clf else None)
			for node, clf in learn_taxonomy_selection_function(alignment, full_taxonomy, save_cross_val_data, penalty_v, solver_v, procs=procs)
		]
		with open(taxfunc_file, "wb") as clf_out:
			pickle.dump(tax_function, clf_out)
		open(taxfunc_file + ".ok", "w").close()

	logging.info('TIME:Finished learning taxonomy selection function')

	# 6. save the result
	logging.info('MAIN:Saving database to file')
	save_to_file(classifiers, full_taxonomy, tax_function, use_cmalign, output, hmm_file_path=hmm_file_path, protein_fasta_input=protein_fasta_input)
	logging.info('TIME:Finished saving database')

	logging.info('MAIN:Finished')


"""

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
	#with tempfile.NamedTemporaryFile(delete=False, mode="w") as outfile:
	#	os.chmod(outfile.name, 0o644)
		print("gene", "predicted", "prob", "ground_truth", "removed_level", sep="\t", file=outfile)
		for gene, predicted, prob, ground_truth, removed_level in tax_function:
			predicted, prob, ground_truth = ("/".join(s) for s in (predicted, ["{:.2f}".format(pr) for pr in prob], ground_truth))
			print(gene, predicted, prob, ground_truth, removed_level, sep="\t", file=outfile)
	#	try:
	#		outfile.flush()
	#		os.fsync(outfile.fileno())
	#	except:
	#		print("[E::main] Error: failed to save the cross validation results", file=sys.stderr)
	#try:
	#	shutil.move(outfile.name, save_cross_val_data)
	#except:
	#	print("[E::main] Error: failed to save the cross validation results\n" + \
	#		  f"[E::main] you can find the file here: \n{outfile.name}\n", file=sys.stderr)

	tax_function = estimate_function(tax_function)

	classifiers = pickle.load(open(args.classifiers, "rb"))
	logging.info('TIME:Finished loading/estimating')

	# 4. save the result
	logging.info('MAIN:Saving database to file')
	save_to_file(classifiers, full_taxonomy, tax_function, args.use_cmalign, f"{output}.stagDB", hmm_file_path=args.hmm_file_path, protein_fasta_input=args.protein_fasta_input)
	logging.info('TIME:Finished saving database')
		

	logging.info('MAIN:Finished')





"""
# create taxonomy selection function ===========================================
# This function define a function that is able to identify to which taxonomic
# level a new gene should be assigned to.
def learn_taxonomy_selection_function(alignment, full_taxonomy, save_cross_val_data, penalty_v, solver_v, procs=None):
	# find number of levels

	# do the cross validation for each level
	all_calc_functions = list()
	for level in range(n_levels):
		all_calc_functions += learn_function(level, alignment, full_taxonomy, penalty_v, solver_v, procs=procs)
	# do the cross val. for the last level (using the genes)
	all_calc_functions += learn_function(n_levels, alignment, full_taxonomy, penalty_v, solver_v, gene_level=True, procs=procs)

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
"""




if __name__ == "__main__":
	main()
