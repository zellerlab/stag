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
