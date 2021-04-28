import sys
import time
import os
import tempfile
import shutil
import contextlib

import numpy as np
import h5py

from . import __version__ as tool_version
import stag.align as align


def load_genome_DB(database, tool_version, verbose):
    dirpath = tempfile.mkdtemp()
    shutil.unpack_archive(database, dirpath, "gztar")
    list_files = [f for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, f))]
    for f in ("threshold_file.tsv", "hmm_lengths_file.tsv", "concatenated_genes_STAG_database.HDF5"):
        if f not in list_files:
            raise ValueError(f"[E::align] Error: {f} is missing.")

    with open(os.path.join(dirpath, "threshold_file.tsv")) as threshold_in:
        gene_thresholds = dict(line.rstrip().split("\t") for line in threshold_in)
        gene_order = list(gene_thresholds.keys())

    with open(os.path.join(dirpath, "hmm_lengths_file.tsv")) as hmm_lengths_in:
        ali_lengths = dict(line.rstrip().split("\t") for line in hmm_lengths_in)

    list_files.remove("threshold_file.tsv")
    list_files.remove("hmm_lengths_file.tsv")
    list_files.remove("concatenated_genes_STAG_database.HDF5")
    return list_files, dirpath, gene_thresholds, gene_order, ali_lengths, os.path.join(dirpath, "concatenated_genes_STAG_database.HDF5")


def load_db(hdf5_DB_path, protein_fasta_input=None, aligned_sequences=None, dir_output=None):
    with h5py.File(hdf5_DB_path, 'r') as db_in:
        params_out = open(os.path.join(dir_output, "parameters.tsv"), "w") if dir_output else contextlib.nullcontext()
        with params_out:
            # zero: tool version -------------------------------------------------------
            db_tool_version = db_in['tool_version'][0]
            use_proteins = db_in['align_protein'][0] # bool
            use_cmalign = db_in['use_cmalign'][0] # bool
            if dir_output:
                params_out.write("Tool version: "+str(db_tool_version)+"\n")
                params_out.write("Use proteins for the alignment: "+str(use_proteins)+"\n")
                params_out.write("Use cmalign instead of hmmalign: "+str(use_cmalign)+"\n")

        # check that it is the correct database, for 'classify', we need a single
        # gene
        if db_in['db_type'][0] != "single_gene":
            sys.stderr.write("[E::main] Error: this database is not designed to run with stag classify\n")
            sys.exit(1)
        # check if we used proteins
        if not aligned_sequences and not dir_output:
            if protein_fasta_input and not db_in['align_protein'][0]:
                # some proteins are provided in the classify but the db was constructed without using the proteins
                raise ValueError("Protein provided, but the database was constructed on genes.\n")
            elif not protein_fasta_input and db_in['align_protein'][0]:
                # the classify do not have proteins but the db was constructed WITH the proteins
                raise ValueError("Missing protein file (the database was constructed aligning proteins).\n")

        # first, we save a temporary file with the hmm file ------------------------
        hmm_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
        with hmm_file:
            os.chmod(hmm_file.name, 0o644)
            hmm_file.write(db_in['hmm_file'][0])
            hmm_file.flush()
            os.fsync(hmm_file.fileno())

        if dir_output:
            shutil.move(hmm_file.name, os.path.join(dir_output, "hmmfile.hmm"))

        # second if we need to use cm_align ----------------------------------------
        use_cmalign = db_in['use_cmalign'][0] # bool

        # third: taxonomy ----------------------------------------------------------
        taxonomy = {key: list(db_in['taxonomy/{}'.format(key)]) for key in db_in['taxonomy']}
        if dir_output:
            tax_out = open(os.path.join(dir_output, "node_hierarchy.tsv"), "w")
            with tax_out:
                print("Node", "Children", sep="\t", file=tax_out)
                for key, values in taxonomy.items():
                    print(key, *map(str, values), sep="\t", file=tax_out)

        # fourth: tax_function -----------------------------------------------------
        tax_function = {str(key): np.array(db_in['tax_function/{}'.format(key)], dtype=np.float64) 
                        for key in db_in['tax_function']}
        if dir_output:
            tax_func_out = open(os.path.join(dir_output, "taxonomy_function.tsv"), "w")
            with tax_func_out:
                for key, value in tax_function.items():
                    print(key, value, sep="\t", file=tax_func_out)

        # fifth: the classifiers ---------------------------------------------------
        classifiers = dict()
        class_out = open(os.path.join(dir_output, "classifiers_weights.tsv"), "w") if dir_output else contextlib.nullcontext()
        with class_out:
            for key in db_in['classifiers']:
                classifier = db_in['classifiers/{}'.format(key)]
                if not isinstance(classifier[0], str):
                    classifiers[key] = np.array(classifier, dtype=np.float64) 
                else:
                    classifiers[key] = "no_negative_examples"
                if dir_output:
                    print(key, *classifiers[key], sep="\t", file=class_out)

    return hmm_file.name, use_cmalign, taxonomy, tax_function, classifiers, db_tool_version


def save_to_file(classifiers, full_taxonomy, tax_function, use_cmalign, output, hmm_file_path=None, protein_fasta_input=None):

    string_dt = h5py.special_dtype(vlen=str)
    with h5py.File(output, "w") as h5p_out:
        # zero: tool version -------------------------------------------------------
        h5p_out.create_dataset('tool_version', data=np.array([str(tool_version)], "S100"), dtype=string_dt)
        # and type of database
        h5p_out.create_dataset('db_type', data=np.array(["single_gene"], "S100"), dtype=string_dt)
        # was the alignment done at the protein level?
        h5p_out.create_dataset('align_protein', data=np.array([bool(protein_fasta_input)]), dtype=bool)
        # first we save the hmm file -----------------------------------------------
        hmm_string = "".join(line for line in open(hmm_file_path)) if hmm_file_path else "NA"
        h5p_out.create_dataset('hmm_file', data=np.array([hmm_string], "S" + str(len(hmm_string) + 100)), dtype=string_dt, compression="gzip")
        # second, save the use_cmalign info ----------------------------------------
        h5p_out.create_dataset('use_cmalign', data=np.array([use_cmalign]), dtype=bool)
        # third, we save the taxonomy ---------------------------------------------
        h5p_out.create_group("taxonomy")
        for node, _ in full_taxonomy.get_all_nodes(get_root=True):
            h5p_out.create_dataset(f"taxonomy/{node}", data=np.array(list(full_taxonomy[node].children.keys()), "S10000"), dtype=string_dt, compression="gzip")
        # fourth, the taxonomy function --------------------------------------------
        h5p_out.create_group("tax_function")
        for c in tax_function:
            # we append the intercept at the head (will have position 0)
            vals = np.append(tax_function[c].intercept_, tax_function[c].coef_)
            h5p_out.create_dataset("tax_function/" + str(c), data=vals, dtype=np.float64, compression="gzip")
        # fifth, save the classifiers ----------------------------------------------
        h5p_out.create_group("classifiers")
        for c in classifiers:
            if classifiers[c] != "no_negative_examples":
                vals = np.append(classifiers[c].intercept_, classifiers[c].coef_)
                h5p_out.create_dataset("classifiers/" + c, data=vals, dtype=np.float64, compression="gzip", compression_opts=8)
            else:
                # in this case, it always predict 1, we save it as an array of
                # with the string "no_negative_examples"
                h5p_out.create_dataset("classifiers/" + c, data=np.array(["no_negative_examples"], "S40"), dtype=string_dt, compression="gzip")

        h5p_out.flush()
