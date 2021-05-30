import sys
import time
import os
import tempfile
import shutil
import contextlib

import numpy as np
import pandas as pd

import pickle
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
    with open(hdf5_DB_path, "rb") as f:
        ALL_DATA = pickle.load(f)


    # we save a temporary file with the hmm file -------------------------------
    hmm_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    with hmm_file:
        os.chmod(hmm_file.name, 0o644)
        hmm_file.write(ALL_DATA["hmm_file"])
        hmm_file.flush()
        os.fsync(hmm_file.fileno())

    # we check if saving the database somewhere --------------------------------
    if dir_output:
        params_out = open(os.path.join(dir_output, "parameters.tsv"), "w") if dir_output else contextlib.nullcontext()

    return hmm_file.name, ALL_DATA["use_cmalign"], ALL_DATA["taxonomy"], ALL_DATA["tax_function"], ALL_DATA["classifiers"], ALL_DATA["tool_version"]


def save_to_file(classifiers, full_taxonomy, tax_function, use_cmalign, output, all_LMNN, thresholds_NN, centroid_seq, species_to_tax, hmm_file_path=None, protein_fasta_input=None):
    # we create a dict with all the data we need to save
    ALL_DATA = dict()
    # save all
    ALL_DATA["tool_version"] = str(tool_version)
    ALL_DATA["db_type"] = "single_gene"
    ALL_DATA["align_protein"] = bool(protein_fasta_input)
    ALL_DATA["use_cmalign"] = use_cmalign

    hmm_string = "".join(line for line in open(hmm_file_path)) if hmm_file_path else "NA"
    ALL_DATA["hmm_file"] = hmm_string

    ALL_DATA["taxonomy"] = full_taxonomy

    tax_function_this = dict()
    for c in tax_function:
        vals = np.append(tax_function[c].intercept_, tax_function[c].coef_)
        tax_function_this[c] = vals
    ALL_DATA["tax_function"] = tax_function_this

    classifiers_this = dict()
    for c in classifiers:
        if classifiers[c] != "no_negative_examples":
            vals = np.append(classifiers[c].intercept_, classifiers[c].coef_)
            classifiers_this[c] = vals
        else:
            classifiers_this[c] = "no_negative_examples"
    ALL_DATA["classifiers"] = classifiers_this

    all_LMNN_this = dict()
    for c in all_LMNN:
        all_LMNN_this[c] = all_LMNN[c].components_.T
    ALL_DATA["LMNN"] = all_LMNN_this

    ALL_DATA["thresholds_NN"] = thresholds_NN
    ALL_DATA["centroid_seq"] = centroid_seq
    ALL_DATA["species_to_tax"] = species_to_tax

    # save ------------------------
    with open(output, "wb") as f:
        pickle.dump(ALL_DATA, f)
