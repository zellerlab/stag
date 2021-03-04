import tempfile
import shutil
import os
import sys

# ==============================================================================
# UNZIP THE DATABASE
# ==============================================================================
def load_genome_DB(database, tool_version, verbose):
    dirpath = tempfile.mkdtemp()
    shutil.unpack_archive(database, dirpath,"gztar")
    list_files = [f for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, f))]
    # there must be a file with the name 'threshold_file.tsv'
    if not "threshold_file.tsv" in list_files:
        sys.stderr.write("[E::align] Error: threshold_file.tsv is missing.\n")
        sys.exit(1)
    if not "hmm_lengths_file.tsv" in list_files:
        sys.stderr.write("[E::align] Error: hmm_lengths_file.tsv is missing.\n")
        sys.exit(1)
    if not "concatenated_genes_STAG_database.HDF5" in list_files:
        sys.stderr.write("[E::align] Error: concatenated_genes_STAG_database.HDF5 is missing.\n")
        sys.exit(1)
    # we load the thresholds and gene order
    gene_order = list()
    gene_thresholds = dict()
    o = open(os.path.join(dirpath, "threshold_file.tsv"))
    for line in o:
        vals = line.rstrip().split("\t")
        gene_thresholds[vals[0]] = vals[1]
        gene_order.append(vals[0])
    o.close()

    # we load the gene legths
    ali_lengths = dict()
    o = open(os.path.join(dirpath, "hmm_lengths_file.tsv"))
    for line in o:
        vals = line.rstrip().split("\t")
        ali_lengths[vals[0]] = vals[1]
    o.close()


    # we remove the threshold file from the list of genes
    list_files.remove("threshold_file.tsv")
    list_files.remove("hmm_lengths_file.tsv")
    list_files.remove("concatenated_genes_STAG_database.HDF5")
    return list_files,dirpath,gene_thresholds,gene_order,ali_lengths,os.path.join(dirpath, "concatenated_genes_STAG_database.HDF5")

