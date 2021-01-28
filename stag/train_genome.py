"""
Scripts that trains the database for the genome
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import sys
import time
import os
import tempfile
import shutil
import subprocess
import shlex
import errno
import h5py
import re
import tarfile

from stag.helpers import check_file_exists
from stag.classify import classify

# find the length of the alignments --------------------------------------------
def find_length_ali(gene_db, fasta_input, protein_fasta_input):
    #outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
    #os.chmod(outfile.name, 0o644)

    return classify(gene_db, fasta_input=fasta_input, protein_fasta_input=protein_fasta_input) #, save_ali_to_file=outfile.name)

    #CMD = "stag classify -d "+gene_db
    #CMD = CMD + " -i "+temp_fasta
    #CMD = CMD + " -p "+temp_fasta2
    #CMD = CMD + " -S "+outfile.name

    #split_CMD = shlex.split(CMD)
    #stag_CMD = subprocess.Popen(split_CMD)

    #return_code = stag_CMD.wait()
    #if return_code:
    #    sys.stderr.write("[E::align] Error. cannot find length of the alignments\n\n")
    #    sys.exit(1)

    #len_ali = len(next(open(outfile.name)).rstrip().split("\t")) - 1
    #os.remove(outfile.name)

    #return len_ali

#===============================================================================
#                                      MAIN
#===============================================================================
def train_genome(output, list_genes, gene_thresholds, threads, verbose, concat_stag_db):
    # temp file where to save the result ---------------------------------------
    outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(outfile.name, 0o644)

    # we need a file with the thresholds ---------------------------------------
    check_file_exists(gene_thresholds,isfasta = False)
    genes_threhold_file = list()
    o = open(gene_thresholds)
    for i in o:
        vals = i.rstrip().split("\t")
        genes_threhold_file.append(vals[0])
    o.close()
    for name in list_genes.split(","):
        if not name.split("/")[-1] in genes_threhold_file:
            sys.stderr.write("[E::main] Error: ")
            sys.stderr.write("gene "+name.split("/")[-1]+" is missing from the threshold file (-T)\n")
            sys.exit(1)


    # we create a tar.gz with all the genes ------------------------------------
    tar = tarfile.open(outfile.name, "w:gz")
    for name in list_genes.split(","):
        check_file_exists(name,isfasta = False)
        try:
            name_file = os.path.basename(name)
            if name_file == "threshold_file.tsv":
                sys.stderr.write("[E::main] Error: gene databases cannot have name 'threshold_file.tsv'. Please, choose another name.\n")
                sys.exit(1)
            if name_file == "hmm_lengths_file.tsv":
                sys.stderr.write("[E::main] Error: gene databases cannot have name 'hmm_lengths_file.tsv'. Please, choose another name.\n")
                sys.exit(1)
            if name_file == "concatenated_genes_STAG_database.HDF5":
                sys.stderr.write("[E::main] Error: gene databases cannot have name 'concatenated_genes_STAG_database.HDF5'. Please, choose another name.\n")
                sys.exit(1)
            if len(name_file.split("##")) > 1:
                sys.stderr.write("Error with: "+name_file+"\n")
                sys.stderr.write("[E::main] Error: gene databases cannot have in the name '##'. Please, choose another name.\n")
                sys.exit(1)
            tar.add(name, name_file)
        except:
            sys.stderr.write("[E::main] Error: when adding "+name+" to the database\n")
            sys.exit(1)

    # we add the file with the thresholds to the tar.gz
    tar.add(gene_thresholds, "threshold_file.tsv")


    # we need to find the length of the alignments -----------------------------
    len_f = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(len_f.name, 0o644)
    # temp fasta file
    temp_fasta = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_fasta.name, 0o644)
    temp_fasta.write(">test\nAAA\n")
    temp_fasta.flush()
    # protein
    temp_fasta2 = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_fasta2.name, 0o644)
    temp_fasta2.write(">test\nA\n")
    temp_fasta2.flush()
    for gene_db in list_genes.split(","):
        len_this = find_length_ali(gene_db,temp_fasta.name,temp_fasta2.name)
        len_f.write(os.path.basename(gene_db) + "\t" + str(len_this) + "\n")
        len_f.flush()

    os.remove(temp_fasta.name)
    # we add the file with the lengths to the tar.gz
    tar.add(len_f.name, "hmm_lengths_file.tsv")


    # add file with stag DB of the concatenated alis ---------------------------
    tar.add(concat_stag_db, "concatenated_genes_STAG_database.HDF5")


    # close tar file -----------------------------------------------------------
    tar.close()
    # close
    try:
        outfile.flush()
        os.fsync(outfile.fileno())
        outfile.close()
    except:
        sys.stderr.write("[E::main] Error: failed to save the result\n")
        sys.exit(1)
    try:
        #os.rename(outfile.name,output) # atomic operation
        shutil.move(outfile.name,output) #It is not atomic if the files are on different filsystems.
    except:
        sys.stderr.write("[E::main] Error: failed to save the resulting database\n")
        sys.stderr.write("[E::main] you can find the file here:\n"+outfile.name+"\n")
        sys.exit(1)
