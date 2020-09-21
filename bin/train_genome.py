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

# position of the script -------------------------------------------------------
path_this = os.path.realpath(__file__)
path_array = path_this.split("/")
stag_path = "/".join(path_array[0:-2]) + "/stag"


# function that checks if a file exists ----------------------------------------
def check_file_exists(file_name, isfasta = False):
    try:
        o = open(file_name,"r")
        # if fasta file, then check that it starts with ">"
        if isfasta:
            if not(o.readline().startswith(">")):
                sys.stderr.write(f"{bco.Red}{bco.Bold}[E::main] Error: {bco.ResetAll}")
                sys.stderr.write("Not a fasta file: "+file_name+"\n")
                sys.stderr.write("          Fasta file is expected to start with '>'\n")
                o.close()
                sys.exit(1)
        o.close()
    except Exception as e:
        sys.stderr.write(f"{bco.Red}{bco.Bold}[E::main] Error: {bco.ResetAll}")
        sys.stderr.write("Cannot open file: "+file_name+"\n")
        sys.stderr.write(str(e)+"\n")
        sys.exit(1)

# find the length of the alignments --------------------------------------------
def find_length_ali(gene_db,temp_fasta,temp_fasta2):
    outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(outfile.name, 0o644)

    CMD = stag_path + " classify -d "+gene_db
    CMD = CMD + " -i "+temp_fasta
    CMD = CMD + " -p "+temp_fasta2
    CMD = CMD + " -S "+outfile.name

    split_CMD = shlex.split(CMD)
    stag_CMD = subprocess.Popen(split_CMD)

    return_code = stag_CMD.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. cannot find length of the alignments\n\n")
        sys.exit(1)

    o = open(outfile.name,"r")
    len_ali = len(o.readline().rstrip().split("\t")) - 1
    o.close()
    os.remove(outfile.name)

    return len_ali

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
            sys.stderr.write("gene "+name.split("/")[-1]+" is missing from the threshold file (-e)\n")
            sys.exit(1)


    # we create a tar.gz with all the genes ------------------------------------
    tar = tarfile.open(outfile.name, "w:gz")
    for name in list_genes.split(","):
        check_file_exists(name,isfasta = False)
        try:
            name_file = os.path.basename(name)
            if name_file == "threshold_file.tsv":
                sys.stderr.write("[E::main] Error: gene databases cannot have name 'threshold_file.tsv'. Please, choose anothe name.\n")
                sys.exit(1)
            if name_file == "hmm_lengths_file.tsv":
                sys.stderr.write("[E::main] Error: gene databases cannot have name 'hmm_lengths_file.tsv'. Please, choose anothe name.\n")
                sys.exit(1)
            if name_file == "concatenated_genes_STAG_database.HDF5":
                sys.stderr.write("[E::main] Error: gene databases cannot have name 'concatenated_genes_STAG_database.HDF5'. Please, choose anothe name.\n")
                sys.exit(1)
            if len(name_file.split("##")) > 1:
                sys.stderr.write("Error with: "+name_file+"\n")
                sys.stderr.write("[E::main] Error: gene databases cannot have in the name '##'. Please, choose anothe name.\n")
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
        len_f.write(name_file + "\t" + str(len_this) + "\n")
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
