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

#===============================================================================
#                                      MAIN
#===============================================================================
def train_genome(output, list_genes, gene_thresholds, threads, verbose):
    # temp file where to save the result -----------------------------------
    outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(outfile.name, 0o644)

    # we need a file with the thresholds -----------------------------------
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

    # we create a tar.gz with all the genes --------------------------------
    tar = tarfile.open(outfile.name, "w:gz")
    for name in list_genes.split(","):
        check_file_exists(name,isfasta = False)
        try:
            name_file = os.path.basename(name)
            tar.add(name, name_file)
        except:
            sys.stderr.write("[E::main] Error: when adding "+name+" to the database\n")
            sys.exit(1)

    # we add the file with the thresholds to the tar.gz
    tar.add(gene_thresholds, os.path.basename(gene_thresholds))
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
