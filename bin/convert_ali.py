"""
Scripts to convert between alignments
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import shutil
import time
import subprocess
import shlex
import os
import errno
import sys
import tempfile
import numpy as np
import re

#===============================================================================
#                                 FUNCTIONS
#===============================================================================

# ------------------------------------------------------------------------------
# function to check if a specific tool exists
def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

# ------------------------------------------------------------------------------
# function to convert a fasta file with multiple lines into a one line separated
# by a "\t"
# Example:
# >test_fasta_header
# ATTGCGATTTCT
# CGGTATCGGTAT
# CGGTTA
### TO:
# >test_fasta_header\tATTGCGATTTCTCGGTATCGGTATCGGTTA
def merge_fasta(filein):
    # filein is a stream of data (from hmmalign)
    seq = ""
    for line_b in filein:
        line = line_b.decode("utf-8").rstrip()
        if line.startswith(">"):
            if seq != "":
                yield seq[1:] # we skip the 0 character, which is ">"
            seq = line+"\t"
        else:
            seq = seq + line
    # give back the last sequence
    if seq != "":
        yield seq[1:] # we skip the 0 character, which is ">"

# ------------------------------------------------------------------------------
# function to convert the nucleotide alignment into 1-hot encoding.
# Note that we select only the nucleotides that corresponds to the inner state
# of the HMM.
encoding_dic = {
               "A":"0\t0\t0\t0\t1",
               "C":"0\t0\t0\t1\t0",
               "G":"0\t0\t1\t0\t0",
               "T":"0\t1\t0\t0\t0",
               "U":"0\t1\t0\t0\t0",
               "others":"1\t0\t0\t0\t0"
               }
encoding_dic_numpy = {
               "A":[False,False,False,False,True],
               "C":[False,False,False,True,False],
               "G":[False,False,True,False,False],
               "T":[False,True,False,False,False],
               "U":[False,True,False,False,False],
               "others":[True,False,False,False,False]
               }

def convert_alignment(merged_fasta,verbose):
    n_aligned_characters = 0
    converted_ali = merged_fasta.split("\t")[0] # first value is the gene_id
    for character in merged_fasta.split("\t")[1]:
        # 1-hot encoding
        # the ACGTU are converted, everything else that is upper case, is considered
        # as a gap ('-').
        # for example also 'N' is converted to "-" -> "1,0,0,0,0"
        # Note that the upper case letters and "-" represents alignment to the
        # hidden state of the HMM.
        five_vals = ""
        if not character.islower():
            try:
                five_vals = encoding_dic[character]
                # if it doesnt break, we count it as an aligned character
                n_aligned_characters = n_aligned_characters + 1
            except KeyError:
                five_vals = encoding_dic["others"]
            five_vals = "\t"+five_vals
        # if it was lower case character, then five_vals = ""
        converted_ali = converted_ali + five_vals
    return converted_ali, n_aligned_characters

def convert_alignment_numpy(merged_fasta,verbose):
    n_aligned_characters = 0
    gene_id = merged_fasta.split("\t")[0]
    converted_ali = list()
    for character in merged_fasta.split("\t")[1]:
        # 1-hot encoding
        # the ACGTU are converted, everything else that is upper case, is considered
        # as a gap ('-').
        # for example also 'N' is converted to "-" -> "1,0,0,0,0"
        # Note that the upper case letters and "-" represents alignment to the
        # hidden state of the HMM.
        if not character.islower():
            try:
                converted_ali.extend(encoding_dic_numpy[character])
                # if it doesnt break, we count it as an aligned character
                n_aligned_characters = n_aligned_characters + 1
            except KeyError:
                converted_ali.extend(encoding_dic_numpy["others"])
    to_return = dict()
    to_return[gene_id] = np.array(converted_ali,dtype=bool)
    return to_return, n_aligned_characters

# function that read genes and return them as one line -------------------------
def yield_genes(seq_file):
    o = open(seq_file,"r")
    seq = ""
    for line in o:
        if line.startswith(">"):
            if seq != "":
                yield seq[1:] # we skip the 0 character, which is ">"
            seq = line.rstrip()+"\t"
        else:
            seq = seq + line.rstrip()
    o.close()
    # give back the last sequence
    if seq != "":
        yield seq[1:] # we skip the 0 character, which is ">"

# function that transform a protein MSA to a nucleotide MSA --------------------
# if check_length is True, then we check that
# len(protein) == len(gene)*3 OR len(protein)-3 == len(gene)*3
def proteinAl_2_geneAl(protein_alignment, gene_sequence, check_length):
    gene_id = gene_sequence.split("\t")[0]

    protein_alignment = protein_alignment.split("\t")[1]
    gene_sequence = gene_sequence.split("\t")[1]

    # check that the lenght is correct
    correct_length = True
    only_AA_from_ali = re.sub(r'\-', '', protein_alignment)
    if check_length:
        if len(only_AA_from_ali)*3 != len(gene_sequence)-3:
            # some proteins do not have the '*' at the end
            if len(only_AA_from_ali)*3 != len(gene_sequence):
                sys.stderr.write("Error, length of genes/alignment is not correct")
                sys.stderr.write(" (protein: "+str(len(only_AA_from_ali)*3)+", gene: "+str(len(gene_sequence))+")\n")
                sys.stderr.write(" ID: "+gene_id+"\n")
                correct_length = False
                return None

    # convert alignment
    pos_gene = 0
    al_gene = ""
    for i in protein_alignment:
        found = False
        if i == "-":
            al_gene = al_gene + "---"
            found = True
        if i.isupper():
            al_gene = al_gene + gene_sequence[pos_gene] + gene_sequence[pos_gene+1] + gene_sequence[pos_gene+2]
            pos_gene = pos_gene + 3
            found = True
        if i.islower():
            found = True
            # since we have to remove the lower case letter, we do not
            # add those to the alignment, but we anyway increase pos_gene
            pos_gene = pos_gene + 3
        if not found:
            sys.stderr.write("Error, character not identified\n")

    return gene_id + "\t" + al_gene

# function that identify the type of input -------------------------------------
# either 1-hot encoding or fasta
# return either "1-hot" or "fasta"
def find_input_type(file_in, verbose):
    o = open(file_in)
    line1 = o.readline()
    line2 = o.readline()
    o.close()

    last_val = line2.rstrip().split("\t")[-1]
    if last_val == "1" or  last_val == "0":
        return "1-hot"
    else:
        return "fasta"


# ------------------------------------------------------------------------------
# main function
def convert_ali(file_in, file_out, verbose):
    # First, we need to understand the type of input
    input_type = find_input_type(file_in, verbose)
    print(input_type)
