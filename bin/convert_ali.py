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
        line = line_b.rstrip()
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

# function that identify the type of input -------------------------------------
# either 1-hot encoding or fasta
# return either "1-hot" or "fasta"
def find_input_type(file_in, verbose):
    o = open(file_in,"r")
    line1 = o.readline()
    line2 = o.readline()
    o.close()

    last_val = line2.rstrip().split("\t")[-1]
    if last_val == "1" or  last_val == "0":
        return "1-hot"
    else:
        return "fasta"

# function to convert to 1-hot encoding ----------------------------------------
def convert_to_1_hot(file_in, file_out, verbose):
    # open the temporary file where to save the result
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_file.name, 0o644)
    # go through the input file
    o = open(file_in,"r")
    for line in merge_fasta(o):
        converted_line, n_aligned_characters = convert_alignment(line,verbose)
        temp_file.write(converted_line+"\n")
    o.close()

# ------------------------------------------------------------------------------
# main function
def convert_ali(file_in, file_out, verbose):
    # First, we need to understand the type of input
    input_type = find_input_type(file_in, verbose)
    # Second we convert
    if input_type == "fasta":
        convert_to_1_hot(file_in, file_out, verbose)
    else:
        convert_to_fasta(file_in, file_out, verbose)
