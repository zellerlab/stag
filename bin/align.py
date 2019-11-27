"""
Scripts to align sequences and transoform them into 1-hot encoding
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import shutil
import time
import subprocess
import shlex
import os
import errno
import sys

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
def merged_fasta(filein):
    # filein is a stream of data (from hmmalign)
    seq = ""
    for line_b in filein:
        line = line_b.decode("utf-8").rstrip()
        if line.startswith(">"):
            if seq != "":
                yield seq
            seq = line+"\t"
        else:
            seq = seq + line
    # give back the last sequence
    if seq != "":
        yield seq

# ------------------------------------------------------------------------------
# function to convert the nucleotide alignment into 1-hot encoding.
# Note that we select only the nucleotides that corresponds to the inner state
# of the HMM.
encoding_dic = {
               "A":"00001",
               "C":"00010",
               "G":"00100",
               "T":"01000",
               "U":"01000",
               "others":"10000"
               }
def convert_alignment(merged_fasta,verbose):
    converted_ali = merged_fasta.split("\t")[0]+"\t"
    for i in merged_fasta.split("\t")[1]:
        # 1-hot encoding
        # the ACGTU are converted, everything else that is upper case, is considered
        # as a gap ('-').
        # for example also 'N' is converted to "-" -> "1,0,0,0,0"
        # Note that the upper case letters and "-" represents alignment to the
        # hidden state of the HMM.
        i_c = ""
        if not i.islower():
            if i in encoding_dic:
                i_c = encoding_dic[i]
            else:
                i_c = encoding_dic["others"]
        converted_ali = converted_ali + i_c
    return converted_ali
