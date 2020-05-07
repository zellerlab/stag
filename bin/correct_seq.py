"""
Scripts to correct sequences that are in the wrong orientation
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
# main function
def correct_seq(seq_file, hmm_file, use_cmalign, n_threads, verbose, min_perc_state, res_file):
    """Correct sequences orientation
    Parameters
    ----------
     seq_file:       file with the nucleotide sequences [string]
     hmm_file:       file with the hmm model [string]
     use_cmalign:    if True, we use cmalign. If false, we use hmmalign [bool]
     n_threads:      number of threads to use for cmalign (hmmalign can run only
                     on one thread) [string/int]
     verbose:        how much info to print [int]
     min_perc_state  internal states coverage percentge threshold for the sequences
     res_file:       where to save the result.
    Returns
    -------
     It will print to res_file the correct sequences
    """

    # we create a file with the reverse complement

    # check that the tools are available
    if use_cmalign:
        if not is_tool("cmalign"):
            sys.stderr.write("[E::align] Error: cmalign is not in the path. Please install Infernal.\n")
            sys.exit(1)
    else:
        if not is_tool("hmmalign"):
            sys.stderr.write("[E::align] Error: hmmalign is not in the path. Please install HMMER3.\n")
            sys.exit(1)

    if not is_tool("esl-reformat"):
        sys.stderr.write("[E::align] Error: esl-reformat is not in the path. Please install Easel.\n")
        sys.exit(1)

    # prepare the command to run
    cmd = "hmmalign "
    if use_cmalign:
        cmd = "cmalign --cpu "+str(n_threads)+" "

    # add the file
    cmd = cmd + hmm_file +" "+ seq_file

    if verbose > 4:
        sys.stderr.write("Command used to align the sequences: "+cmd+"\n")
