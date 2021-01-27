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
# function that creates a file with reverse complement
def rev_complement(seq_file, verbose):
    # we use seqtk to reverse complement
    if not is_tool("seqtk"):
        sys.stderr.write("[E::align] Error: seqtk is not in the path. Please install seqtk.\n")
        sys.exit(1)
    # temp file
    if verbose > 2: sys.stderr.write("Create file with reverse complement...")
    rev_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    cmd = "seqtk seq -r "+seq_file
    if verbose > 4:
        sys.stderr.write("\nCommand used to reverse complement: "+cmd+" > "+rev_file.name+"\n")
    CMD = shlex.split(cmd)

    parse_cmd = subprocess.Popen(CMD,stdout=rev_file,)
    rev_file.flush()
    os.fsync(rev_file.fileno())
    rev_file.close()
    return_code = parse_cmd.wait()
    if return_code:
        sys.stderr.write("\n[E::align] Error. seqtk failed\n")
        sys.exit(1)
    if verbose > 2: sys.stderr.write("done\n")
    return rev_file.name


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
                yield seq
            seq = line+"\t"
        else:
            seq = seq + line
    # give back the last sequence
    if seq != "":
        yield seq

# function that calculate the number of internal states per sequence -----------
def calc_al(fasta_file, hmm_file, use_cmalign, n_threads, verbose):
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
    cmd = cmd + hmm_file +" "+ fasta_file

    if verbose > 4:
        sys.stderr.write("Command used to align the sequences: "+cmd+"\n")

    # we call the command
    CMD = shlex.split(cmd)
    align_cmd = subprocess.Popen(CMD,stdout=subprocess.PIPE,)

    # parse the alignment
    cmd2 = "esl-reformat a2m -"
    CMD2 = shlex.split(cmd2)
    parse_cmd = subprocess.Popen(CMD2,stdin=align_cmd.stdout,stdout=subprocess.PIPE,)

    all_lines = dict()
    for line in merge_fasta(parse_cmd.stdout):
        id = line.split("\t")[0]
        # calculate the number of internal state covered
        mat_i_s = 0 # internal states that match (even mismatch is counted I guess), they are upper case letters
        deletions = 0 # number of deletions (they are "-")
        insetions = 0 # insertions are lower case letters
        for i in line.split("\t")[1]:
            if i == "-":
                deletions = deletions + 1
            else:
                if i.isupper():
                    mat_i_s = mat_i_s + 1
                if i.islower():
                    insetions = insetions + 1

        all_lines[id] = ( mat_i_s/(mat_i_s+deletions) ) * 100

    # check that hmmalign/cmalign finished correctly
    align_cmd.stdout.close()
    return_code = align_cmd.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. hmmalig/cmalign failed\n")
        sys.exit(1)
    # check that converting the file worked correctly
    parse_cmd.stdout.close()
    return_code = parse_cmd.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. esl-reformat failed\n")
        sys.exit(1)

    return all_lines




# ------------------------------------------------------------------------------
# find the one that align the best and save the result
def save_best_seq(seq_al, rev_al, seq_file, rev_file, min_perc_state, output, verbose):
    # prepare output
    if not(output is None):
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(outfile.name, 0o644)
    else:
        outfile = sys.stdout

    if verbose > 2: sys.stderr.write("Select sequences...")
    removed_seq_count = 0 # seq lower than min_perc_state cutoff
    rotated_seq_count = 0 # sequenced that need to be reverse complement
    original_seq_count = 0

    # original sequences + count different types
    o = open(seq_file,"r")
    for l in o:
        l = l.rstrip()
        if l.startswith(">"):
            if seq_al[l] < min_perc_state and rev_al[l] < min_perc_state:
                removed_seq_count = removed_seq_count + 1
                print_this = False
            else:
                if seq_al[l] >= rev_al[l]:
                    original_seq_count = original_seq_count + 1
                    print_this = True
                else:
                    rotated_seq_count = rotated_seq_count + 1
                    print_this = False
        if print_this:
            outfile.write(l+"\n")
    o.close()

    # reversed sequences
    o = open(rev_file,"r")
    for l in o:
        l = l.rstrip()
        if l.startswith(">"):
            if rev_al[l] > seq_al[l] and rev_al[l] > min_perc_state:
                print_this = True
            else:
                print_this = False
        if print_this:
            outfile.write(l+"\n")
    o.close()

    # write info
    if verbose > 2: sys.stderr.write("done\n")
    if verbose > 2: sys.stderr.write("Sequences in correct orientation: "+str(original_seq_count)+"\n")
    if verbose > 2: sys.stderr.write("Reverse-complemented sequences: "+str(rotated_seq_count)+"\n")
    if verbose > 2: sys.stderr.write("Dropped sequences (below threshold): "+str(removed_seq_count)+"\n")

    # close file with result
    if not(output is None):
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
            sys.stderr.write("[E::main] Error: failed to save the profile\n")
            sys.stderr.write("[E::main] you can find the file here:\n"+outfile.name+"\n")
            sys.exit(1)

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

    # 1. create reverse complement file
    rev_file = rev_complement(seq_file, verbose)

    # 2. align both files, and find the percentage of internal states covered
    if verbose > 2: sys.stderr.write("Align the two files...")
    seq_al = calc_al(seq_file, hmm_file, use_cmalign, n_threads, verbose)
    rev_al = calc_al(rev_file, hmm_file,use_cmalign, n_threads, verbose)
    if verbose > 2: sys.stderr.write("done\n")

    # 3. find the one that align the best and save the result
    save_best_seq(seq_al, rev_al, seq_file, rev_file, min_perc_state, res_file, verbose)

    # 4. close and remove reverse complementa fasta file
    os.remove(rev_file)
