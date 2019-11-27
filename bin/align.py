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
import tempfile

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

# ------------------------------------------------------------------------------
# main function as a generator
def align_generator(seq_file, hmm_file, use_cmalign, n_threads, verbose):
    """Align sequences and transform them into 1-hot encoding, ready for
       classification.
    Parameters
    ----------
     seq_file :    file with the nucleotide sequences [string]
     hmm_file :    file with the hmm model [string]
     use_cmalign : if True, we use cmalign. If false, we use hmmalign [bool]
     n_threads:    number of threads to use for cmalign (hmmalign can run only
                   on one thread) [string/int]
     verbose:      how much info to print [int]
    Returns
    -------
     Returns a generator with:
     'fasta_id\taligned_sequence'
    """

    # check that the tools are available
    if use_cmalign:
        if not is_tool("cmalign"):
            sys.stderr.write("[E::align] Error: cmalign is not in the path. Please install Infernal.\n")
            sys.exit(1)
    else:
        if not is_tool("hmmalign"):
            sys.stderr.write("[E::align] Error: hmmalign is not in the path. Please install HMMER3.\n")
            sys.exit(1)

    # prepare the command to run
    cmd = "hmmalign --outformat A2M "
    if use_cmalign:
        cmd = "cmalign --cpu "+str(n_threads)+" --outformat A2M "

    cmd = cmd + hmm_file +" "+ seq_file

    if verbose > 4:
        sys.stderr.write("Command used to align the sequences: "+cmd+"\n")

    # run the command
    CMD = shlex.split(cmd)
    align_cmd = subprocess.Popen(CMD,stdout=subprocess.PIPE,)

    # parse the result and return/save to file
    for line in merged_fasta(align_cmd.stdout):
        converted_line = convert_alignment(line,verbose)
        yield converted_line

    # check that hmmalign/cmalign finished correctly
    align_cmd.stdout.close()
    return_code = align_cmd.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. hmmalig/cmalign failed\n")
        sys.exit(1)

# ------------------------------------------------------------------------------
# main function
def align_file(seq_file, hmm_file, use_cmalign, n_threads, verbose, res_file):
    """Align sequences and transform them into 1-hot encoding, ready for
       classification.
    Parameters
    ----------
     seq_file :    file with the nucleotide sequences [string]
     hmm_file :    file with the hmm model [string]
     use_cmalign : if True, we use cmalign. If false, we use hmmalign [bool]
     n_threads:    number of threads to use for cmalign (hmmalign can run only
                   on one thread) [string/int]
     verbose:      how much info to print [int]
     res_file:     where to save the result.
    Returns
    -------
     It will save the aligned sequences to the specified file.
    """

    # check that the tools are available
    if use_cmalign:
        if not is_tool("cmalign"):
            sys.stderr.write("[E::align] Error: cmalign is not in the path. Please install Infernal.\n")
            sys.exit(1)
    else:
        if not is_tool("hmmalign"):
            sys.stderr.write("[E::align] Error: hmmalign is not in the path. Please install HMMER3.\n")
            sys.exit(1)

    # prepare the command to run
    cmd = "hmmalign --outformat A2M "
    if use_cmalign:
        cmd = "cmalign --cpu "+str(n_threads)+" --outformat A2M "

    cmd = cmd + hmm_file +" "+ seq_file

    if verbose > 4:
        sys.stderr.write("Command used to align the sequences: "+cmd+"\n")

    # run the command
    CMD = shlex.split(cmd)
    align_cmd = subprocess.Popen(CMD,stdout=subprocess.PIPE,)

    # open the temporary file where to save the result
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_file.name, 0o644)

    # parse the result and save to temp_file
    for line in merged_fasta(align_cmd.stdout):
        converted_line = convert_alignment(line,verbose)
        temp_file.write(converted_line+"\n")

    # if we save the result to a file, then we close it now
    try:
        temp_file.flush()
        os.fsync(temp_file.fileno())
        temp_file.close()
    except:
        if verbose>4: sys.stderr.write("[E::align] Error when saving the resulting file\n")
        sys.exit(1)
    # move temp file to the final destination
    try:
        #os.rename(bam_temp_file.name,args.profile_bam_file) # atomic operation
        shutil.move(temp_file.name,res_file) #It is not atomic if the files are on different filsystems.
    except:
        sys.stderr.write("[E::align] The resulting file couldn't be save in the final destination. You can find the file here:\n"+temp_file.name+"\n")
        sys.exit(1)

    # check that hmmalign/cmalign finished correctly
    align_cmd.stdout.close()
    return_code = align_cmd.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. hmmalig/cmalign failed\n")
        sys.exit(1)
