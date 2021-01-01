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
    n_char = 0
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
            n_char = n_char + 1
            try:
                five_vals = encoding_dic[character]
                # if it doesnt break, we count it as an aligned character
                n_aligned_characters = n_aligned_characters + 1
            except KeyError:
                five_vals = encoding_dic["others"]
            five_vals = "\t"+five_vals
        # if it was lower case character, then five_vals = ""
        converted_ali = converted_ali + five_vals
    return converted_ali, n_aligned_characters/n_char

def convert_alignment_numpy(merged_fasta,verbose):
    n_aligned_characters = 0
    n_char = 0
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
            n_char = n_char + 1
            try:
                converted_ali.extend(encoding_dic_numpy[character])
                # if it doesnt break, we count it as an aligned character
                n_aligned_characters = n_aligned_characters + 1
            except KeyError:
                converted_ali.extend(encoding_dic_numpy["others"])
    to_return = dict()
    to_return[gene_id] = np.array(converted_ali,dtype=bool)
    return to_return, n_aligned_characters/n_char

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

# ------------------------------------------------------------------------------
# main function as a generator
def align_generator(seq_file, protein_file, hmm_file, use_cmalign, n_threads, verbose, return_numpy, min_perc_state):
    """Align sequences and transform them into 1-hot encoding, ready for
       classification.
    Parameters
    ----------
     seq_file:     file with the nucleotide sequences [string]
     protein_file:  file with the protein sequences [string or None]
     hmm_file:     file with the hmm model [string]
     use_cmalign:  if True, we use cmalign. If false, we use hmmalign [bool]
     n_threads:    number of threads to use for cmalign (hmmalign can run only
                   on one thread) [string/int]
     verbose:      how much info to print [int]
     return_numpy: True if you want to return a numpy array instead of a string
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
    if not is_tool("esl-reformat"):
        sys.stderr.write("[E::align] Error: esl-reformat is not in the path. Please install Easel.\n")
        sys.exit(1)

    # prepare the command to run
    cmd = "hmmalign "
    if use_cmalign:
        cmd = "cmalign --cpu "+str(n_threads)+" "

    if protein_file == None:
        cmd = cmd + hmm_file +" "+ seq_file
    else:
        cmd = cmd + hmm_file +" "+ protein_file

    if verbose > 4:
        sys.stderr.write("Command used to align the sequences: "+cmd+"\n")

    # run the command
    CMD = shlex.split(cmd)
    align_cmd = subprocess.Popen(CMD,stdout=subprocess.PIPE,)

    # command to parse the alignment from STOCKHOLM to fasta format
    cmd2 = "esl-reformat a2m -"
    CMD2 = shlex.split(cmd2)
    parse_cmd = subprocess.Popen(CMD2,stdin=align_cmd.stdout,stdout=subprocess.PIPE,)

    # parse the result and return/save to file - NORMAL ------------------------
    if protein_file == None:
        for line in merge_fasta(parse_cmd.stdout):
            if return_numpy:
                converted_line, perc_aligned_characters = convert_alignment_numpy(line,verbose)
            else:
                converted_line, perc_aligned_characters = convert_alignment(line,verbose)
            if perc_aligned_characters >= min_perc_state:
                yield converted_line

    # parse the result and return/save to file - WITH PROTEINS -----------------
    if protein_file != None:
        for protein_line, gene_line in zip(merge_fasta(parse_cmd.stdout), yield_genes(seq_file)):
            line = proteinAl_2_geneAl(protein_line, gene_line, True)
            if return_numpy:
                converted_line, perc_aligned_characters = convert_alignment_numpy(line,verbose)
            else:
                converted_line, perc_aligned_characters = convert_alignment(line,verbose)
            if perc_aligned_characters >= min_perc_state:
                yield converted_line



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

# ------------------------------------------------------------------------------
# main function
def align_file(seq_file, protein_file, hmm_file, use_cmalign, n_threads, verbose, res_file, min_perc_state):
    """Align sequences and transform them into 1-hot encoding, ready for
       classification.
    Parameters
    ----------
     seq_file :    file with the nucleotide sequences [string]
     protein_file:  file with the protein sequences [string or None]
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

    if not is_tool("esl-reformat"):
        sys.stderr.write("[E::align] Error: esl-reformat is not in the path. Please install Easel.\n")
        sys.exit(1)

    # prepare the command to run
    cmd = "hmmalign "
    if use_cmalign:
        cmd = "cmalign --cpu "+str(n_threads)+" "

    if protein_file == None:
        cmd = cmd + hmm_file +" "+ seq_file
    else:
        cmd = cmd + hmm_file +" "+ protein_file

    if verbose > 4:
        sys.stderr.write("Command used to align the sequences: "+cmd+"\n")

    # run the command
    CMD = shlex.split(cmd)
    align_cmd = subprocess.Popen(CMD,stdout=subprocess.PIPE,)

    # command to parse the alignment from STOCKHOLM to fasta format
    cmd2 = "esl-reformat a2m -"
    CMD2 = shlex.split(cmd2)
    parse_cmd = subprocess.Popen(CMD2,stdin=align_cmd.stdout,stdout=subprocess.PIPE,)

    # open the temporary file where to save the result
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_file.name, 0o644)

    # parse the result and return/save to file - NORMAL ------------------------
    if protein_file == None:
        for line in merge_fasta(parse_cmd.stdout):
            converted_line, perc_aligned_characters = convert_alignment(line,verbose)
            if perc_aligned_characters >= min_perc_state:
                temp_file.write(converted_line+"\n")

    # parse the result and return/save to file - WITH PROTEINS -----------------
    if protein_file != None:
        for protein_line, gene_line in zip(merge_fasta(parse_cmd.stdout), yield_genes(seq_file)):
            line = proteinAl_2_geneAl(protein_line, gene_line, True)
            converted_line, perc_aligned_characters = convert_alignment(line,verbose)
            if perc_aligned_characters >= min_perc_state:
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
    # check that converting the file worked correctly
    parse_cmd.stdout.close()
    return_code = parse_cmd.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. esl-reformat failed\n")
        sys.exit(1)
