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

from stag.helpers import is_tool, read_fasta

#===============================================================================
#                                 FUNCTIONS
#===============================================================================

# ------------------------------------------------------------------------------
# function to convert the nucleotide alignment into 1-hot encoding.
# Note that we select only the nucleotides that corresponds to the inner state
# of the HMM.
encoding_dic = {
    "A": [0, 0, 0, 0, 1],
    "C": [0, 0, 0, 1, 0],
    "G": [0, 0, 1, 0, 0],
    "T": [0, 1, 0, 0, 0],
    "U": [0, 1, 0, 0, 0],
    "others": [1, 0, 0, 0, 0]
}

def convert_alignment(alignment, verbose, as_numpy=False):
    n_aligned_characters, n_char = 0, 0
    converted_ali = list()
    for character in alignment:
        # 1-hot encoding
        # the ACGTU are converted, everything else that is upper case, is considered
        # as a gap ('-').
        # for example also 'N' is converted to "-" -> "1,0,0,0,0"
        # Note that the upper case letters and "-" represents alignment to the
        # hidden state of the HMM.
        if not character.islower():
            n_char += 1
            encoded_block = encoding_dic.get(character) #, encoding_dic["others"])
            if encoded_block: #not encoded_block[0]:
                # others' high bit = 1
                n_aligned_characters += 1
            else:
                encoded_block = encoding_dic["others"]
            converted_ali.extend(encoded_block)
    #if as_numpy:
    #    converted_ali = np.array(list(map(bool, converted_ali)), dtype=bool)
    return np.array(converted_ali, dtype=bool), n_aligned_characters / n_char * 100

# function that transform a protein MSA to a nucleotide MSA --------------------
# if check_length is True, then we check that
# len(protein) == len(gene)*3 OR len(protein)-3 == len(gene)*3
def protein2gene_alignment(gene_id, protein_alignment, gene_sequence, check_length=True):
    # check that the length is correct
    only_AA_from_ali = re.sub(r'\-', '', protein_alignment)
    if check_length:
        expected_gene_length = len(only_AA_from_ali) * 3
        # check if lengths of gene and protein sequence match, with or without stop codon
        if len(gene_sequence) != expected_gene_length and len(gene_sequence) - 3 != expected_gene_length:
            sys.stderr.write("Error, length of genes/alignment is not correct")
            sys.stderr.write(" (protein: "+str(len(only_AA_from_ali)*3)+", gene: "+str(len(gene_sequence))+")\n")
            sys.stderr.write(" ID: "+gene_id+"\n")
            return None

    # convert alignment
    pos_gene, al_gene = 0, list()
    for res in protein_alignment:
        found = False
        if res == "-":
            al_gene.append("---")
            found = True
        elif res.isupper():
            al_gene.extend(gene_sequence[pos_gene:pos_gene + 3])
            pos_gene += 3
            found = True
        elif res.islower():
            found = True
            # since we have to remove the lower case letter, we do not
            # add those to the alignment, but we anyway increase pos_gene
            pos_gene += 3
        if not found:
            sys.stderr.write("Error, character not identified\n")

    return "".join(al_gene)

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
     (fasta_id, aligned_sequence) tuples
    """

    # number of sequences that pass and sont pass the filter
    n_pass, n_not_pass = 0, 0
    # check that the tools are available
    if use_cmalign and not is_tool("cmalign"):
        sys.stderr.write("[E::align] Error: cmalign is not in the path. Please install Infernal.\n")
        sys.exit(1)
    elif not is_tool("hmmalign"):
        sys.stderr.write("[E::align] Error: hmmalign is not in the path. Please install HMMER3.\n")
        sys.exit(1)
    if not is_tool("esl-reformat"):
        sys.stderr.write("[E::align] Error: esl-reformat is not in the path. Please install Easel.\n")
        sys.exit(1)

    # prepare the command to run
    cmd = "hmmalign "
    if use_cmalign:
        cmd = "cmalign --cpu "+str(n_threads)+" "

    if not protein_file:
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

    if protein_file:
        seq_stream = zip(read_fasta(parse_cmd.stdout, head_start=1),
                         read_fasta(open(seq_file), is_binary=False, head_start=1))
    else:
        seq_stream = read_fasta(parse_cmd.stdout, head_start=1)

    for item in seq_stream:
        if protein_file:
            (pid, pseq), (gid, gseq) = item
            if pid != gid:
                sys.stderr.write("[E::align] Error. protein and gene identifiers {} {} don't match.".format(pid, gid))
                sys.exit(1)
            gseq = protein2gene_alignment(gid, pseq, gseq, check_length=True)
        else:
            gid, gseq = item

        converted_ali, perc_aligned_characters = convert_alignment(gseq, verbose, as_numpy=return_numpy)
        if perc_aligned_characters >= min_perc_state:
            n_pass += 1
            yield gid, converted_ali
        else:
            n_not_pass += 1

    # check that hmmalign/cmalign finished correctly
    align_cmd.stdout.close()
    return_code = align_cmd.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. hmmalign/cmalign failed\n")
        sys.exit(1)
    # check that converting the file worked correctly
    parse_cmd.stdout.close()
    return_code = parse_cmd.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. esl-reformat failed\n")
        sys.exit(1)

    # print the number of sequences that were filtered
    if verbose > 3:
        sys.stderr.write(" Number of sequences that pass the filter: "+str(n_pass)+"\n")
        sys.stderr.write(" Number of sequences that do not pass the filter: "+str(n_not_pass)+"\n")

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

    # open the temporary file where to save the result
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_file.name, 0o644)
    with temp_file:
        for gid, ali in align_generator(seq_file, protein_file, hmm_file, use_cmalign,
                                        n_threads, verbose, False, min_perc_state):
            print(gid, *map(int, ali), sep="\t", file=temp_file)

        # if we save the result to a file, then we close it now
        try:
            temp_file.flush()
            os.fsync(temp_file.fileno())
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
