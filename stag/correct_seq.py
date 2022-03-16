"""
Scripts to correct sequences that are in the wrong orientation
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import shutil
import subprocess
import shlex
import os
import sys
import tempfile

from stag.helpers import is_tool, read_fasta


# ------------------------------------------------------------------------------
# function that creates a file with reverse complement
def rev_complement(seq_file, verbose):
    # we use seqtk to reverse complement
    if not is_tool("seqtk"):
        sys.stderr.write("[E::align] Error: seqtk is not in the path. Please install seqtk.\n")
        sys.exit(1)
    # temp file
    if verbose > 2:
        sys.stderr.write("Create file with reverse complement...")
    rev_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    cmd = f"seqtk seq -r {seq_file}"
    if verbose > 4:
        sys.stderr.write(f"\nCommand used to reverse complement: {cmd} > {rev_file.name}\n")
    
    parse_cmd = subprocess.Popen(shlex.split(cmd), stdout=rev_file)
    rev_file.flush()
    os.fsync(rev_file.fileno())
    rev_file.close()
    return_code = parse_cmd.wait()
    if return_code:
        sys.stderr.write("\n[E::align] Error. seqtk failed\n")
        sys.exit(1)
    if verbose > 2:
        sys.stderr.write("done\n")
    return rev_file.name


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
        cmd = f"cmalign --cpu {n_threads} "

    cmd = f"{cmd} {hmm_file} {fasta_file}"

    if verbose > 4:
        sys.stderr.write(f"Command used to align the sequences: {cmd}\n")

    # we call the command
    CMD = shlex.split(cmd)
    align_cmd = subprocess.Popen(CMD, stdout=subprocess.PIPE)

    # parse the alignment
    cmd2 = "esl-reformat a2m -"
    CMD2 = shlex.split(cmd2)
    parse_cmd = subprocess.Popen(CMD2, stdin=align_cmd.stdout, stdout=subprocess.PIPE)

    all_lines = {}
    for sid, seq in read_fasta(parse_cmd.stdout, head_start=0):
        # calculate the number of internal state covered
        mat_i_s = 0  # internal states that match (even mismatch is counted I guess), they are upper case letters
        deletions = 0  # number of deletions (they are "-")
        insertions = 0  # insertions are lower case letters
        for c in seq:
            if c == "-":
                deletions += 1
            elif c.isupper():
                mat_i_s += 1
            elif c.islower():
                insertions += 1

        all_lines[sid] = (mat_i_s / (mat_i_s + deletions)) * 100

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

    return all_lines


# ------------------------------------------------------------------------------
# find the one that align the best and save the result
def save_best_seq(seq_al, rev_al, seq_file, rev_file, min_perc_state, output, verbose):
    # prepare output
    if output is not None:
        outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(outfile.name, 0o644)
    else:
        outfile = sys.stdout

    if verbose > 2:
        sys.stderr.write("Select sequences...")
    removed_seq_count = 0  # seq lower than min_perc_state cutoff
    rotated_seq_count = 0  # sequenced that need to be reverse complement
    original_seq_count = 0

    with outfile:
        # original sequences + count different types
        for sid, seq in read_fasta(seq_file, head_start=0, is_binary=False):
            if seq_al[sid] < min_perc_state and rev_al[sid] < min_perc_state:
                removed_seq_count += 1
            elif seq_al[sid] >= rev_al[sid]:
                original_seq_count += 1
                print(sid, seq, sep="\n", file=outfile)
            else:
                rotated_seq_count += 1

        # reversed sequences
        for sid, seq in read_fasta(rev_file, head_start=0, is_binary=False):
            if rev_al[sid] > seq_al[sid] and rev_al[sid] > min_perc_state:
                print(sid, seq, sep="\n", file=outfile)

        # write info
        if verbose > 2:
            sys.stderr.write("done\n")
            sys.stderr.write(f"Sequences in correct orientation: {original_seq_count}\n")
            sys.stderr.write(f"Reverse-complemented sequences: {rotated_seq_count}\n")
            sys.stderr.write(f"Dropped sequences (below threshold): {removed_seq_count}\n")

        if output is not None:
            try:
                outfile.flush()
                os.fsync(outfile.fileno())
            except Exception:
                sys.stderr.write("[E::main] Error: failed to save the result\n")
                sys.exit(1)

    # close file with result
    if output is not None:
        try:
            shutil.move(outfile.name, output)  # It is not atomic if the files are on different filesystems.
        except Exception:
            sys.stderr.write("[E::main] Error: failed to save the profile\n")
            sys.stderr.write(f"[E::main] you can find the file here:\n{outfile.name}\n")
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
    if verbose > 2:
        sys.stderr.write("Align the two files...")
    seq_al = calc_al(seq_file, hmm_file, use_cmalign, n_threads, verbose)
    rev_al = calc_al(rev_file, hmm_file, use_cmalign, n_threads, verbose)
    if verbose > 2:
        sys.stderr.write("done\n")

    # 3. find the one that align the best and save the result
    save_best_seq(seq_al, rev_al, seq_file, rev_file, min_perc_state, res_file, verbose)

    # 4. close and remove reverse complementa fasta file
    os.remove(rev_file)
