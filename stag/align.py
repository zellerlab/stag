"""
Scripts to align sequences and transoform them into 1-hot encoding
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import shutil
import subprocess
import shlex
import os
import sys
import tempfile
import re

import numpy as np

from stag.helpers import is_tool, read_fasta


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
    converted_ali = []
    for character in alignment:
        # 1-hot encoding
        # the ACGTU are converted, everything else that is upper case, is considered
        # as a gap ('-').
        # for example also 'N' is converted to "-" -> "1,0,0,0,0"
        # Note that the upper case letters and "-" represents alignment to the
        # hidden state of the HMM.
        if not character.islower():
            n_char += 1
            encoded_block = encoding_dic.get(character)
            if encoded_block:
                # others' high bit = 1
                n_aligned_characters += 1
            else:
                encoded_block = encoding_dic["others"]
            converted_ali += encoded_block

    return np.array(converted_ali, dtype=bool), n_aligned_characters / n_char * 100


# function that transform a protein MSA to a nucleotide MSA --------------------
# if check_length is True, then we check that
# len(protein) == len(gene)*3 OR len(protein)-3 == len(gene)*3
def protein2gene_alignment(gene_id, protein_alignment, gene_sequence, check_length=True):
    # check that the length is correct
    only_AA_from_ali = re.sub(r'\-', '', protein_alignment)
    if check_length:
        expected_gene_length = len(only_AA_from_ali) * 3
        # check if lengths of gene and protein sequence match, with or without stop codon
        if len(gene_sequence) != expected_gene_length and len(gene_sequence) - 3 != expected_gene_length:
            sys.stderr.write("Error, length of genes/alignment is not correct")
            sys.stderr.write(" (protein: "+str(len(only_AA_from_ali)*3)+", gene: "+str(len(gene_sequence))+")\n")
            sys.stderr.write(" ID: "+gene_id+"\n")
            return None

    # convert alignment
    pos_gene, al_gene = 0, []
    for res in protein_alignment:
        found = False
        if res == "-":
            al_gene.append("---")
            found = True
        elif res.isupper():
            al_gene += gene_sequence[pos_gene:pos_gene + 3]
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
        raise ValueError("[E::align] Error: cmalign is not in the path. Please install Infernal.")
    elif not is_tool("hmmalign"):
        raise ValueError("[E::align] Error: hmmalign is not in the path. Please install HMMER3.")
    if not is_tool("esl-reformat"):
        raise ValueError("[E::align] Error: esl-reformat is not in the path. Please install Easel.")

    aligner = f"cmalign --cpu {n_threads}" if use_cmalign else "hmmalign"
    seq_input = protein_file if protein_file else seq_file
    align_cmd = f"{aligner} {hmm_file} {seq_input}"

    if verbose > 4:
        print(f"Command used to align the sequences: {align_cmd}", file=sys.stderr)

    # run the command
    CMD = shlex.split(align_cmd)
    align_cmd = subprocess.Popen(CMD, stdout=subprocess.PIPE,)

    # command to parse the alignment from STOCKHOLM to fasta format
    cmd2 = "esl-reformat a2m -"
    CMD2 = shlex.split(cmd2)
    parse_cmd = subprocess.Popen(CMD2, stdin=align_cmd.stdout, stdout=subprocess.PIPE,)

    if protein_file:
        seq_stream = zip(read_fasta(parse_cmd.stdout, head_start=1),
                         read_fasta(seq_file, is_binary=False, head_start=1))
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
        raise ValueError("[E::align] Error. hmmalign/cmalign failed.")
    # check that converting the file worked correctly
    parse_cmd.stdout.close()
    return_code = parse_cmd.wait()
    if return_code:
        raise ValueError("[E::align] Error. esl-reformat failed.")

    # print the number of sequences that were filtered
    if verbose > 3:
        print(f" Number of sequences that pass the filter: {n_pass}", file=sys.stderr)
        print(f" Number of sequences that do not pass the filter: {n_not_pass}", file=sys.stderr)


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

    temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_file.name, 0o644)
    encoder = None
    with temp_file:
        for gid, ali in align_generator(
            seq_file, protein_file, hmm_file, use_cmalign,
            n_threads, verbose, False, min_perc_state
        ):
            if encoder is None:
                encoder = AlignmentEncoder(ali)
                print(encoder.ncols, encoder.npads, sep="\t", file=temp_file)

            print(gid, *encoder.encode(ali), sep="\t", file=temp_file)

        try:
            temp_file.flush()
            os.fsync(temp_file.fileno())
        except Exception:
            raise ValueError("[E::align] Error when saving the resulting file.")

    try:
        shutil.move(temp_file.name, res_file)
    except Exception:
        raise ValueError(
            "[E::align] The resulting file couldn't be saved. "
            f"You can find the file here:\n{temp_file.name}."
        )


class AlignmentEncoder:
    def __init__(self, aln_row):
        self.ncols = len(aln_row)
        self.npads = (32 - self.ncols % 32) if self.ncols % 32 else 0

    def encode(self, aln_row):
        split_row = np.array_split(aln_row > 0, np.arange(32, len(aln_row), 32))
        split_row[-1] = np.append(split_row[-1], np.zeros(self.npads) > 0)
        encoded = np.sum(
            np.apply_along_axis(lambda x:(x * (1 << np.arange(0, 32))[::-1]), 1, split_row),  # noqa: E231
            axis=1
        )
        return encoded
