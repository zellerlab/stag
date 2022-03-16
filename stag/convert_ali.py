"""
Scripts to convert between alignments
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import shutil
import os
import sys
import tempfile

from stag.helpers import read_fasta


# ------------------------------------------------------------------------------
# function to convert the nucleotide alignment into 1-hot encoding.
# Note that we select only the nucleotides that corresponds to the inner state
# of the HMM.
encoding_dic = {
    "A": "0\t0\t0\t0\t1",
    "C": "0\t0\t0\t1\t0",
    "G": "0\t0\t1\t0\t0",
    "T": "0\t1\t0\t0\t0",
    "U": "0\t1\t0\t0\t0",
    "others": "1\t0\t0\t0\t0"
}

decoding_dic = {
    '00001': "A",
    '00010': "C",
    '00100': "G",
    '01000': "T",
    '10000': "-",
}


def convert_alignment(fasta_record, verbose):
    n_aligned_characters = 0
    converted_ali = [fasta_record[0]]
    for character in fasta_record[1]:
        # 1-hot encoding
        # the ACGTU are converted, everything else that is upper case, is considered
        # as a gap ('-').
        # for example also 'N' is converted to "-" -> "1,0,0,0,0"
        # Note that the upper case letters and "-" represents alignment to the
        # hidden state of the HMM.
        
        if not character.islower():
            five_vals = encoding_dic.get(character)
            if five_vals is None:
                five_vals = encoding_dic["others"]
            else:
                n_aligned_characters += 1
            converted_ali.append(five_vals)
    
    return "\t".join(converted_ali), n_aligned_characters


# function that identify the type of input -------------------------------------
# either 1-hot encoding or fasta
# return either "1-hot" or "fasta"
def find_input_type(file_in, verbose):
    o = open(file_in, "r")
    _ = o.readline()
    line2 = o.readline()
    o.close()

    last_val = line2.rstrip().split("\t")[-1]

    return "1-hot" if last_val in "01" else "fasta"


# function to convert to 1-hot encoding ----------------------------------------
# given an aligned fasta input
def convert_to_1_hot(file_in, file_out, verbose):
    # open the temporary file where to save the result
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_file.name, 0o644)
    # go through the input file
    with open(file_in, "rb") as _in:
        for line in read_fasta(_in, head_start=1):
            converted_line, _ = convert_alignment(line, verbose)
            temp_file.write(converted_line + "\n")

    # we saved the result to a temp file, then we close it now
    try:
        temp_file.flush()
        os.fsync(temp_file.fileno())
        temp_file.close()
    except Exception:
        if verbose > 4:
            sys.stderr.write("[E::align] Error when saving the resulting file\n")
        sys.exit(1)
    # move temp file to the final destination
    try:
        shutil.move(temp_file.name, file_out)  # It is not atomic if the files are on different filsystems.
    except Exception:
        sys.stderr.write(
            "[E::align] The resulting file couldn't be save in the final destination. "
            f"You can find the file here:\n{temp_file.name}\n"
        )
        sys.exit(1)


# function to convert 1 single line of 1-hot encoding to a fasta ---------------
def back_to_fasta(line):
    # split by "\t"
    vals = line.rstrip().split("\t")
    # first value is the fasta id
    fasta_res = ">" + vals[0] + "\n"
    # now we parse the values
    vals = vals[1:]
    for pos in range(0, len(vals), 5):
        this_val = "".join(vals[pos:pos+5])
        fasta_res = fasta_res + decoding_dic[this_val]
    fasta_res = fasta_res + "\n"
    return fasta_res


# function to convert back to fasta ali ----------------------------------------
# given a 1-hot encoding input
def convert_to_fasta(file_in, file_out, verbose):
    # open the temporary file where to save the result
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_file.name, 0o644)
    # go through the input file
    o = open(file_in, "r")
    for line in o:
        fasta_line = back_to_fasta(line)
        temp_file.write(fasta_line)
    o.close()

    # we saved the result to a temp file, then we close it now
    try:
        temp_file.flush()
        os.fsync(temp_file.fileno())
        temp_file.close()
    except Exception:
        if verbose > 4:
            sys.stderr.write("[E::align] Error when saving the resulting file\n")
        sys.exit(1)
    # move temp file to the final destination
    try:
        shutil.move(temp_file.name, file_out)  # It is not atomic if the files are on different filsystems.
    except Exception:
        sys.stderr.write(
            "[E::align] The resulting file couldn't be save in the final destination. "
            f"You can find the file here:\n{temp_file.name}\n"
        )
        sys.exit(1)


def convert_ali(file_in, file_out, verbose):
    # First, we need to understand the type of input
    input_type = find_input_type(file_in, verbose)
    # Second we convert
    if input_type == "fasta":
        convert_to_1_hot(file_in, file_out, verbose)
    else:
        convert_to_fasta(file_in, file_out, verbose)
