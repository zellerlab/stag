import logging

import pandas as pd
import numpy as np

# Function to identify the rownames and number of columns in an alignment
def find_raw_names_ncol(file_name):
    gene_names = list()
    with open(file_name) as f:
        for line in f:
            gene_names.append(line[:line.find("\t")])
        return gene_names, line.count("\t")

# function to load an alignment produced by the "align" option =================
# Input:
#  - a file created by "align"
# Output:
#  - a panda object
# as a note, numpy.loadtxt is way slower than pandas read.csv
# It works also on .gz files
def load_alignment_from_file(file_name, safe_mode=False):
    gene_names, align_length = find_raw_names_ncol(file_name)
    alignment = pd.DataFrame(False, index=gene_names, columns=range(align_length))
    with open(file_name) as f:
        if safe_mode:
            for pos, line in enumerate(f):
                align = [int(c) for c in line.split("\t")[1:]]
                if len(align) != align_length or any((c != 0 and c != 1) for c in align):
                    raise ValueError(f"Malformatted alignment in line {pos}:\n{gene}\t{''.join(align)}")
                alignment.iloc[pos] = np.array([c == 1 for c in align])
        else:
            for pos, line in enumerate(f):
                alignment.iloc[pos] = np.array([int(c) == 1 for c in line.split("\t")[1:]])

    logging.info(f'   LOAD_AL: Number of genes: {len(list(alignment.index.values))}')
    alignment = alignment.drop_duplicates()
    logging.info(f'   LOAD_AL: Number of genes, after removing duplicates: {len(list(alignment.index.values))}')
    return alignment
