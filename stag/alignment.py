import logging

import pandas as pd
import numpy as np

# Function to identify the rownames and number of columns in an alignment
def find_raw_names_ncol(file_name):
    gene_names = list()
    with open(file_name) as f:
        for line in f:
            gene_names.append(line[:line.find("\t")].replace("/", "-"))
        return gene_names, line.count("\t")

# function to load an alignment produced by the "align" option =================
# Input:
#  - a file created by "align"
# Output:
#  - a panda object
# as a note, numpy.loadtxt is way slower than pandas read.csv
# It works also on .gz files
def load_alignment_from_file_2(file_name, safe_mode=False):
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
                alignment.iloc[pos] = np.array([c == "1" for c in line.split("\t")[1:]])

    logging.info(f'   LOAD_AL: Number of genes: {len(list(alignment.index.values))}')
    alignment = alignment.drop_duplicates()
    logging.info(f'   LOAD_AL: Number of genes, after removing duplicates: {len(list(alignment.index.values))}')
    return alignment


class MultipleAlignment:
    def __init__(self, fn):
        self.alignment = None
        self.ncols = None
        self.npads = None
        self._load_alignment(fn)

    def _load_alignment(self, fn):
        with open(fn) as _in:
            self.ncols, self.npads = map(int, next(_in).strip().split("\t"))
            self.alignment = pd.read_csv(_in, header=None, index_col=0, sep="\t")
            logging.info(f'   LOAD_AL: Number of genes: {len(list(self.alignment.index.values))}')
            self.alignment = self.alignment.drop_duplicates()
            logging.info(f'   LOAD_AL: Number of genes, after removing duplicates: {len(list(self.alignment.index.values))}')

    def get_index(self):
        return list(alignment.index.values))


    def get_auxiliary_negative_examples(self, positive_examples, negative_examples, required_negatives, min_positives=5):
        positives = self.alignment.loc[positive_examples, : ].to_numpy()

        n_positives = len(positives)
        for i in range(n_positives, min_positives):
            positives = np.vstack((
                positives, 
                random.choice(range(n_positives))
            ))

        # find possible genes to augment the negative set
        clade_indices = set()
        negative_candidates = list(set(self.alignment.index.values).difference(positive_examples + negative_examples))
        if negative_candidates:
            # at the highest level, it is not possible to add negative examples
            check_negatives = alignment.loc[negative_candidates, : ].to_numpy()
            n_negatives = len(negative_candidates)

            # choose n random positive clades for comparison
            check_positives = positives[random.sample(range(len(positives)), min_positives), ]

            random_clades = list()
            for clade in range(min_positives):

                m_for_diff = np.tile(check_positives[clade,], (n_negatives, 1))

                # differences = np.sum(np.bitwise_xor(m_for_diff, check_negatives), axis=1)
                differences = np.array([
                    sum(count_bits(int(c)) for c in b)
                    for b in np.bitwise_xor(m_for_diff, check_negatives)
                ], dtype=np.int64)

                non_zero = np.sum(differences != 0)

                differences = np.where(differences == 0, np.nan, differences)

                corr_ord = np.argsort(differences)[:non_zero + 1]

                random_clades.append(list(corr_ord))

            for indices in zip(*random_clades):
                clade_indices.update(indices)
                if len(clade_indices) > required_negatives:
                    break

        return [negative_candidates[i] for i in clade_indices]


    def get_rows(self, rows):
        raise NotImplementedError()
        return self.alignment.loc[ rows, : ].to_numpy()






