import random
import logging
import time

import pandas as pd
import numpy as np


def count_bits(n):
    # https://stackoverflow.com/questions/9829578/fast-way-of-counting-non-zero-bits-in-positive-integer
    n = (n & 0x5555555555555555) + ((n & 0xAAAAAAAAAAAAAAAA) >> 1)
    n = (n & 0x3333333333333333) + ((n & 0xCCCCCCCCCCCCCCCC) >> 2)
    n = (n & 0x0F0F0F0F0F0F0F0F) + ((n & 0xF0F0F0F0F0F0F0F0) >> 4)
    n = (n & 0x00FF00FF00FF00FF) + ((n & 0xFF00FF00FF00FF00) >> 8)
    n = (n & 0x0000FFFF0000FFFF) + ((n & 0xFFFF0000FFFF0000) >> 16)
    n = (n & 0x00000000FFFFFFFF) + ((n & 0xFFFFFFFF00000000) >> 32)  # This last & isn't strictly necessary.
    return n


class EncodedAlignment:
    # this has some unsafe methods (__init__, filter_alignment), but I don't have time right now.
    # "fast" decoding via
    # https://stackoverflow.com/questions/22227595/convert-integer-to-binary-array-with-suitable-padding
    # https://stackoverflow.com/questions/45604688/apply-function-on-each-row-row-wise-of-a-numpy-array
    # https://pypi.org/project/bitarray/ <- another possibility
    def __init__(self, fn=None, other_aln=None, ncols=None, npads=None, drop_duplicates=True):
        self.alignment = None
        self.ncols = None
        self.npads = None
        if fn:
            self._load_alignment(fn, drop_duplicates=drop_duplicates)
        elif other_aln is not None:
            self.alignment = other_aln
            self.ncols = ncols
            self.npads = npads

    def filter_alignment(self, rows):
        return EncodedAlignment(
            other_aln=self.alignment.loc[list(rows), :],
            ncols=self.ncols,
            npads=self.npads
        )

    def _load_alignment(self, fn, drop_duplicates=True):
        with open(fn) as _in:
            self.ncols, self.npads = map(int, next(_in).strip().split("\t"))
            self.alignment = pd.read_csv(_in, header=None, index_col=0, sep="\t")
            logging.info(f'   LOAD_AL: Number of genes: {len(list(self.alignment.index.values))}')
            if drop_duplicates:
                self.alignment = self.alignment.drop_duplicates()
                logging.info(
                    '   LOAD_AL: Number of genes, after removing duplicates: '
                    f'{len(list(self.alignment.index.values))}'
                )

    def get_index(self):
        return list(self.alignment.index.values)

    def get_additional_negative_examples(
        self, positive_examples, negative_examples, required_negatives, min_positives=5
    ):
        positives = self.alignment.loc[positive_examples, :].to_numpy()

        n_positives = len(positives)
        for i in range(n_positives, min_positives):
            positives = np.vstack((
                positives,
                positives[random.choice(range(n_positives)), ]
            ))

        # find possible genes to augment the negative set
        clade_indices = set()
        negative_candidates = list(set(self.alignment.index.values).difference(positive_examples + negative_examples))
        if negative_candidates:
            # at the highest level, it is not possible to add negative examples
            check_negatives = self.alignment.loc[negative_candidates, :].to_numpy()
            n_negatives = len(negative_candidates)

            # choose n random positive clades for comparison
            check_positives = positives[random.sample(range(len(positives)), min_positives), ]

            random_clades = []
            for clade in range(min_positives):

                m_for_diff = np.tile(check_positives[clade, ], (n_negatives, 1))

                differences = np.apply_along_axis(
                    lambda x:np.count_nonzero(x & (1 << np.arange(32))[:, None] != 0),  # noqa: E231
                    1,
                    np.bitwise_xor(m_for_diff, check_negatives)
                )

                # differences = np.sum(
                # 	np.apply_along_axis(
                # 		lambda x:(((x[:, None] & (1 << np.arange(32))[::-1]) > 0).flatten()[:-self.npads] == 1),
                # 		1,
                # 		np.bitwise_xor(m_for_diff, check_negatives)
                # 	),
                # 	axis=1
                # )

                # original, ~90s for nodes with missing negatives
                #  differences = np.sum(np.bitwise_xor(m_for_diff, check_negatives), axis=1)
                # differences = np.sum(np.vstack([  # ~250s for nodes with missing negatives
                # 	(((row[:,None] & (1 << np.arange(32))[::-1])) > 0).flatten()[:-self.npads] == 1
                # 	for row in np.bitwise_xor(m_for_diff, check_negatives)
                # ]))

                # differences = np.array([ # ~2000s(!) for nodes with missing negatives
                # 	sum(count_bits(int(c)) for c in b)
                # 	for b in np.bitwise_xor(m_for_diff, check_negatives)
                # ], dtype=np.int64)

                non_zero = np.sum(differences != 0)

                differences = np.where(differences == 0, np.nan, differences)

                corr_ord = np.argsort(differences)[:non_zero + 1]

                random_clades.append(list(corr_ord))

            for indices in zip(*random_clades):
                clade_indices.update(indices)
                if len(clade_indices) > required_negatives:
                    break

        return [negative_candidates[i] for i in clade_indices]

    def get_rows(self, node, rows):
        # logging.info(f"Unpacking: {node} {len(rows)} rows...")
        # t0 = time.time()

        # subset = self.alignment.loc[rows, : ].to_numpy()
        subset = np.apply_along_axis(
            lambda x: (((x[:, None] & (1 << np.arange(32))[::-1])) > 0).flatten()[:-self.npads],
            1,
            self.alignment.loc[rows, :].to_numpy()
        )
        # logging.info(f"Unpacked {node}: {len(rows)} rows in {time.time() - t0:.3f}s.")

        return subset

    def get_rows_old(self, node, rows):
        logging.info(f"Unpacking: {node} {len(rows)} rows...")
        t0 = time.time()

        alignment = []
        for row in self.alignment.loc[rows, :].to_numpy():
            alignment.append(
                (((row[:, None] & (1 << np.arange(32))[::-1])) > 0).flatten()[:-self.npads]
            )
        logging.info(f"Unpacked {node}: {len(rows)} rows in {time.time() - t0:.3f}s.")

        return np.vstack(alignment) if alignment else np.array([])
