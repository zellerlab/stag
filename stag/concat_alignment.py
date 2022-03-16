import argparse
import os
import sys

import numpy as np

from stag.alignment import EncodedAlignment
from stag.align import AlignmentEncoder



def main():

	ap = argparse.ArgumentParser()
	ap.add_argument("alignments", nargs="*")
	ap.add_argument("--suffix", type=str, default=".ali")

	args = ap.parse_args()


	alignments = {
		os.path.basename(aln).replace(args.suffix, ""): EncodedAlignment(aln, drop_duplicates=False)
		for aln in args.alignments
	}
	

	common_index = set()
	for aln in alignments.values():
		common_index.update(aln.alignment.index)

	# print(*alignments.items(), sep="\n")

	# print(len(common_index))
	
	encoder = None
	for gene in sorted(common_index):
		bin_row = list()
		for i, aln in enumerate(alignments.values()):
			try:
				aln_row = aln.get_rows(gene, (gene,))
			except KeyError:
				bin_row += [0] * aln.ncols
			else:
				bin_row += list(aln_row[0])
		bin_row = np.array(bin_row)
		if encoder is None:
			encoder = AlignmentEncoder(bin_row)
			print(encoder.ncols, encoder.npads, sep="\t")
		print(gene, *encoder.encode(bin_row), sep="\t")


if __name__ == "__main__":
	main()

