"""
Scripts that checks the input files for the database
"""

import sys

# there are two input files to check:
# - the taxonomy file
# - the fasta file

# ------------------------------------------------------------------------------
# 1. check taxonomy
def check_taxonomy(tax_path):
    # 0. check that file exists
    try:
        o = open(tax_path,"r")
    except:
        sys.stderr.write("ERROR. Couldn't open taxonomy file\n")
        return True

    # 1. check the number of levels is consistent
    number_of_taxonomic_levels = len(o.readline().rstrip().split("\t"))
    o.seek(0)

    sys.stderr.write("Detected "+str(number_of_taxonomic_levels-1)+" taxonomic levels\n")
    if number_of_taxonomic_levels < 2:
        sys.stderr.write(" ERROR: We need at least one level (Like: 'gene_ID\\tlevel1\\tlevel2')\n")
        sys.stderr.write("        For first line we detected only gene id: '"+o.readline().rstrip()+"'\n")
        return True

    # variable to save all levels ids (for test 2)
    tax_ids = dict()
    for i in range(number_of_taxonomic_levels-1):
        tax_ids[i] = set()

    # variable to test that there is only one pssible parent for a given id (test 3)
    parent = dict()

    # variable with gene ids (test4)
    gene_ids = list()

    sys.stderr.write("Check number of taxonomy levels...")
    found_error1 = False
    for i in o:
        vals = i.rstrip().split("\t")
        if len(vals) != number_of_taxonomic_levels:
            sys.stderr.write("\n ERROR: Line with different number of tax levels ("+str(len(vals))+" instead of "+str(number_of_taxonomic_levels)+"): "+i)
            found_error1 = True
        else:
            # set up for test 4
            gene_ids.append(vals[0])
            # set up for test 2
            for l in range(number_of_taxonomic_levels-1):
                tax_ids[l].add(vals[l+1])
            # set up for test 3
            for l in range(number_of_taxonomic_levels-2):
                current_n = vals[l+2]
                parent_n = vals[l+1]
                if current_n not in parent:
                    parent[current_n] = set()
                parent[current_n].add(parent_n)
    o.close()
    if not found_error1:
        sys.stderr.write("               correct")

    # 2. check that the names are unique for one level (i.e. two levels do not
    #    have the same id)
    sys.stderr.write("\nCheck if the names are unique across levels...")
    found_error2 = False
    for i in range(number_of_taxonomic_levels-2):
        for j in range(i+1,number_of_taxonomic_levels-1):
            intersect_s = tax_ids[i].intersection(tax_ids[j])
            if len(intersect_s) != 0:
                for v in intersect_s:
                    sys.stderr.write("\n ERROR: '"+str(v)+"' is present in both level "+str(i)+" and "+str(j)+"\n")
                    found_error2 = True
    if not found_error2:
        sys.stderr.write("   correct")

    # 3. check that there is no “convergent evolution”.
    #    For example, the same genus cannot appear in two different families.
    sys.stderr.write("\nCheck if there are multiple parents...")
    found_error3 = False
    for c in parent:
        if len(parent[c]) > 1:
            sys.stderr.write("\n ERROR: Node '"+c+"' has multiple parents: "+str(parent[c]))
            found_error3 = True
    if not found_error3:
        sys.stderr.write("           correct")

    # 4. check the gene ids (that they are unique)
    sys.stderr.write("\nFound "+str(len(gene_ids))+" genes (lines)\n")
    found_error4 = False
    gene_ids_unique = set(gene_ids)
    if len(gene_ids_unique) != len(gene_ids):
        found_error4 = True
        sys.stderr.write(" ERROR: There are only "+str(len(gene_ids_unique))+" unique gene ids\n")

    # if there is any error, return True
    return (found_error1 or found_error2 or found_error3 or found_error4),gene_ids_unique

# ------------------------------------------------------------------------------
# 2. check sequences
def check_sequences(file_name):
    sys.stderr.write("Check that the sequences are in fasta format...")
    try:
        o = open(file_name,"r")
    except:
        sys.stderr.write("\n ERROR: cannot open file\n")
        return True

    # check that it is a fasta file
    try:
        if not(o.readline().startswith(">")):
            sys.stderr.write("\n ERROR: Not a fasta file\n")
            o.close()
            return True
    except:
        o.close()
        sys.stderr.write("\n ERROR: Not a fasta file\n")
        return True

    sys.stderr.write("  correct\n")
    return False # if we arrive here, there were no errors

# ------------------------------------------------------------------------------
# 3. check correspondence between fasta file and sequence file
def check_correspondence(file_name, gene_ids_from_tax):
    sys.stderr.write("Check correspondences between tax gene ids,\n                              and sequences...")
    try:
        o = open(file_name,"r")
    except:
        sys.stderr.write("\n ERROR: cannot open file\n")
        return True

    found_error = False
    try:
        for ll in o:
            if ll.startswith(">"):
                gene_id = ll.rstrip()[1:]
                if gene_id not in gene_ids_from_tax:
                    found_error = True
                    sys.stderr.write("\n ERROR: '"+gene_id+"' not in the taxonomy\n")
        o.close()
    except:
        o.close()
        sys.stderr.write("\n ERROR: Not a fasta file\n")
        return True

    if not found_error:
        sys.stderr.write("   correct\n")

    return found_error

#===============================================================================
#                                      MAIN
#===============================================================================

def check_input_files(seq_file, tax_file, hmm_file, cmalign):
    # 1. check taxonomy alone
    sys.stderr.write("------ CHECK TAXONOMY FILE:\n")
    found_error_tax, gene_ids = check_taxonomy(tax_file)

    # 2. check that the seq file is a proper fasta file
    sys.stderr.write("\n------ CHECK FASTA FILE:\n")
    found_error_seq = check_sequences(seq_file)

    # 3. check correspondences between tax and fasta file
    sys.stderr.write("\n------ CHECK CORRESPONDENCES:\n")
    found_error_corr = check_correspondence(seq_file, gene_ids)

    sys.stderr.write("\n")

    if (found_error_tax or found_error_seq):
        sys.exit(1)
