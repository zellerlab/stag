"""
Scripts that checks the input files for the database
"""

import sys
import os
import subprocess
import tempfile
import shlex

from stag.helpers import is_tool, linearise_fasta

# there are two input files to check:
# - the taxonomy file
# - the fasta file

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


# ------------------------------------------------------------------------------
# 1. check taxonomy
def check_taxonomy(tax_path):
    # 0. check that file exists
    try:
        o = open(tax_path,"r")
    except:
        sys.stderr.write(f"{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR.{bcolors.ENDC} ")
        sys.stderr.write("Couldn't open taxonomy file\n")
        return True

    # 1. check the number of levels is consistent
    first_line = o.readline().rstrip().split("\t")
    first_line = [first_line[0]]+first_line[1].split(";")
    o.seek(0)
    number_of_taxonomic_levels = len(first_line)

    sys.stderr.write("Detected "+str(number_of_taxonomic_levels-1)+" taxonomic levels\n")
    if number_of_taxonomic_levels < 2:
        sys.stderr.write(f"{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("We need at least one level (Like: 'gene_ID\\tlevel1;level2')\n")
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

    # full tax per gene
    full_taxonomy = dict()

    sys.stderr.write("Check number of taxonomy levels.......................")
    found_error1 = False
    for i in o:
        vals = i.rstrip().replace("/","-").split("\t")
        vals = [vals[0]]+vals[1].split(";")
        if len(vals) != number_of_taxonomic_levels:
            sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
            sys.stderr.write("Line with different number of tax levels ("+str(len(vals))+" instead of "+str(number_of_taxonomic_levels)+"): "+i)
            found_error1 = True
        else:
            # save taxonomy
            full_taxonomy[vals[0]] = vals[-1]
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
        sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}")

    # 2. check that the names are unique for one level (i.e. two levels do not
    #    have the same id)
    sys.stderr.write("\nCheck if the names are unique across levels...........")
    found_error2 = False
    for i in range(number_of_taxonomic_levels-2):
        for j in range(i+1,number_of_taxonomic_levels-1):
            intersect_s = tax_ids[i].intersection(tax_ids[j])
            if len(intersect_s) != 0:
                for v in intersect_s:
                    sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
                    sys.stderr.write("'"+str(v)+"' is present in both level "+str(i)+" and "+str(j)+"\n")
                    found_error2 = True
    if not found_error2:
        sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}")

    # 3. check that there is no “convergent evolution”.
    #    For example, the same genus cannot appear in two different families.
    sys.stderr.write("\nCheck if there are multiple parents...................")
    found_error3 = False
    for c in parent:
        if len(parent[c]) > 1:
            sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
            sys.stderr.write("Node '"+c+"' has multiple parents: "+str(parent[c]))
            found_error3 = True
    if not found_error3:
        sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}")

    # 4. check the gene ids (that they are unique)
    sys.stderr.write("\nFound "+str(len(gene_ids))+" genes (lines)\n")
    found_error4 = False
    gene_ids_unique = set(gene_ids)
    if len(gene_ids_unique) != len(gene_ids):
        found_error4 = True
        sys.stderr.write(f"{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("There are only "+str(len(gene_ids_unique))+" unique gene ids\n")

    # if there is any error, return True
    return (found_error1 or found_error2 or found_error3 or found_error4),gene_ids_unique,full_taxonomy

# ------------------------------------------------------------------------------
# 2. check sequences
def check_sequences(file_name):
    sys.stderr.write("Check that the sequences are in fasta format..........")
    try:
        o = open(file_name,"r")
    except:
        sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("cannot open file\n")
        return True

    # check that it is a fasta file
    try:
        if not(o.readline().startswith(">")):
            sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
            sys.stderr.write("Not a fasta file\n")
            o.close()
            return True
    except:
        o.close()
        sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("Not a fasta file\n")
        return True

    sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}\n")

    # check duplicates ---------------------------------------------------------
    duplicates_info = dict()
    o = open(file_name,"r")
    gene_id = ""
    n_genes = 0
    for i in o:
        if i.startswith(">"):
            if gene_id != "":
                if not seq in duplicates_info:
                    duplicates_info[seq] = list()
                duplicates_info[seq].append(gene_id)
            gene_id = i.rstrip()
            seq = ""
            n_genes = n_genes + 1
        else:
            seq = seq + i.rstrip()
    if not seq in duplicates_info:
        duplicates_info[seq] = list()
    duplicates_info[seq].append(gene_id)
    o.close()

    sys.stderr.write("Number of genes: "+str(n_genes) + "\n")
    sys.stderr.write("Number of unique genes: "+str(len(duplicates_info)) + "\n")

    return False, duplicates_info # if we arrive here, there were no errors (False means no error)

# ------------------------------------------------------------------------------
# 2.b if there is a protein file, then we check that it is correct and that it
#     maps correctly to the gene file
def check_protein_file(seq_file, protein_file):
    # general check of the protein file
    sys.stderr.write("Check that the protein sequences are in fasta format..")
    try:
        o = open(protein_file,"r")
    except:
        sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("cannot open file\n")
        return True
    # check that it is a fasta file
    try:
        if not(o.readline().startswith(">")):
            sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
            sys.stderr.write("Not a fasta file\n")
            o.close()
            return True
    except:
        o.close()
        sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("Not a fasta file\n")
        return True

    sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}\n")

    # find length of genes
    sys.stderr.write("Load gene file: ")
    gene_lengths = list()
    gene_ids = list()
    o = open(seq_file,"r")
    cont = -1
    for i in o:
        if i.startswith(">"):
            cont = cont + 1
            gene_lengths.append(0)
            gene_ids.append(i.rstrip())
        else:
            gene_lengths[cont] = gene_lengths[cont] + len(i.rstrip())
    o.close()
    sys.stderr.write("   found "+str(len(gene_lengths))+" genes\n")

    # find length of proteins
    sys.stderr.write("Load protein file: ")
    protein_lengths = list()
    protein_ids = list()
    o = open(protein_file,"r")
    cont = -1
    for i in o:
        if i.startswith(">"):
            cont = cont + 1
            protein_lengths.append(0)
            protein_ids.append(i.rstrip())
        else:
            protein_lengths[cont] = protein_lengths[cont] + len(i.rstrip())
    o.close()
    sys.stderr.write("found "+str(len(protein_lengths))+" proteins\n")

    # check number of sequences is the same:
    if len(gene_lengths) != len(protein_lengths):
        sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("different number of sequences\n")
        return True


    # check that the lengths make sense
    sys.stderr.write("Check the gene/protein match lengths..................")
    found_error = False
    count = -1
    for g_len, p_len in zip(gene_lengths, protein_lengths):
        count = count + 1
        if g_len != p_len*3:
            if (g_len-3) != p_len*3:
                found_error = True
                sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
                sys.stderr.write("different lengths for gene: "+gene_ids[cont]+"; protein: "+protein_ids[cont]+"\n")

    if not found_error:
        sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}\n")
    return found_error

# ------------------------------------------------------------------------------
# 3. check correspondence between fasta file and sequence file
def check_correspondence(file_name, gene_ids_from_tax, duplicates_info, full_taxonomy, warning_file_check_input):
    sys.stderr.write("Check correspondences of gene ids to the tax ids......")
    try:
        o = open(file_name,"r")
    except:
        sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("cannot open file\n")
        return True

    found_error = False
    try:
        for ll in o:
            if ll.startswith(">"):
                gene_id = ll.rstrip()[1:]
                if gene_id not in gene_ids_from_tax:
                    found_error = True
                    sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
                    sys.stderr.write("'"+gene_id+"' not in the taxonomy\n")
        o.close()
    except:
        o.close()
        sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("Not a fasta file\n")
        return True

    if not found_error:
        sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}\n")


    # check that genes with same sequence have the same taxonomy ---------------
    if warning_file_check_input != None:
        warning_f = open(warning_file_check_input,"w")
        warning_f.write("-- Check taxonomy of genes with same sequence --\n")

    sys.stderr.write("Check taxonomy of genes with same sequence............")
    found_error2 = False
    for i in duplicates_info:
        if len(duplicates_info[i]) > 1:
            species_0 = full_taxonomy[duplicates_info[i][0][1:]]
            for j in duplicates_info[i]:
                if full_taxonomy[j[1:]] != species_0:
                    found_error2 = True
                    if warning_file_check_input != None:
                        warning_f.write(str(duplicates_info[i])+"\n")
                    else:
                        sys.stderr.write(f"\n{bcolors.WARNING}{bcolors.BOLD}{bcolors.UNDERLINE}   WARNING:{bcolors.ENDC} ")
                        sys.stderr.write(str(duplicates_info[i])+"\n")

    if warning_file_check_input != None:
        warning_f.close()

    if not found_error2:
        sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}\n")
    else:
        sys.stderr.write(f"\n{bcolors.WARNING}{bcolors.BOLD}{bcolors.UNDERLINE} WARNING:{bcolors.ENDC} ")
        sys.stderr.write("Some genes have same sequence, but different taxonomy.\n")
    return (found_error or found_error2)


# ------------------------------------------------------------------------------
# 4. check that the tool is in the path and check the alignment
def check_tool(seq_file, hmm_file, use_cmalign):
    if use_cmalign:
        sys.stderr.write("Check that 'cmalign' is in the path...................")
        if not is_tool("cmalign"):
            sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
            sys.stderr.write("cmalign is not in the path. Please install Infernal.\n")
            return True
    else:
        sys.stderr.write("Check that 'hmmalign' is in the path..................")
        if not is_tool("hmmalign"):
            sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
            sys.stderr.write("hmmalign is not in the path. Please install HMMER3.\n")
            return True
    # if we arrive here, then the tool is in the path:
    sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}\n")

    # check esl-reformat
    sys.stderr.write("Check that 'esl-reformat' is in the path..............")
    if not is_tool("esl-reformat"):
        sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
        sys.stderr.write("esl-reformat is not in the path. Please install Easel.\n")
        return True
    # if we arrive here, then the tool is in the path:
    sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}\n")



    # check that the file is correct -------------------------------------------
    # we create a temporary file with the first tree fasta sequences:
    sys.stderr.write("Try to run alignment tool.............................")
    sys.stderr.flush()
    temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(temp_file.name, 0o644)

    o = open(seq_file,"r")
    count = 0
    for line in o:
        if line.startswith(">"):
            count = count + 1
        if count < 4:
            temp_file.write(line)

    try:
        temp_file.flush()
        os.fsync(temp_file.fileno())
        temp_file.close()
    except:
        if verbose>4:
            sys.stderr.write(f"\n{bcolors.FAIL}{bcolors.BOLD}{bcolors.UNDERLINE} ERROR:{bcolors.ENDC} ")
            sys.stderr.write("Error when saving the temp file\n")
        sys.exit(1)

    # we create the command to call
    cmd = "hmmalign "
    if use_cmalign:
        cmd = "cmalign "

    cmd = cmd + hmm_file +" "+ temp_file.name

    # we call the command
    CMD = shlex.split(cmd)
    align_cmd = subprocess.Popen(CMD,stdout=subprocess.PIPE,)

    # parse the alignment
    cmd2 = "esl-reformat a2m -"
    CMD2 = shlex.split(cmd2)
    parse_cmd = subprocess.Popen(CMD2,stdin=align_cmd.stdout,stdout=subprocess.PIPE,)

    all_lines = list()
    for line in linearise_fasta(parse_cmd.stdout, head_start=1):
        all_lines.append(line)

    align_cmd.stdout.close()
    return_code = align_cmd.wait()
    if return_code:
        os.remove(temp_file.name)
        return True
    # check that converting the file worked correctly
    parse_cmd.stdout.close()
    return_code = parse_cmd.wait()
    if return_code:
        os.remove(temp_file.name)
        return True

    # remove temporary file
    os.remove(temp_file.name)
    # if we arrive here, then the tool is in the path:
    sys.stderr.write(f"{bcolors.OKGREEN}{bcolors.BOLD}{bcolors.UNDERLINE}correct{bcolors.ENDC}\n")




    # check alignment quality --------------------------------------------------
    sys.stderr.write("\nCheck alignment quality:\n")

    # number of internal HMM states
    n_internal_states = 0
    for i in all_lines[0].split("\t")[1]:
        # gap (deletions) are counted
        if i == "-":
            n_internal_states = n_internal_states + 1
        else:
            # and capital letters
            if i.isupper():
                n_internal_states = n_internal_states + 1
    sys.stderr.write(" Internal states: "+str(n_internal_states)+"\n")

    count = 0
    for al in all_lines:
        count = count + 1
        sys.stderr.write("\n Sequence "+str(count)+":\n")
        # count occurences
        mat_i_s = 0 # internal states that match (even mismatch is counted I guess), they are upper case letters
        deletions = 0 # number of deletions (they are "-")
        insetions = 0 # insertions are lower case letters
        for i in all_lines[count-1].split("\t")[1]:
            if i == "-":
                deletions = deletions + 1
            else:
                if i.isupper():
                    mat_i_s = mat_i_s + 1
                if i.islower():
                    insetions = insetions + 1
        # print
        sys.stderr.write("   Internal states matches: "+str(mat_i_s)+" ("+str(round(mat_i_s/n_internal_states * 100))+"%)\n")
        sys.stderr.write("   Deletions: "+str(deletions)+" ("+str(round(deletions/n_internal_states * 100))+"%)\n")
        sys.stderr.write("   Insertions: "+str(insetions)+"\n")



#===============================================================================
#                                      MAIN
#===============================================================================

def check_input_files(seq_file, protein_file, tax_file, hmm_file, cmalign, warning_file_check_input):
    # 1. check taxonomy alone
    sys.stderr.write(f"{bcolors.OKBLUE}{bcolors.BOLD}------ CHECK TAXONOMY FILE:{bcolors.ENDC}\n")
    found_error_tax, gene_ids, full_taxonomy = check_taxonomy(tax_file)

    # 2. check that the seq file is a proper fasta file
    sys.stderr.write(f"{bcolors.OKBLUE}{bcolors.BOLD}------ CHECK FASTA FILE:{bcolors.ENDC}\n")
    found_error_seq, duplicates_info = check_sequences(seq_file)

    # 2.b check correspondence between protein and gene file
    found_error_prot = False
    if protein_file != None:
        sys.stderr.write(f"{bcolors.OKBLUE}{bcolors.BOLD}------ CHECK PROTEIN AND GENE FILE:{bcolors.ENDC}\n")
        found_error_prot = check_protein_file(seq_file, protein_file)

    # 3. check correspondences between tax and fasta file
    sys.stderr.write(f"{bcolors.OKBLUE}{bcolors.BOLD}------ CHECK CORRESPONDENCES:{bcolors.ENDC}\n")
    found_error_corr = check_correspondence(seq_file, gene_ids, duplicates_info, full_taxonomy, warning_file_check_input)

    # 4. test tool and alignment
    sys.stderr.write(f"{bcolors.OKBLUE}{bcolors.BOLD}------ CHECK TOOL:{bcolors.ENDC}\n")
    if protein_file != None:
        found_error_tool = check_tool(protein_file, hmm_file, cmalign)
    else:
        found_error_tool = check_tool(seq_file, hmm_file, cmalign)

    sys.stderr.write("\n")

    if (found_error_tax or found_error_seq or found_error_prot or found_error_corr or found_error_tool):
        sys.exit(1)
