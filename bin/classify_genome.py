"""
Scripts that classify a genome
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import sys
import time
import os
import tempfile
import shutil
import subprocess
import shlex
import errno
import h5py
import re

# ------------------------------------------------------------------------------
# function to check if a specific tool exists
def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

# ------------------------------------------------------------------------------
# dev null
try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

# ==============================================================================
# UNZIP THE DATABASE
# ==============================================================================
def load_genome_DB(database, tool_version, verbose):
    dirpath = tempfile.mkdtemp()
    shutil.unpack_archive(database, dirpath,"gztar")
    list_files = [f for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, f))]
    # there must be a file with the name 'threshold_file.tsv'
    if not "threshold_file.tsv" in list_files:
        sys.stderr.write("[E::align] Error: threshold_file.tsv is missing.\n")
        sys.exit(1)
    # we load the thresholds
    gene_thresholds = dict()
    o = open(os.path.join(dirpath, "threshold_file.tsv"))
    for line in o:
        vals = line.rstrip().split("\t")
        gene_thresholds[vals[0]] = vals[1]
    o.close()
    # we remove the threshold file from the list of genes
    list_files.remove("threshold_file.tsv")
    return list_files,dirpath,gene_thresholds

# ==============================================================================
# RUN PRODIGAL
# ==============================================================================
# run prodigal on one genome
def run_prodigal(genome, verbose):
    # we need two files, one for the proteins and one for the genes
    genes = tempfile.NamedTemporaryFile(delete=False, mode="w")
    proteins = tempfile.NamedTemporaryFile(delete=False, mode="w")
    # prodigal command
    prodigal_command = "prodigal -i "+genome+" -d "+genes.name
    prodigal_command = prodigal_command + " -a "+proteins.name
    # run prodigal
    if not is_tool("prodigal"):
        sys.stderr.write("[E::align] Error: prodigal is not in the path.\n")
        sys.exit(1)
    CMD2 = shlex.split(prodigal_command)
    parse_cmd = subprocess.Popen(CMD2, stdout=DEVNULL,stderr=subprocess.PIPE)
    # we save stderr if necessary
    all_stderr = ""
    for line in parse_cmd.stderr:
        line = line.decode('ascii')
        all_stderr = all_stderr + line
    return_code = parse_cmd.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. prodigal failed\n\n")
        sys.stderr.write(all_stderr)
        sys.exit(1)

    # we re-name the header of the fasta files ---------------------------------
    # we expect to have the same number of genes and proteins, and also that the
    #
    parsed_genes = tempfile.NamedTemporaryFile(delete=False, mode="w")
    parsed_proteins = tempfile.NamedTemporaryFile(delete=False, mode="w")

    o = open(genes.name,"r")
    count = 0
    for i in o:
        if i.startswith(">"):
            parsed_genes.write(">"+genome+"_"+str(count)+"\n")
            count = count + 1
        else:
            parsed_genes.write(i)
    o.close()

    o = open(proteins.name,"r")
    count = 0
    for i in o:
        if i.startswith(">"):
            parsed_proteins.write(">"+genome+"_"+str(count)+"\n")
            count = count + 1
        else:
            parsed_proteins.write(i)
    o.close()

    # remove old files
    os.remove(genes.name)
    os.remove(proteins.name)

    # close new files
    parsed_proteins.flush()
    os.fsync(parsed_proteins.fileno())
    parsed_proteins.close()
    parsed_genes.flush()
    os.fsync(parsed_genes.fileno())
    parsed_genes.close()

    return [parsed_genes.name, parsed_proteins.name]

# run prodigal on all the genomes listed in fasta_input
def run_prodigal_genomes(genomes_file_list, verbose):
    result = dict()
    for g in genomes_file_list:
        result[g] = run_prodigal(g, verbose)
    return result


# ==============================================================================
# EXTRACT THE MARKER GENES
# ==============================================================================
# find gene ids that we can use (run hmmsearch)
def extract_gene_from_one_genome(file_to_align, hmm_file, gene_threshold):
    # INFO: genes_path, proteins_path [where to save the result]
    # we run hmmsearch
    temp_hmm = tempfile.NamedTemporaryFile(delete=False, mode="w")
    hmm_cmd = "hmmsearch --tblout "+temp_hmm.name+" "+hmm_file+" "+file_to_align

    CMD = shlex.split(hmm_cmd)
    hmm_CMD = subprocess.Popen(CMD, stdout=DEVNULL,stderr=subprocess.PIPE)
    # we save stderr if necessary
    all_stderr = ""
    for line in hmm_CMD.stderr:
        line = line.decode('ascii')
        all_stderr = all_stderr + line
    return_code = hmm_CMD.wait()
    if return_code:
        sys.stderr.write("[E::align] Error. hmmsearch failed\n\n")
        sys.stderr.write(all_stderr)
        sys.exit(1)

    # in temp_hmm.name there is the result from hmm ----------------------------
    # we select which genes/proteins we need to extract from the fasta files
    # produced by prodigal
    sel_genes = dict()
    o = open(temp_hmm.name,"r")
    for line in o:
        if not line.startswith("#"):
            vals = re.sub(" +"," ",line.rstrip()).split(" ")
            gene_id = vals[0]
            e_val = vals[4]
            score = vals[5]
            if float(score) > float(gene_threshold):
                sel_genes[gene_id] = score
    o.close()

    # remove file with the result from the hmm
    if os.path.isfile(temp_hmm.name): os.remove(temp_hmm.name)

    # return
    return sel_genes



# for one marker gene, we extract all the genes/proteins from all genomes
def extract_genes(mg_name, hmm_file, use_protein_file, genomes_pred, gene_threshold, all_genes_raw):
    # we go throught the genome and find the genes that pass the filter
    genes_pass_filter = dict()
    for g in genomes_pred:
        if not (g in all_genes_raw):
            all_genes_raw[g] = dict()
        if mg_name in all_genes_raw[g]:
            sys.stderr.write("Error. gene already present\n")
        # which file do we use for the hmmsearch?
        if use_protein_file:
            file_to_align = genomes_pred[g][1]
        else:
            file_to_align = genomes_pred[g][0]
        # call function that uses hmmsearch
        all_genes_raw[g][mg_name] = extract_gene_from_one_genome(file_to_align, hmm_file, gene_threshold)

def select_genes(all_genes_raw, keep_all_genes):
    return_dict = dict()
    # all_genes_raw: genome1: MG1: geneA: 276
    #                              geneB: 243
    #                         MG2: geneC: 589
    #                genome2: MG1: geneX: 267
    #                              geneY: 212
    #                         MG2: geneZ: 459
    #                              geneY: 543
    for genome in all_genes_raw:
        return_dict[genome] = dict()
        # we first check if there is any gene that is in multiple mgs:
        gene_sel = dict()
        for mg in all_genes_raw[genome]:
            for g in all_genes_raw[genome][mg]:
                if not (g in gene_sel):
                    gene_sel[g] = float(all_genes_raw[genome][mg][g])
                else:
                    if float(all_genes_raw[genome][mg][g]) > float(gene_sel[g]):
                        gene_sel[g] = float(all_genes_raw[genome][mg][g])
        # in gene_sel there is the gene id -> highest score
        # example: geneX->267; geneY->543; geneZ->459

        # now we select the correct genes and decide if keep one or many
        for mg in all_genes_raw[genome]:
            return_dict[genome][mg] = list()
            # if we keep all genes
            if keep_all_genes:
                for g in all_genes_raw[genome][mg]:
                    if float(all_genes_raw[genome][mg][g]) == float(gene_sel[g]):
                        return_dict[genome][mg].append(g)
            # if we keep only one gene per marker gene
            if not keep_all_genes:
                max_v = 0
                sel_gene = ""
                for g in all_genes_raw[genome][mg]:
                    if float(all_genes_raw[genome][mg][g]) == float(gene_sel[g]):
                        if float(all_genes_raw[genome][mg][g]) > max_v:
                            max_v = float(all_genes_raw[genome][mg][g])
                            sel_gene = g
                if max_v != 0:
                    return_dict[genome][mg].append(sel_gene)
    return return_dict

# function that extract the genes and proteins based on the IDs from
# "selected_genes", only for one marker gene
def extract_genes_from_fasta(mg, selected_genes, genomes_pred, verbose, use_protein_file):
    genes = tempfile.NamedTemporaryFile(delete=False, mode="w")
    n_genes = 0
    if use_protein_file:
        proteins = tempfile.NamedTemporaryFile(delete=False, mode="w")
        n_proteins = 0

    for genome in selected_genes:
        if not(mg in selected_genes[genome]):
            sys.stderr.write("Warning: missing marker gene in genome "+genome+"\n")
        else:
            # for genes
            o = open(genomes_pred[genome][0])
            print_this = False
            for i in o:
                if i.startswith(">"):
                    if i[1:].rstrip() in selected_genes[genome][mg]:
                        print_this = True
                        n_genes = n_genes + 1
                        # we print a different header
                        genes.write(i.rstrip()+"##"+mg+"\n")
                    else:
                        print_this = False
                else:
                    if print_this:
                        genes.write(i)
            o.close()
            # for proteins
            if use_protein_file:
                o = open(genomes_pred[genome][1])
                print_this = False
                for i in o:
                    if i.startswith(">"):
                        if i[1:].rstrip() in selected_genes[genome][mg]:
                            print_this = True
                            n_proteins = n_proteins + 1
                            # we print a different header
                            proteins.write(i.rstrip()+"##"+mg+"\n")
                        else:
                            print_this = False
                    else:
                        if print_this:
                            proteins.write(i)
                o.close()

    if verbose > 3:
        sys.stderr.write(" Found "+str(n_genes)+" genes\n")
    if use_protein_file:
        if n_genes != n_proteins:
            sys.stderr.write("Error: Number of genes and proteins is different")

    # if there are no genes, we remove the files and return None
    if n_genes == 0:
        os.remove(genes.name)
        gene_file_name = None
        protein_file_name = None
        if use_protein_file:
            os.remove(proteins.name)

    else:
        gene_file_name = genes.name
        if use_protein_file:
            protein_file_name = proteins.name
        else:
            protein_file_name = "no_protein"

    return gene_file_name, protein_file_name


# extract the marker genes from the genes/proteins produced from prodigal
# for multiple genomes and multiple MGs
def fetch_MGs(database_files, database_path, genomes_pred, keep_all_genes, gene_thresholds, verbose):
    all_genes_raw = dict()
    mg_info_use_protein = dict()
    for mg in database_files:
        # for each MG, we extract the hmm and if using proteins or not ---------
        path_mg = os.path.join(database_path, mg)
        f = h5py.File(path_mg, 'r')
        # first, we save a temporary file with the hmm file
        hmm_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(hmm_file.name, 0o644)
        hmm_file.write(f['hmm_file'][0])
        hmm_file.flush()
        os.fsync(hmm_file.fileno())
        hmm_file.close()
        # second, check if we need to use proteins
        use_protein_file = False
        if f['align_protein'][0]:
            use_protein_file = True
            mg_info_use_protein[mg] = True
        else:
            mg_info_use_protein[mg] = False
        f.close()

        # run hmmsearch for each genome and find which genes pass the filter
        extract_genes(mg, hmm_file.name, use_protein_file, genomes_pred, gene_thresholds[mg], all_genes_raw)
        # the result is saved in all_genes_raw (passed as input)

        # remove hmm file
        os.remove(hmm_file.name)

    # now we need to select the marker genes and extracted them from the prodigal
    # fasta files
    # all_genes_raw: genome1: MG1: geneA: 276
    #                              geneB: 243
    #                         MG2: geneC: 589
    #                genome2: MG1: geneX: 267
    #                              geneY: 212
    #                         MG2: geneZ: 459
    #                              geneY: 543
    # NOTE: the same gene can apper in two different marker genes. Hence, we
    # need to assign it to only one
    # For example 'geneY' is in both MG1 and MG2 for genome 2.
    # we will select MG2 because the score is higher
    selected_genes = select_genes(all_genes_raw, keep_all_genes)
    #                               keep_all_genes=F    keep_all_genes=T
    # selected_genes: genome1: MG1:          (geneA)       (geneA,geneB)
    #                          MG2:          (geneC)             (geneC)
    #                 genome2: MG1:          (geneX)             (geneX)
    #                          MG2:          (geneY)       (geneZ,geneY)
    #                     dict dict          list

    all_predicted = dict()
    for mg in database_files:
        fna_path, faa_path = extract_genes_from_fasta(mg, selected_genes, genomes_pred, verbose, mg_info_use_protein[mg])
        all_predicted[mg] = [fna_path, faa_path]

    return all_predicted




# ==============================================================================
# TAXONOMICALLY ANNOTATE MARKER GENES
# ==============================================================================
# position of the script -------------------------------------------------------
path_this = os.path.realpath(__file__)
path_array = path_this.split("/")
stag_path = "/".join(path_array[0:-2]) + "/stag"

# we run stag classify, for each marker gene
def annotate_MGs(MGS, database_files, database_base_path):
    all_classifications = dict()
    for mg in MGS:
        if MGS[mg][0] != None:
            # it means that there are some genes to classify
            CMD = stag_path + " classify -d "+database_base_path+"/"+mg
            # check that the database is correct
            if not os.path.isfile(database_base_path+"/"+mg):
                sys.stderr.write("Error: file for gene database is missing")
            CMD = CMD + " -i "+MGS[mg][0]
            if MGS[mg][1] != "no_protein":
                # it means that we align proteins
                CMD = CMD + " -p "+MGS[mg][1]
            # we run stag CMD
            split_CMD = shlex.split(CMD)
            stag_CMD = subprocess.Popen(split_CMD, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            # save stderr for the message is necessary
            all_stderr = ""
            for line in stag_CMD.stderr:
                line = line.decode('ascii')
                all_stderr = all_stderr + line
            # save stdout with the resutls
            n_genes = -1
            for line in stag_CMD.stdout:
                n_genes = n_genes + 1
                if n_genes != 0:
                    # we skip the header
                    vals = line.decode('ascii').split("\t")
                    all_classifications[vals[0]] = vals[1].rstrip()
            # check errors
            return_code = stag_CMD.wait()
            if return_code:
                sys.stderr.write("[E::align] Error. stag classify failed\n\n")
                sys.stderr.write(all_stderr)
                sys.exit(1)

    return all_classifications



# ==============================================================================
# MERGE TAXONOMY OF SINGLE GENES
# ==============================================================================
def merge_genes_predictions(genomes_file_list, mgs_list, all_classifications, verbose, threads, output, long_out, keep_all_genes):
    # we parse "all_classifications"
    parsed_data = dict()
    for i in all_classifications:
        vals = i.rstrip().split("##")
        genome = "_".join(vals[0].split("_")[0:-1])
        if not genome in parsed_data:
            parsed_data[genome] = dict()
        if not vals[1] in parsed_data[genome]:
            parsed_data[genome][vals[1]] = list()
        parsed_data[genome][vals[1]].append(all_classifications[i])

    # we go throught the genomes and analyse them
    for g in genomes_file_list:
        if not g in parsed_data:
            # it means that genes were not idenitified or not classified
            print(g+"\tNone")
            break

        # annotation for each tax level
        levels_annotation = dict()
        for mg in parsed_data[g]:
            for tax in parsed_data[g][mg]:
                l = 0
                if tax != "":
                    for taxa in tax.split(";"):
                        l = l + 1
                        if not l in levels_annotation:
                            levels_annotation[l] = dict()
                        if not taxa in levels_annotation[l]:
                            levels_annotation[l][taxa] = 1
                        else:
                            levels_annotation[l][taxa] = levels_annotation[l][taxa] + 1

        # levels_annotation: 1: Bacteria: 3
        #                    2: Firmicutes: 2
        #                       Bacteroidetes: 1
        # Find best annotaiton per level
        genome_annotation = list()
        for l in levels_annotation:
            max = 0
            sel_taxa = ""
            for taxa in levels_annotation[l]:
                if levels_annotation[l][taxa] > max:
                    max = levels_annotation[l][taxa]
                    sel_taxa = taxa
            genome_annotation.append(sel_taxa)
        print(g+"\t"+";".join(genome_annotation))



#===============================================================================
#                                      MAIN
#===============================================================================
def classify_genome(database, genomes_file_list, verbose, threads, output, long_out, tool_version, keep_all_genes):
    # ZERO: we need to check that the genome files do not contain "##"
    for g in genomes_file_list:
        if len(g.split("##")) > 1:
            sys.stderr.write("Error with: "+g+"\n")
            sys.stderr.write("[E::main] Error: file cannot have in the name '##'. Please, choose anothe name.\n")
            sys.exit(1)

    # FIRST: unzip the database ------------------------------------------------
    if verbose > 2:
        sys.stderr.write("Unzip the database\n")
    database_files, temp_dir, gene_thresholds = load_genome_DB(database, tool_version, verbose)

    # SECOND: run prodigal on the fasta genome ---------------------------------
    if verbose > 2:
        sys.stderr.write("Run prodigal\n")
    genomes_pred = run_prodigal_genomes(genomes_file_list, verbose)
    # genomes_pred is a dictionary where the keys are the genome paths and the
    # values are lists. First value of the list is the path to the gene file and
    # second the path to the protein file

    # THIRD: find the marker genes from the predicted genes --------------------
    if verbose > 2:
        sys.stderr.write("Extract the marker genes\n")
    MGS = fetch_MGs(database_files, temp_dir, genomes_pred, keep_all_genes, gene_thresholds, verbose)
    # MGS = {'COG12':['path/to/genes','path/to/proteins'],
    #        'COG18':['path/to/genes','path/to/proteins'],
    #        ...}
    # not that the header of the sequences is:
    # genome_path + "_" + a_number + "##" + marker_gene_path
    # Example:
    # "User/Desktop/test_genome.fna_342##COG0012"

    # note: MGS = {'COG12':[None,None]
    # means that there are no genes detected to classify

    # note: MGS = {'COG12':["path/to/genes","no_protein"]
    # means that the alignment should be done at the level of the genes and not
    # proteins

    # we save in the outdir the file with the MG sequences
    os.mkdir(output+"/MG_sequences")
    for m in MGS:
        try:
            shutil.move(MGS[m][0],output+"/MG_sequences/"+m+".fna")
            shutil.move(MGS[m][1],output+"/MG_sequences/"+m+".faa")
            MGS[m][0] = output+"/MG_sequences/"+m+".fna"
            MGS[m][1] = output+"/MG_sequences/"+m+".faa"
        except:
            sys.stderr.write("[E::main] Error: failed to save the marker gene sequences\n")
            sys.exit(1)


    # FOURTH: classify the marker genes ----------------------------------------
    if verbose > 2:
        sys.stderr.write("Taxonomically annotate marker genes\n")
    all_classifications = annotate_MGs(MGS, database_files, temp_dir)
    # all_classifications is a dict: 'genome_id_NUMBER##cog_id': taxonomy
    #
    # Example:
    # '/Users/alex/Dropbox/genomeA_356##COG0012': "Bacteria;Firmicutes"
    # '/Users/alex/Dropbox/genomeA_51##COG0012': "Bacteria;"
    # '/Users/alex/Dropbox/genomeA_784##COG0018': "Bacteria;Firmicutes;Bacilli"
    # '/Users/alex/Dropbox/genomeBB_1853##COG0012': "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales"
    # '/Users/alex/Dropbox/genomeBB_862##COG0172': "Bacteria;Bacteroidetes;Bacteroidia"

    # we remove the temp dir ---------------------------------------------------
    shutil.rmtree(temp_dir)
    # and the result from prodigal
    for i in genomes_pred:
        if os.path.isfile(genomes_pred[i][0]): os.remove(genomes_pred[i][0])
        if os.path.isfile(genomes_pred[i][1]): os.remove(genomes_pred[i][1])
    # and the file with the marker genes
    #for m in MGS:
    #    if MGS[m][0] != None:
    #        if os.path.isfile(MGS[m][0]): os.remove(MGS[m][0])
    #        if os.path.isfile(MGS[m][1]): os.remove(MGS[m][1])


    # FIFTH: join prediction ---------------------------------------------------
    if verbose > 2:
        sys.stderr.write("Join taxonomy of different genes\n")
    merge_genes_predictions(genomes_file_list, list(database_files), all_classifications, verbose, threads, output, long_out, keep_all_genes)
