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
        #filter lines
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
def extract_genes_from_one_genome(file_to_align, hmm_file, gene_threshold):
    # INFO: genes_path, proteins_path [where to save the result]
    # we run hmmsearch
    temp_hmm = tempfile.NamedTemporaryFile(delete=False, mode="w")
    hmm_cmd = "hmmsearch --tblout "+temp_hmm.name+" "+hmm_file+" "+file_to_align

    CMD = shlex.split(hmm_cmd)
    hmm_CMD = subprocess.Popen(CMD, stdout=DEVNULL,stderr=subprocess.PIPE)
    # we save stderr if necessary
    all_stderr = ""
    for line in hmm_CMD.stderr:
        #filter lines
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
def extract_genes(mg_name, hmm_file, use_protein_file, genomes_pred, gene_threshold):
    # two temp files that will contain all the MGs (of one type) for all genomes
    genes = tempfile.NamedTemporaryFile(delete=False, mode="w")
    if use_protein_file:
        proteins = tempfile.NamedTemporaryFile(delete=False, mode="w")
        proteins_n = proteins.name
    else:
        proteins_n = ""
        proteins = ""
    # we go throught the genome and find the genes that pass the filter
    genes_pass_filter = dict()
    for g in genomes_pred:
        if use_protein_file:
            file_to_align = genomes_pred[g][1]
        else:
            file_to_align = genomes_pred[g][0]
        genes_pass_filter[g] = extract_genes_from_one_genome(file_to_align, hmm_file, gene_threshold)
    return genes_pass_filter

# extract the marker genes from the genes/proteins produced from prodigal
# for multiple genomes and multiple MGs
def fetch_MGs(database_files, database_path, genomes_pred, keep_all_genes, gene_thresholds):
    all_genes_raw = dict()
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
        f.close()

        # run hmmsearch for each genome and find which genes pass the filter
        all_genes_raw[mg] = extract_genes(mg, hmm_file.name, use_protein_file, genomes_pred, gene_thresholds[mg])
        # the result is a dict: genome -> dict with genes -> score

        # remove hmm file
        os.remove(hmm_file.name)

    # now we need to select the marker genes and extracted them from the prodigal
    # fasta files
    # all_genes_raw: MG1: genome1: geneA: 276
    #                              geneB: 243
    #                     genome2: geneX: 267
    #                              geneY: 212
    #                MG2: genome1: geneC: 589
    #                     genome2: geneZ: 459
    #                              geneY: 543
    # NOTE: the same gene can apper in two different marker genes. Hence, we
    # need to assign it to only one
    # For example 'geneY' is in both MG1 and MG2 for genome 2.
    # we will select MG2 because the score is higher
    print(all_genes_raw)


    all_predicted = dict()
    fna_path, faa_path = extract_genes(mg, hmm_file.name, use_protein_file, genomes_pred, keep_all_genes, gene_thresholds[mg])
    all_predicted[mg] = [fna_path, faa_path]
    return all_predicted

#===============================================================================
#                                      MAIN
#===============================================================================
def classify_genome(database, genomes_file_list, verbose, threads, output, long_out, tool_version, keep_all_genes):
    # FIRST: unzip the database
    database_files, temp_dir, gene_thresholds = load_genome_DB(database, tool_version, verbose)

    # SECOND: run prodigal on the fasta genome
    genomes_pred = run_prodigal_genomes(genomes_file_list, verbose)
    # genomes_pred is a dictionary where the keys are the genome paths and the
    # values are lists. First value of the list is the path to the gene file and
    # second the path to the protein file

    # THIRD: find the marker genes from the predicted genes
    MGS = fetch_MGs(database_files, temp_dir, genomes_pred, keep_all_genes, gene_thresholds)
    # MGS = {'COG12':['path/to/genes','path/to/proteins'],
    #        'COG18':['path/to/genes','path/to/proteins'],}
    # if 'path/to/proteins' == "", then the alignment is with genes

    # FOURTH: classify the marker genes

    # we remove the temp dir
    shutil.rmtree(temp_dir)
    # and the result from prodigal
    for i in genomes_pred:
        if os.path.isfile(genomes_pred[i][0]): os.remove(genomes_pred[i][0])
        if os.path.isfile(genomes_pred[i][1]): os.remove(genomes_pred[i][1])


    # FIFTH: join prediction
