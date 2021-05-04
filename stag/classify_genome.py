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
import json
import pathlib

import contextlib
import multiprocessing as mp

from stag.helpers import is_tool, read_fasta
from stag.classify import classify
from stag import __version__ as tool_version
from stag.databases import load_genome_DB

try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

def validate_genome_files(files):
    if any("##" in f for f in files):
        offender = [f for f in files if "##" in f][0]
        err = "Error with: {}\n".format(offender)
        raise ValueError(f"{err}[E::main] Error: file cannot have in the name '##'. Please, choose another name.\n")

def cleanup_prodigal(files):
    for genes, proteins in files:
        try:
            [os.remove(f) for f in (genes, proteins)]
        except:
            pass


def cleanup_prodigal(files):
    for genes, proteins in files:
        try:
            [os.remove(f) for f in (genes, proteins)]
        except:
            pass

def run_prodigal(genome):
    if not is_tool("prodigal"):
        raise ValueError("[E::align] Error: prodigal is not in the path.\n")

    # we need two files, one for the proteins and one for the genes
    genes = tempfile.NamedTemporaryFile(delete=False, mode="w")
    proteins = tempfile.NamedTemporaryFile(delete=False, mode="w")
    # prodigal command
    prodigal_cmd = "prodigal -i {genome} -d {gene_file} -a {protein_file}".format(
        genome=genome, gene_file=genes.name, protein_file=proteins.name
    )
    cmd = shlex.split(prodigal_cmd)
    parse_cmd = subprocess.Popen(cmd, stdout=DEVNULL,stderr=subprocess.PIPE)
    # we save stderr if necessary
    all_stderr = ""
    for line in parse_cmd.stderr:
        line = line.decode('ascii')
        all_stderr = all_stderr + line
    return_code = parse_cmd.wait()
    if return_code:
        raise ValueError(f"[E::align] Error. prodigal failed\n\n{all_stderr}")

    # we re-name the header of the fasta files ---------------------------------
    # we expect to have the same number of genes and proteins, and also that the
    def copy_fasta(fasta_file, seqid, is_binary=True, head_start=0):
        with tempfile.NamedTemporaryFile(delete=False, mode="w") as fasta_out, open(fasta_file) as fasta_in:
            for index, (sid, seq) in enumerate(read_fasta(fasta_in, is_binary=is_binary), start=1):
                print(">{seqid}_{index}".format(**locals()), seq, sep="\n", file=fasta_out)
            fasta_out.flush()
            os.fsync(fasta_out.fileno())
            return fasta_out.name, index

    parsed_genes, gene_count = copy_fasta(genes.name, genome, is_binary=False)
    parsed_proteins, protein_count = copy_fasta(proteins.name, genome, is_binary=False)

    os.remove(genes.name)
    os.remove(proteins.name)

    return parsed_genes, parsed_proteins

# run prodigal on all the genomes listed in fasta_input
def run_prodigal_genomes(genome_files):
    return {genome: run_prodigal(genome) for genome in genome_files}

# ==============================================================================
# EXTRACT THE MARKER GENES
# ==============================================================================
# find gene ids that we can use (run hmmsearch)
def extract_gene_from_one_genome(file_to_align, hmm_file, gene_threshold,mg_name):
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
        raise ValueError(f"[E::align] Error. hmmsearch failed\n\nMG: {mg_name}\nCALL: {hmm_cmd}\n\n{all_stderr}")

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
        all_genes_raw[g][mg_name] = extract_gene_from_one_genome(file_to_align, hmm_file, gene_threshold,mg_name)

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
    n_genes, n_proteins = 0, 0
    genes = tempfile.NamedTemporaryFile(delete=False, mode="w")
    if use_protein_file:
        proteins = tempfile.NamedTemporaryFile(delete=False, mode="w")
    else:
        proteins = contextlib.nullcontext()

    def filter_sequences(fasta_in, fasta_out, whitelist, mg):
        n_written = 0
        for sid, seq in read_fasta(fasta_in, head_start=1, is_binary=False):
            if sid in whitelist:
                n_written += 1
                print(">{sid}##{mg}".format(sid=sid, mg=mg), seq, sep="\n", file=fasta_out)
        return n_written

    with genes, proteins:
        for genome in selected_genes:
            if not mg in selected_genes[genome]:
                sys.stderr.write("Warning: missing marker gene in genome "+genome+"\n")
            else:
                with open(genomes_pred[genome][0]) as fna_in:
                    n_genes += filter_sequences(fna_in, genes, selected_genes[genome][mg], mg)
                if use_protein_file:
                    with open(genomes_pred[genome][1]) as faa_in:
                        n_proteins += filter_sequences(faa_in, proteins, selected_genes[genome][mg], mg)

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
            protein_file_name = None #"no_protein"

    return gene_file_name, protein_file_name


# extract the marker genes from the genes/proteins produced from prodigal
# for multiple genomes and multiple MGs
def fetch_MGs(database_files, database_path, genomes_pred, keep_all_genes, gene_thresholds, verbose):
    all_genes_raw = dict()
    mg_info_use_protein = dict()
    for mg in database_files:
        # for each MG, we extract the hmm and if using proteins or not ---------
        path_mg = os.path.join(database_path, mg)

        with h5py.File(path_mg, 'r') as db_in, tempfile.NamedTemporaryFile(delete=False, mode="w") as hmm_file:
            os.chmod(hmm_file.name, 0o644)
            hmm_file.write(db_in['hmm_file'][0])
            hmm_file.flush()
            os.fsync(hmm_file.fileno())

            use_protein_file = mg_info_use_protein[mg] = bool(db_in['align_protein'][0])

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

def annotate_MGs(MGS, database_files, database_base_path, dir_ali, procs=2):

    found_marker_genes = {
        mg: (fna, faa) 
        for mg, (fna, faa) in MGS.items()
        if os.path.exists(fna) and pathlib(fna).stat().st_size
    }

    if not found_marker_genes:
        raise ValueError("No marker genes found!")

    for mg in found_marker_genes:
        db = os.path.join(database_base_path, mg)
        if not os.path.isfile(db):
            raise ValueError(f"Error: file for gene database {db} is missing")

    pool = mp.Pool(processes=procs)

    results = (
        pool.apply_async(
            classify,
            args=(os.path.join(database_base_path,mg),),
            kwds={"fasta_input": fna, "protein_fasta_input": faa,
                  "save_ali_to_file": os.path.join(dir_ali, mg),
                  "internal_call": True}
        )
        for mg, (fna, faa) in found_marker_genes.items()
    )

    d = dict()
    for p in results:
        _, predictions = p.get()
        d.update(predictions)
    return d
    #return dict(p.get()[1] for p in results)

def merge_gene_predictions(genome_files, mgs_list, all_classifications, verbose, threads, output, long_out, keep_all_genes, full_genomes=True):
    outdir = os.path.join(output, "genes_predictions")
    pathlib.Path(outdir).mkdir(exist_ok=True, parents=True)

    print(*all_classifications.items(), sep="\n")
    print("**********")
    # we parse "all_classifications"
    merged_predictions = dict()
    for marker_gene, lineage in all_classifications.items():
        genome, mg_id = marker_gene.rstrip().split("##")
        sep = "_" if "_" in genome else "."
        genome = genome.split(sep)
        genome = sep.join(genome[:-1] if len(genome) > 1 else genome)
        merged_predictions.setdefault(genome, list()).append("\t".join([marker_gene.rstrip(), mg_id, lineage]))
    print(*merged_predictions.items(), sep="\n")

    for genome, predictions in merged_predictions.items():
        genome_filename = os.path.basename(genome)
        with open(os.path.join(outdir, os.path.basename(genome)), "w") as merged_out:
            print(*predictions, sep="\n", file=merged_out, flush=True)

def concat_alignments(genome_files, ali_dir, gene_order, ali_lengths, full_genomes=True):
    all_genes = dict()
    for pos, mg in enumerate(gene_order):
        mg_alignment_file = os.path.join(ali_dir, mg)
        if os.path.exists(mg_alignment_file):
            with open(mg_alignment_file) as align_in:
                for line in align_in:
                    genome, *alignment = line.strip().split("\t")
                    sep = "_" if "_" in genome else "."
                    genome = genome.split("##")[0].split(sep)
                    genome = sep.join(genome[:-1] if len(genome) > 1 else genome)
                    all_genes.setdefault(genome, ["\t".join(["0" for i in range(int(ali_lengths[mg]))]) for mg in gene_order])
                    all_genes[genome][pos] = "\t".join(alignment)

    print("YYY")
    for key, value in all_genes.items():
        print(key, len(value))

    # we create a temp file and save the sequences
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as concat_ali_f:
        os.chmod(concat_ali_f.name, 0o644)
        for genome, alignment in all_genes.items():
            print(genome, *alignment, sep="\t", file=concat_ali_f, flush=True)

    return concat_ali_f.name


def store_marker_sequences(marker_sequences, outdir, genome_files_available=True):

    def store_seqs(source, target, copy_function):
        try:
            if not source:
                open(target, "w").close()
            else:
                copy_function(os.path.abspath(source), target)
        except Exception as e:
            raise ValueError(f"[E::main] Error: failed to save the marker gene sequences\n{e}")
        return target

    print(marker_sequences)
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=True)
    copy_function = shutil.move if genome_files_available else os.symlink

    for marker, (fna, faa) in marker_sequences.items():
        fna = store_seqs(fna, os.path.join(outdir, f"{marker}.fna"), copy_function)
        faa = store_seqs(faa, os.path.join(outdir, f"{marker}.faa"), copy_function)
        marker_sequences[marker] = [fna, faa]


def classify_genome(database, genome_files=None, marker_genes=None, verbose=None,
                    threads=1, output=None, long_out=False, keep_all_genes=False):

    # FIRST: unzip the database ------------------------------------------------
    if verbose > 2:
        sys.stderr.write("Unzip the database\n")
    database_files, temp_dir, gene_thresholds, gene_order, ali_lengths, concat_ali_stag_db = load_genome_DB(database, tool_version, verbose)
    genomes_pred = dict()
    print(*ali_lengths.items(), sep="\n")

    if marker_genes:
        MGS = json.load(open(marker_genes[0]))
    elif genome_files:
        # SECOND: run prodigal on the fasta genome ---------------------------------
        if verbose > 2:
            sys.stderr.write("Run prodigal\n")
        genomes_pred = run_prodigal_genomes(genome_files)
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

        # check if all genes are empty
        if not any(genes for genes, _ in MGS.values()):
            shutil.rmtree(temp_dir)
            cleanup_prodigal(genomes_pred.values())
            raise ValueError("[W::main] Warning: no marker genes identified\n          Stopping annotation.\n")

    store_marker_sequences(MGS, os.path.join(output, "MG_sequences"), genome_files_available=bool(genome_files))

    # FOURTH: classify the marker genes ----------------------------------------
    if verbose > 2:
        sys.stderr.write("Taxonomically annotate single marker genes\n")

    # when doing the classification, we also create the alignment files
    align_dir = os.path.join(output, "MG_ali")
    pathlib.Path(align_dir).mkdir(exist_ok=True, parents=True)

    all_classifications = annotate_MGs(MGS, database_files, temp_dir, align_dir, procs=threads)
    # all_classifications is a dict: 'genome_id_NUMBER##cog_id': taxonomy
    #
    # Example:
    # '/Users/alex/Dropbox/genomeA_356##COG0012': "Bacteria;Firmicutes"
    # '/Users/alex/Dropbox/genomeA_51##COG0012': "Bacteria;"
    # '/Users/alex/Dropbox/genomeA_784##COG0018': "Bacteria;Firmicutes;Bacilli"
    # '/Users/alex/Dropbox/genomeBB_1853##COG0012': "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales"
    # '/Users/alex/Dropbox/genomeBB_862##COG0172': "Bacteria;Bacteroidetes;Bacteroidia"

    # join prediction ----------------------------------------------------------
    input_files = genome_files if genome_files else marker_genes
    merge_gene_predictions(input_files, list(database_files), all_classifications,
                           verbose, threads, output, long_out, keep_all_genes,
                           full_genomes=bool(genome_files))

    # FIFTH: classify the concatenation of the MGs, which represents the -------
    #         annotation for the genome ----------------------------------------
    if verbose > 2:
        sys.stderr.write("Taxonomically annotate genomes\n")
    # First, create a concatenated alignment. The alignments were created in the
    # 4th step
    file_ali = concat_alignments(input_files, align_dir, gene_order, ali_lengths, full_genomes=bool(genome_files))

    # Second, classify the alignments
    classify(concat_ali_stag_db, aligned_sequences=file_ali, output=os.path.join(output, "genome_annotation"), long_out=long_out)

    # we remove the file with the concatenated alignment
    print(file_ali)
    #os.remove(file_ali)

    # we remove the temp dir ---------------------------------------------------
    shutil.rmtree(temp_dir)

    if genome_files:
        cleanup_prodigal(genomes_pred.values())
