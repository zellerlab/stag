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
import pathlib

import contextlib

from stag.helpers import is_tool, read_fasta
from stag.classify import classify

try:
    from subprocess import DEVNULL
except ImportError:
    DEVNULL = open(os.devnull, 'wb')

class MarkerGeneExtractor:
    pass

class ProdigalWrapper:
    def __init__(self):
        if not is_tool("prodigal"):                                                     
            raise ValueError("[E::align] Error: prodigal is not in the path.\n")
    @staticmethod
    def run_prodigal(genome_fasta, gene_file, protein_file):
        prodigal_cmd = f"prodigal -i {genome_fasta} -d {gene_file} -a {protein_file}"
        parse_cmd = subprocess.Popen(prodigal_cmd, shell=True, stdout=DEVNULL, stderr=subprocess.PIPE)
        all_stderr = [line.decode("ascii") for line in parse_cmd.stderr] 
        return_code = parse_cmd.wait()
        if return_code:
            msg = ["[E::align] Error. prodigal failed", ""] + all_stderr
            raise ValueError("\n".join(msg))
    @staticmethod
    def cleanup_prodigal(files):
        for genes, proteins in files:
            try:
                [os.remove(f) for f in (genes, proteins)]
            except:
                pass


class MarkerGeneAnnotator:
    @staticmethod
    def copy_fasta(fasta_in, fasta_out):
        for index, (sid, seq) in enumerate(read_fasta(fasta_in, is_binary=False)):
            print(">{genome}_{index}".format(**locals()), seq, sep="\n", file=fasta_out)
        fasta_out.flush()
        os.fsync(fasta_out.fileno())
        return index + 1                                                                    
    @staticmethod
    def annotate_genomes(genomes):
        return {genome: MarkerGeneAnnotator(genome).run()
                for genome in genomes}

    def __init__(self, genome_fasta):
        self.genome = genome_fasta        
    def run(self):
        genes = tempfile.NamedTemporaryFile(delete=False, mode="w")
        proteins = tempfile.NamedTemporaryFile(delete=False, mode="w")
        ProdigalWrapper().run_prodigal(self.genome, genes.name, proteins.name)

        parsed_genes = tempfile.NamedTemporaryFile(delete=False, mode="w")
        parsed_proteins = tempfile.NamedTemporaryFile(delete=False, mode="w")

        with parsed_genes, open(genes.name) as genes_in:
            n_genes = MarkerGeneAnnotator.copy_fasta(genes_in, parsed_genes)
        with parsed_proteins, open(proteins.name) as proteins_in:
            n_proteins = MarkerGeneAnnotator.copy_fasta(proteins_in, parsed_proteins)

        os.remove(genes.name)
        os.remove(proteins.name)

        return parsed_genes.name, parsed_proteins.name

def save_sequences(seq_source, prefix, suffix, outdir):
    target = None
    try:
        if not seq_source:
            open(target, "w").close()
        else:
            target = os.path.join(outdir, "{prefix}.{suffix}".format(**locals()))
            shutil.move(seq_source, target)
    except Exception as e:
        msg = ["[E::main] Error: failed to save the marker gene sequences", str(e)]
        raise ValueError("\n".join(msg))
        

def save_marker_genes(marker_genes, outdir):
    outdir = os.path.join(outdir, "MG_sequences")
    pathlib.Path(outdir).mkdir(exist_ok=True, parents=True)

    for mg, (genes, proteins) in marker_genes.items():
        marker_genes[mg] = [
            save_sequences(genes, mg, "fna", outdir),
            save_sequences(proteins, mg, "faa", outdir)
        ]
        

# ==============================================================================
# EXTRACT THE MARKER GENES
# ==============================================================================
# find gene ids that we can use (run hmmsearch)
def extract_gene_from_one_genome(file_to_align, hmm_file, gene_threshold, mg_name):
    # INFO: genes_path, proteins_path [where to save the result]
    # we run hmmsearch
    temp_hmm = tempfile.NamedTemporaryFile(delete=False, mode="w")
    cmd = "hmmsearch --tblout {} {} {}".format(temp_hmm.name, hmm_file, file_to_align)

    hmm_CMD = subprocess.Popen(cmd, stdout=DEVNULL, stderr=subprocess.PIPE, shell=True)
    # we save stderr if necessary
    all_stderr = list(line.decode("ascii") for line in hmm_CMD.stderr)
    return_code = hmm_CMD.wait()
    if return_code:
        msg = [
            "[E::align] Error. hmmsearch failed", "",
            "MG: {}".format(mg_name),
            "CALL: {}".format(hmm_CMD), ""
        ] + all_stderr
        raise ValueError("\n".join(msg))

    # in temp_hmm.name there is the result from hmm ----------------------------
    # we select which genes/proteins we need to extract from the fasta files
    # produced by prodigal
    gene_threshold = float(gene_threshold) # should be float before!
    sel_genes = dict()
    with open(temp_hmm.name, "r") as hmm_in:
        for line in hmm_in:
            if not line.startswith("#"):
                vals = re.sub(" +", " ", line.rstrip()).split(" ")
                gene_id, e_val, score, = vals[0], vals[4], float(vals[5])
                if score > gene_threshold:
                    sel_genes[gene_id] = score

    os.remove(temp_hmm.name)
    return sel_genes


def select_genes(all_genes_raw, keep_all_genes):
    selected_genes = dict()
    # all_genes_raw: genome1: MG1: geneA: 276
    #                              geneB: 243
    #                         MG2: geneC: 589
    #                genome2: MG1: geneX: 267
    #                              geneY: 212
    #                         MG2: geneZ: 459
    #                              geneY: 543
    for genome, (marker_candidates, _) in all_genes_raw.items():
        # we first check if there is any gene that is in multiple mgs:
        best_score = dict()
        for mg, genes in marker_candidates.items():
            for gene, score in genes.items():
                score = genes[gene] = float(score) # this should be float before!
                best_score[gene] = max(score, best_score.get(gene, 0.0))
        # in gene_sel there is the gene id -> highest score
        # example: geneX->267; geneY->543; geneZ->459

        # now we select the correct genes and decide if keep one or many
        for mg, genes in marker_candidates.items():
            selected_genes.setdefault(genome, dict())[mg] = list()
            
            best_genes = [
                (score, gene) for gene, score in marker_candidates.items()
                if score == best_score[gene]
            ]
            
            if keep_all_genes:
                selected_genes[genome][mg].extend(gene for gene, _ in best_genes)
            else:
                selected_genes[genome][mg].append(sorted(best_genes, key=lambda x:x[0], reverse=True)[0][1])
                    
    return selected_genes

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
        for genome, marker_genes in selected_genes.items():
            mg_genes = set(marker_genes.get(mg, list()))
            if not mg_genes:
                sys.stderr.write("Warning: missing marker gene in genome "+genome+"\n")        
            else:
                with open(genomes_pred[genome][0]) as fna_in:
                    n_genes += filter_sequences(fna_in, genes, mg_genes, mg)
                if use_protein_file:
                    with open(genomes_pred[genome][1]) as faa_in:
                        n_proteins += filter_sequences(faa_in, proteins, mg_genes, mg)

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


def write_hmm(db_in):
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as hmm_file:
        os.chmod(hmm_file.name, 0o644)
        hmm_file.write(db_in['hmm_file'][0])
        hmm_file.flush()
        os.fsync(hmm_file.fileno())

    return hmm_file.name

# extract the marker genes from the genes/proteins produced from prodigal
# for multiple genomes and multiple MGs
def extract_marker_genes(database_files, database_path, genomes_pred, gene_thresholds):
    all_genes_raw = dict()
    for mg in database_files:
        # for each MG, we extract the hmm and if using proteins or not ---------
        path_mg = os.path.join(database_path, mg)
        with h5py.File(path_mg, 'r') as db_in:
            hmm_file = write_hmm(db_in)
            use_protein_file = bool(db_in['align_protein'][0])

        # run hmmsearch for each genome and find which genes pass the filter
        # the result is saved in all_genes_raw (passed as input)
        # extract_genes(mg, hmm_file, mg_info_use_protein[mg], genomes_pred, gene_thresholds[mg], all_genes_raw)
        for genome, genes_proteins in genomes_pred.items():
            mg_present = all_genes_raw.setdefault(genome, dict()).get(mg)
            if mg_present:
                sys.stderr.write("Error. gene already present\n")
            extracted = extract_gene_from_one_genome(
                genes_proteins[int(use_protein_file)],
                hmm_file, gene_thresholds[mg], mg
            )
            all_genes_raw[genome][mg] = (extracted, use_protein_file)

        os.remove(hmm_file)

    return all_genes_raw

def filter_marker_genes(all_genes_raw, keep_all_genes, genomes_pred, verbose):

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
    for mg, (extracted, use_proteins) in all_genes_raw.items():
        fna_path, faa_path = extract_genes_from_fasta(mg, selected_genes, genomes_pred, verbose, use_proteins)
        all_predicted[mg] = [fna_path, faa_path]

    return all_predicted




# ==============================================================================
# TAXONOMICALLY ANNOTATE MARKER GENES
# ==============================================================================

def annotate_MGs(MGS, database_files, database_base_path, dir_ali):
    all_classifications = dict()
    for mg, (fna, faa) in MGS.items():
        if fna:
            db = os.path.join(database_base_path, mg)
            if not os.path.isfile(db):
                sys.stderr.write("Error: file for gene database {} is missing".format(db))
                sys.exit(1)
            # faa = faa if faa != "no_protein" else None
            align_out = os.path.join(dir_ali, mg)
            _, results = classify(db, fasta_input=fna, protein_fasta_input=faa,
                                  save_ali_to_file=align_out, internal_call=True)
            all_classifications.update(dict(results))

    return all_classifications 
                 

# ==============================================================================
# MERGE TAXONOMY OF SINGLE GENES
# ==============================================================================
def merge_genes_predictions(genomes_file_list, mgs_list, all_classifications, verbose, threads, output, long_out, keep_all_genes):
    # we parse "all_classifications"
    to_print = dict()
    for i in all_classifications:
        vals = i.rstrip().split("##")
        genome = "_".join(vals[0].split("_")[0:-1])
        mg_id = vals[1]
        if not genome in to_print:
            to_print[genome] = ""
        to_print[genome] = to_print[genome] + i.rstrip() + "\t" + mg_id + "\t" + all_classifications[i] + "\n"

    # we go throught the genomes and analyse them
    for g in genomes_file_list:
        if not g in to_print:
            to_print[g] = ""
        genome_file_name = g.split("/")[-1]
        o = open(output+"/genes_predictions/"+genome_file_name,"w")
        o.write(to_print[g])
        o.close()



# ==============================================================================
# CONCAT ALIGNEMENTS
# ==============================================================================
def concat_alis(genomes_file_list, ali_dir, gene_order, ali_lengths):
    # we return a (tmp) file containing the concatenated alignment
    # INPUT:
    #  - list of genomes
    #  - base name of the directory containing the alignments
    #  - order of the genes
    #  - length of the alignments

    # we create the base
    all_genes = dict()
    for ge in genomes_file_list:
        all_genes[ge] = list()
        for mg in gene_order:
            all_genes[ge].append("\t".join(['0'] * int(ali_lengths[mg])))

    # we add the alignments from the real genes
    pos = -1
    for mg in gene_order:
        pos = pos + 1
        if os.path.isfile(ali_dir+mg):
            o = open(ali_dir+mg,"r")
            for line in o:
                genome = "_".join(line.split("##")[0].split("_")[0:-1])
                all_genes[genome][pos] = "\t".join(line.split("\t")[1:]).rstrip()
            o.close()

    # we create a temp file and save the sequences
    concat_ali_f = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(concat_ali_f.name, 0o644)
    for g in genomes_file_list:
        str_g = g.split("/")[-1] + "\t"
        str_g = str_g + "\t".join(all_genes[g])
        concat_ali_f.write(str_g + "\n")
        concat_ali_f.flush()

    return concat_ali_f.name


# ==============================================================================
# ANNOTATE concatenation of the MGs
# ==============================================================================
# the annotation for the genome is based on the annotation of the
def annotate_concat_mgs(stag_db, alignment_file, output_dir):
    _, results = classify(stag_db, aligned_sequences=alignment_file, output=os.path.join(output_dir, "genome_annotation"))


#===============================================================================
#                                      MAIN
#===============================================================================
def classify_genome(database, genomes_file_list, verbose, threads, output_dir, long_out, tool_version, keep_all_genes):
    # FIRST: unzip the database ------------------------------------------------
    if verbose > 2:
        sys.stderr.write("Unzip the database\n")
    database_files, temp_dir, gene_thresholds, gene_order, ali_lengths, concat_ali_stag_db = load_genome_DB(database, tool_version, verbose)

    # SECOND: run prodigal on the fasta genome ---------------------------------
    if verbose > 2:
        sys.stderr.write("Run prodigal\n")
    genomes_pred = MarkerGeneAnnotator.annotate_genomes(genomes_file_list)
    # genomes_pred is a dictionary where the keys are the genome paths and the
    # values are lists. First value of the list is the path to the gene file and
    # second the path to the protein file

    # THIRD: find the marker genes from the predicted genes --------------------
    if verbose > 2:
        sys.stderr.write("Extract the marker genes\n")

    raw_marker_genes = extract_marker_genes(database_files, temp_dir, genomes_pred, gene_thresholds)
    marker_genes = filter_marker_genes(raw_marker_genes, keep_all_genes, genomes_pred, verbose)
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

    if not any(genes for genes, _ in marker_genes.values()):
        shutil.rmtree(temp_dir)
        cleanup_prodigal(genomes_pred.values())
        raise ValueError("[W::main] Warning: no marker genes identified. Stopping annotation.")
    
    save_marker_genes(marker_genes, outdir)

    # FOURTH: classify the marker genes ----------------------------------------
    if verbose > 2:
        sys.stderr.write("Taxonomically annotate single marker genes\n")

    # when doing the classification, we also create the alignment files
    align_dir = os.path.join(output_dir, "MG_ali")
    pathlib.Path(align_dir).mkdir(exist_ok=True, parents=True)
    all_classifications = annotate_MGs(marker_genes, database_files, temp_dir, align_dir)
    # all_classifications is a dict: 'genome_id_NUMBER##cog_id': taxonomy
    #
    # Example:
    # '/Users/alex/Dropbox/genomeA_356##COG0012': "Bacteria;Firmicutes"
    # '/Users/alex/Dropbox/genomeA_51##COG0012': "Bacteria;"
    # '/Users/alex/Dropbox/genomeA_784##COG0018': "Bacteria;Firmicutes;Bacilli"
    # '/Users/alex/Dropbox/genomeBB_1853##COG0012': "Bacteria;Bacteroidetes;Bacteroidia;Bacteroidales"
    # '/Users/alex/Dropbox/genomeBB_862##COG0172': "Bacteria;Bacteroidetes;Bacteroidia"


    # join prediction ----------------------------------------------------------
    prediction_dir = os.path.join(output_dir, "genes_predictions")
    pathlib.Path(prediction_dir).mkdir(exist_ok=True, parents=True)
    merge_genes_predictions(genomes_file_list, list(database_files), all_classifications, verbose, threads, output_dir, long_out, keep_all_genes)

    # FIFTH: classify the concatenation of the MGs, which represents the -------
    #         annotation for the genome ----------------------------------------
    if verbose > 2:
        sys.stderr.write("Taxonomically annotate genomes\n")
    # First, create a concatenated alignment. The alignments were created in the
    # 4th step
    file_ali = concat_alis(genomes_file_list, align_dir, gene_order, ali_lengths)

    # Second, classify the alignments
    annotate_concat_mgs(concat_ali_stag_db, file_ali, output_dir)

    # we remove the file with the concatenated alignment
    os.remove(file_ali)

    # we remove the temp dir ---------------------------------------------------
    shutil.rmtree(temp_dir)
    # and the result from prodigal
    if genomes_file_list:
        ProdigalWrapper.cleanup_prodigal(genomes_pred.values())
