import os
import tempfile
import shutil
import tarfile

from stag.helpers import check_file_exists
from stag.classify import classify

# cschu 2021-04-10: we need to change the alignment format!! -> this is too hacky.
# find the length of the alignments --------------------------------------------
def find_length_ali(gene_db, fasta_input, protein_fasta_input):
    return classify(gene_db, fasta_input=fasta_input,
                    protein_fasta_input=protein_fasta_input, internal_call=True)[0]

def get_dummy_fastas():
    fasta_files = list()
    for seq in ("AAA", "A"):
        with tempfile.NamedTemporaryFile(delete=False, mode="w") as tmp_fasta:
            os.chmod(tmp_fasta.name, 0o644)
            print(">test", seq, sep="\n", file=tmp_fasta, flush=True)
            fasta_files.append(tmp_fasta.name)
    return fasta_files

def get_alignment_lengths(list_genes):
    fna, faa = get_dummy_fastas()
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as length_file:
        os.chmod(length_file.name, 0o644)
        for gene_db in list_genes:
            print(os.path.basename(gene_db), find_length_ali(gene_db, fna, faa), sep="\t", flush=True, file=length_file)
        [os.remove(f) for f in (fna, faa)]
        return length_file.name


def train_genome(output, list_genes, gene_threshold_file, threads, verbose, concat_stag_db):
    check_file_exists(gene_threshold_file, isfasta=False)
    with open(gene_threshold_file) as f:
        gene_thresholds = set(line.strip().split("\t")[0] for line in f if line)

    list_genes = list_genes.split(",")
    missing_thresholds = set(os.path.basename(fn) for fn in list_genes).difference(gene_thresholds)
    if missing_thresholds:
        raise ValueError(f"[E::main] Error: gene {list(missing_thresholds)[0]} is missing from the threshold file (-T)")

    outfile = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(outfile.name, 0o644)
    core_db_files = ("threshold_file.tsv", "hmm_lengths_file.tsv", "concatenated_genes_STAG_database.HDF5")
    with tarfile.open(outfile.name, "w:gz", dereference=True) as genome_tar:
        for fn in list_genes:
            check_file_exists(fn)
            base_fn = os.path.basename(fn)
            if base_fn in core_db_files:
                raise ValueError(f"[E::main] Error: gene databases cannot be named '{base_fn}'. Please choose another name.")
            if "##" in base_fn:
                raise ValueError(f"Error with: {base_fn}\n[E::main] Error: gene database file names cannot contain '##'. Please choose another name.")
            try:
                genome_tar.add(fn, base_fn)
            except:
                raise ValueError(f"[E::main] Error: when adding {fn} to the database")
        for source, target in zip((gene_threshold_file, get_alignment_lengths(list_genes), concat_stag_db), core_db_files):
            genome_tar.add(source, target)

    try:
        outfile.flush()
        os.fsync(outfile.fileno())
        outfile.close()
    except:
        raise ValueError("[E::main] Error: failed to save the result.")
    try:
        shutil.move(outfile.name, output)
    except:
        raise ValueError("[E::main] Error: failed to save the resulting database\n" + \
                         f"[E::main] you can find the file here:\n{outfile.name}")
