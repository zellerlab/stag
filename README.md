<img src="https://github.com/AlessioMilanese/stag/blob/master/pics/stag_logo.png" width="450">

This tool is design to classify metagenomic sequences (marker genes, genomes and amplicon reads) using a Hierarchical Taxonomic Classifier.


Installation
--------------
```bash
git clone https://github.com/AlessioMilanese/stag.git
cd stag
export PATH=`pwd`:$PATH
```

Note: in the following examples we assume that the python script ```stag``` is in the system path.



Dependencies
--------------

The stag classifier requires:
* Python 3 (or higher)
* HMMER3 (or Infernal)
* Easel ([link](https://github.com/EddyRivasLab/easel))
* [seqtk](https://github.com/lh3/seqtk)
* python library:
  * numpy
  * pandas
  * sklearn
  * h5py

If you have [conda](https://conda.io/docs/), you can install the dependencies by creating an environment:
```bash
conda env create -f conda_env_stag.yaml
conda activate stag
```

You can check that the tool works properly with:
```
python test.py
```

Taxonomically annotate unknown sequences
--------------

Given a fasta file (let's say `unknown_seq.fasta`), you can find the taxonomy annotation of these
sequences using:
```
stag classify -d test_db.stagDB -i unknown_seq.fasta
```

The output is:
```
sequence	taxonomy
geneA	d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus
geneB	d__Bacteria
geneC	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria
```

You can either create a database (see [Create a database](https://github.com/AlessioMilanese/stag#create-a-database)), or use one that we already compiled:

| For | Taxonomy used | Run with | Link |
| --- | --- | --- | --- |
| full-length 16S gene | GTDB | `stag classify` | [link](https://www.embl.de/download/zeller/milanese/STAG/databases/16S/GTDB89_16S.stagDB) |

Create a database
--------------
You need three input:
1. a set of reference sequences to use to learn the taxonomy
2. a taxonomy file for the sequences in point 1
3. a hmm file for the sequences in point 1 (check also: [What to do if you don't have a hmm file](https://github.com/AlessioMilanese/stag/wiki/Create-hmm-file))

The reference sequences should be in fasta format ([example file](https://www.embl.de/download/zeller/milanese/STAG/databases/16S/files/GTDB_parsed.fa)):
```
>gene1
ATATGCATTTTACGATATGCA...
>gene2
GCATTATTTCAGGGCTAGGCA...
>gene3
CCGGATTGGGATCAAAAAGCG...
```

The taxonomy file should contain the same ids as the fasta file and the taxonomy
as a tab separated file (where the taxonomy is separated by ";"), like:
```
gene\tKingdom;Phylum;Class;...
```
Example ([example file](https://www.embl.de/download/zeller/milanese/STAG/databases/16S/files/GTDB_parsed.tax)):
```
gene1 d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus
gene2 d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Listeriaceae;g__Listeria;s__Listeria monocytogenes
gene3 d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus suis
```

Example HMM file: ([example file](https://www.embl.de/download/zeller/milanese/STAG/databases/16S/files/bacterial_model.hmm)).

To check that your files are correct, you can run:
```
stag check_input -i <fasta_seqs> -x <taxonomy_file> -a <hmmfile>
```

Once there are no errors, you can run:
```
stag train -i <fasta_seqs> -x <taxonomy_file> -a <hmmfile> -o test_db.stagDB
```
