<img src="https://github.com/AlessioMilanese/stag/blob/master/pics/stag_logo.png" width="450">

This tool is design to classify metagenomic sequences (marker genes, genomes and amplicon reads) using a Hierarchical Taxonomic Classifier.


Pre-requisites
--------------

The stag classifier requires:
* Python 3 (or higher)
* HMMER3 (or Infernal)
* Easel ([link](https://github.com/EddyRivasLab/easel))
* python library:
  * numpy
  * pandas
  * sklearn
  * h5py

In order to use the command ```correct_seq``` you need:
* [seqtk](https://github.com/lh3/seqtk)

If you have [conda](https://conda.io/docs/), you can install the dependencies and create an enviroment:
```bash
conda env create -f conda_env_stag.yaml
conda activate stag
```

Installation
--------------
```bash
git clone https://github.com/AlessioMilanese/stag.git
cd stag
export PATH=`pwd`:$PATH
```

Note: in the following examples we assume that the python script ```stag``` is in the system path.


First: Create a database
--------------
You need three input:
1. a set of reference sequences to use to learn the taxonomy
2. a taxonomy file for the sequences in point 1
3. a hmm file for the sequences in point 1

The reference sequences should be in fasta format:
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
Example:
```
gene1 d__Bacteria;p__Firmicutes;c__Bacilli;o__Staphylococcales;f__Staphylococcaceae;g__Staphylococcus;s__Staphylococcus aureus
gene2 d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Listeriaceae;g__Listeria;s__Listeria monocytogenes
gene3 d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus suis
```

To check that your files are correct, you can run:
```
stag check_input -i <fasta_seqs> -x <taxonomy_file> -a <hmmfile>
```

Once there are no errors, you can run:
```
stag train -i <fasta_seqs> -x <taxonomy_file> -a <hmmfile> -o test_db.htcDB
```


Second: Taxonomically annotate unknown sequences
--------------

Given a fasta file (let's say `unknown_seq.fasta`), you can find the taxonomy annotation of these
sequences using:
```
stag classify -d test_db.htcDB -i unknown_seq.fasta
```
