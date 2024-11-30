<img align="left" src="https://github.com/zellerlab/stag/blob/master/pics/stag_logo.png" width="250">

<br/><br/>
This tool is design to classify metagenomic sequences (marker genes, genomes and amplicon reads) using a Hierarchical Taxonomic Classifier.

Please check also the [wiki](https://github.com/zellerlab/stag/wiki) for more information.

<br/><br/><br/><br/>

Dependencies
--------------

The stag classifier requires:
* Python 3.7 (or higher)
* HMMER3 (or Infernal)
* Easel ([link](https://github.com/EddyRivasLab/easel))
* [seqtk](https://github.com/lh3/seqtk)
* [prodigal](https://github.com/hyattpd/Prodigal) (to predict genes in genomes)
* python library:
  * numpy
  * pandas
  * sklearn
  * h5py = 2.10.0

If you have [conda](https://conda.io/docs/), you can install all the dependencies in `conda_env_stag.yaml`.
See [Installation wiki](https://github.com/zellerlab/stag/wiki/Installation) for more info.


Installation
--------------
```bash
git clone https://github.com/zellerlab/stag.git
cd stag
# if environment is needed
conda env create -f conda_env_stag.yaml
python setup.py bdist_wheel
pip install --no-deps --force-reinstall dist/*.whl
```

Note: in the following examples we assume that the python script `stag` is in the system path.

Execution
---------

```bash
# if environment was installed
conda activate stag
# test the installation
stag test
```


Taxonomically annotate gene sequences
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

You can either create a database (see [Create a database](https://github.com/zellerlab/stag/wiki/Build-STAG-database-for-genes)), or use one that we already compiled:

- [For marker genes](https://github.com/zellerlab/stag/wiki/Classify-genes)
- [For 16S amplicon data](https://github.com/zellerlab/stag/wiki/16S-amplicon-databases)



Taxonomically annotate genomes
--------------

Given a fasta file (let's say `unknown_genome.fasta`), you can find the taxonomy annotation of this genome with:
```
stag classify_genome -i unknown_genome.fasta -d gtdb_30.stagDB -o res_dir
```

The output is saved in the directory `res_dir`. Inside you will find the file `genome_annotation` with the annotation
in the same format as in the gene classification. More information on the other files can be found [here](https://github.com/zellerlab/stag/wiki/Classify-genomes).

To classify multiple genomes, you can use:
```
stag classify_genome -D all/genomes/dir -d gtdb_30.stagDB -o res_dir
```
Where `all/genomes/dir` is a directory, and all fasta files inside the directory will be classified.

Finally, you can find some databases to classify genomes (`gtdb_30.stagDB` in the examples) [here](https://github.com/zellerlab/stag/wiki/Genomes-databases).




Schematic depiction of the STAG workflow. 
--------------

(a) Example taxonomic tree alongside
thirteen (partial) 16S sequences in a multiple sequence alignment (MSA). Four positions in the
MSA are highlighted for the information they contain to distinguish the different clades shown.
For example, at position 455, a ‘G’ distinguishes Enterobacteriaceae from Erwiniaceae. To
leverage this information, STAG trains one LASSO logistic regression classifier for each node in
the tree; coefficients corresponding to aligned bases are shown in b, c and d. For example, a ‘C’ at
position 648 is learnt to facilitate discrimination of Escherichia coli from Escherichia albertii.

To annotate a new sequence, it is first aligned to the MSA constructed during training. Second, the
sequence is classified along the tree, following the path with the highest posterior probabilities
(as returned by the node classifiers). Finally, the taxonomic lineage of the new sequence is
inferred from the probabilities accrued in the previous step; this in particular entails the decision
of which ranks not to assign.

<img src="https://github.com/zellerlab/stag/blob/master/pics/fig1.png">

