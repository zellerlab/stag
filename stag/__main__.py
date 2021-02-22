#!/usr/bin/env python
import os
import sys
import argparse
import shlex
import shutil
import time
import subprocess
import glob
import tempfile
import errno
import tarfile

from . import __version__ as tool_version
from .helpers import bco, print_error, check_file_exists, check_file_doesnt_exists
import stag.align as align
import stag.create_db as create_db
import stag.classify as classify
import stag.check_create_db_input_files as check_create_db_input_files
import stag.correct_seq as correct_seq
import stag.unzip_db as unzip_db
import stag.classify_genome as classify_genome
import stag.train_genome as train_genome
import stag.convert_ali as convert_ali

def handle_error(error, help_f):
    help_f()
    print_error()
    print(error, file=sys.stderr)
    sys.exit(1)

# ------------------------------------------------------------------------------
#       print the help informations
# ------------------------------------------------------------------------------
class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)


def msg(name=None):
    str_msg = f'''
\00
{bco.Cyan}Program:{bco.ResetAll} stag - Supervised Taxonomic Assignment of marker Genes
{bco.Cyan}Version:{bco.ResetAll} '''+tool_version+f'''

{bco.Cyan}Usage:{bco.ResetAll} stag <command> [options]

{bco.Cyan}Command:{bco.ResetAll}
 {bco.LightGreen}-- Single gene{bco.ResetAll}
      {bco.LightBlue}train{bco.ResetAll}        Train a classifier and create a database
      {bco.LightBlue}classify{bco.ResetAll}     Taxonomically annotate a gene

      {bco.LightBlue}align{bco.ResetAll}        Align a sequence to a hmm or infernal model
      {bco.LightBlue}create_db{bco.ResetAll}    Create a database given the aligned sequences
      {bco.LightBlue}check_input{bco.ResetAll}  Check the input for the train command
      {bco.LightBlue}correct_seq{bco.ResetAll}  Correct sequences that are in wrong orientation
      {bco.LightBlue}convert_ali{bco.ResetAll}  Convert between 1-hot-encoding and fasta, and vice versa
      {bco.LightBlue}unzip_db{bco.ResetAll}     Create a directory with the content of a database

 {bco.LightGreen}-- Genome{bco.ResetAll}
      {bco.LightBlue}train_genome{bco.ResetAll}     Merge classifiers of single genes
      {bco.LightBlue}classify_genome{bco.ResetAll}  Taxonomically annotate a genome (predict genes, extract
                       the database marker genes and classify them)

Type stag <command> to print the help for a specific command
        '''
    return str_msg

# ------------------------------------------------------------------------------
def print_menu_align():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} align {bco.LightBlue}-i{bco.ResetAll} <fasta_seqs> {bco.LightBlue}-a{bco.ResetAll} <hmmfile> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  FILE  sequences to be aligned (fasta format) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-p{bco.ResetAll}  FILE  protein sequences, corresponding to -i {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-a{bco.ResetAll}  FILE  hmmfile or cmfile to use as template for the alignment {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  FILE  output file name {bco.LightMagenta}[stdout]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-c{bco.ResetAll}        set if you are using a cmfile\n")
    sys.stderr.write(f"  {bco.LightBlue}-m{bco.ResetAll}  INT   threshold for the number of features per sequence (percentage) {bco.LightMagenta}[0]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-t{bco.ResetAll}  INT   number of threads {bco.LightMagenta}[1]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")
    sys.stderr.write(f"{bco.Cyan}Note:{bco.ResetAll} if -p is provided, then the alignment will be done at the level\nof the proteins and then converted to gene alignment (from -i input).\nThe order of the sequences in -i and -p should be the same.\n\n")
# ------------------------------------------------------------------------------
def print_menu_create_db():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} create_db {bco.LightBlue}-s{bco.ResetAll} <aligned_file> {bco.LightBlue}-x{bco.ResetAll} <taxonomy_file>\n")
    sys.stderr.write(f"                     {bco.LightBlue}-a{bco.ResetAll} <hmmfile> {bco.LightBlue}-o{bco.ResetAll} <output_DB> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-s{bco.ResetAll}  FILE  file with 1-hot encoding MSA (result from stag align) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-a{bco.ResetAll}  FILE  hmmfile or cmfile to used as template for the alignment {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-c{bco.ResetAll}        set if you are using a cmfile\n")
    sys.stderr.write(f"  {bco.LightBlue}-x{bco.ResetAll}  FILE  taxonomy file (tab separated) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  FILE  output file name (HDF5 format) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-f{bco.ResetAll}        force to rewrite output file\n")
    sys.stderr.write(f"  {bco.LightBlue}-C{bco.ResetAll}  FILE  save intermediate cross validation results {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-p{bco.ResetAll}  FILE  protein sequences, if they were used for the alignment {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-e{bco.ResetAll}  STR   penalty for the logistic regression {bco.LightMagenta}[\"l1\"]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-E{bco.ResetAll}  STR   solver for the logistic regrssion {bco.LightMagenta}[\"liblinear\"]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")
# ------------------------------------------------------------------------------
def print_menu_classify():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} classify {bco.LightBlue}-d{bco.ResetAll} <database> [{bco.LightBlue}-i{bco.ResetAll}/{bco.LightBlue}-s{bco.ResetAll}] <seq_file> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-d{bco.ResetAll}  FILE  database created with create_db or train {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  FILE  sequences to taxonomically annotate (fasta format) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-s{bco.ResetAll}  FILE  aligned sequences, can be provided instead of -i {bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-p{bco.ResetAll}  FILE  protein sequences, corresponding to -i {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-S{bco.ResetAll}  FILE  save intermediate alignment file {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  FILE  output file name {bco.LightMagenta}[stdout]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-l{bco.ResetAll}        long output (with more information about the classification)\n")
    sys.stderr.write(f"  {bco.LightBlue}-m{bco.ResetAll}  INT   threshold for the number of features per sequence (percentage) {bco.LightMagenta}[0]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")
# ------------------------------------------------------------------------------
def print_menu_check_input():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} check_input {bco.LightBlue}-i{bco.ResetAll} <fasta_seqs> {bco.LightBlue}-x{bco.ResetAll} <taxonomy_file>\n")
    sys.stderr.write(f"                       {bco.LightBlue}-a{bco.ResetAll} <hmmfile> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  FILE  sequences to be aligned (fasta format) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-p{bco.ResetAll}  FILE  protein sequences, corresponding to -i {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-a{bco.ResetAll}  FILE  hmmfile or cmfile to used as template for the alignment {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-c{bco.ResetAll}        set if you are using a cmfile\n")
    sys.stderr.write(f"  {bco.LightBlue}-x{bco.ResetAll}  FILE  taxonomy file (tab separated) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-w{bco.ResetAll}  FILE  save warning messages to a file {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")
# ------------------------------------------------------------------------------
def print_menu_train():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} train {bco.LightBlue}-i{bco.ResetAll} <fasta_seqs> {bco.LightBlue}-x{bco.ResetAll} <taxonomy_file>\n")
    sys.stderr.write(f"                 {bco.LightBlue}-a{bco.ResetAll} <hmmfile> {bco.LightBlue}-o{bco.ResetAll} <output_DB> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  FILE  sequences to be aligned (fasta format) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-p{bco.ResetAll}  FILE  protein sequences, corresponding to -i {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-a{bco.ResetAll}  FILE  hmmfile or cmfile to used as template for the alignment {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-c{bco.ResetAll}        set if you are using a cmfile\n")
    sys.stderr.write(f"  {bco.LightBlue}-x{bco.ResetAll}  FILE  taxonomy file (tab separated) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  FILE  output file name (HDF5 format) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-f{bco.ResetAll}        force to rewrite output file\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-S{bco.ResetAll}  FILE  save intermediate alignment file {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-C{bco.ResetAll}  FILE  save intermediate cross validation results {bco.LightMagenta}[None]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-t{bco.ResetAll}  INT   number of threads {bco.LightMagenta}[1]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-m{bco.ResetAll}  INT   threshold for the number of features per sequence (percentage) {bco.LightMagenta}[0]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-e{bco.ResetAll}  STR   penalty for the logistic regression {bco.LightMagenta}[\"l1\"]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-E{bco.ResetAll}  STR   solver for the logistic regrssion {bco.LightMagenta}[\"liblinear\"]{bco.ResetAll}\n\n")
    sys.stderr.write(f"{bco.Cyan}Note:{bco.ResetAll} if -p is provided, then the alignment will be done at the level\nof the proteins and then converted to gene alignment (from -i input).\nThe order of the sequences in -i and -p should be the same.\n\n")
# ------------------------------------------------------------------------------
def print_menu_correct_seq():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} correct_seq {bco.LightBlue}-i{bco.ResetAll} <fasta_seqs> {bco.LightBlue}-a{bco.ResetAll} <hmmfile> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  FILE  sequences to be aligned (fasta format) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-a{bco.ResetAll}  FILE  hmmfile or cmfile to use as template for the alignment {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  FILE  output file name {bco.LightMagenta}[stdout]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-c{bco.ResetAll}        set if you are using a cmfile\n")
    sys.stderr.write(f"  {bco.LightBlue}-m{bco.ResetAll}  INT   threshold for the number of features per sequence (percentage) {bco.LightMagenta}[5]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-t{bco.ResetAll}  INT   number of threads {bco.LightMagenta}[1]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")
# ------------------------------------------------------------------------------
def print_menu_train_genome():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} train_genome {bco.LightBlue}-i{bco.ResetAll} <list_gene_DBs> {bco.LightBlue}-T{bco.ResetAll} <gene_thresholds>\n                         {bco.LightBlue}-o{bco.ResetAll} <output_DB> {bco.LightBlue}-C{bco.ResetAll} <concat_genes_DB> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  LIST  list of single gene databases to use (comma separated) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-T{bco.ResetAll}  FILE  hmm treshold for selecting the genes {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-C{bco.ResetAll}  FILE  stag database for the concatenated genes{bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  FILE  output file name (HDF5 format) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-t{bco.ResetAll}  INT   number of threads {bco.LightMagenta}[1]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")
# ------------------------------------------------------------------------------
def print_menu_classify_genome():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} classify_genome {bco.LightBlue}-d{bco.ResetAll} <genome_database> {bco.LightBlue}-o{bco.ResetAll} res_dir\n")
    sys.stderr.write(f"                            [{bco.LightBlue}-i{bco.ResetAll} <fasta_seq>/{bco.LightBlue}-D{bco.ResetAll} <directory>] [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-d{bco.ResetAll}  FILE   database created with train_genome {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  FILE   genome fasta file\n")
    sys.stderr.write(f"  {bco.LightBlue}-D{bco.ResetAll}  DIR    directory containing genome fasta files (only fasta\n             files will be used)\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  DIR    output directory {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-l{bco.ResetAll}         long output (with more information about the classification) {bco.LightMagenta}[False]\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT    verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-r{bco.ResetAll}         use all genes above the filter {bco.LightMagenta}[False]{bco.ResetAll}\n\n")
# ------------------------------------------------------------------------------
def print_menu_convert_ali():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} convert_ali {bco.LightBlue}-i{bco.ResetAll} <file_in> {bco.LightBlue}-o{bco.ResetAll} <file_out> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  FILE   Input file, either a 1-hot-encoding created by stag align,\n")
    sys.stderr.write(f"             or a fasta file of aligned sequences created by hmmalign.\n")
    sys.stderr.write(f"             The input type is detected automatically {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  FILE   A 1-hot-encoding if the input was fasta,\n")
    sys.stderr.write(f"             or a fasta file if the input was 1-hot-encoding {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT    verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")
# ------------------------------------------------------------------------------
def print_menu_unzip_db():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} unzip_db {bco.LightBlue}-d{bco.ResetAll} <database> {bco.LightBlue}-o{bco.ResetAll} <dir_out> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-d{bco.ResetAll}  FILE  database created with create_db or train {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  DIR   create a dir with the unzipped database {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def main(argv=None):

    parser = argparse.ArgumentParser(usage=msg(), formatter_class=CapitalisedHelpFormatter,add_help=False)
    parser.add_argument('command', action="store", default=None, help='mode to run stag',
                        choices=['align','train','classify','create_db','check_input','correct_seq','train_genome',
                                 'classify_genome','test','convert_ali',"unzip_db"])
    parser.add_argument('-o', action="store", dest='output', default=None, help='name of output file')
    parser.add_argument('-t', type=int, action="store", dest='threads', default=1, help='Number of threads to be used.')
    parser.add_argument('-v', action='store', type=int, default=3, dest='verbose', help='Verbose levels', choices=list(range(1,5)))
    parser.add_argument('-c', action='store_true', dest='use_cm_align', help='Set if you want to use cmalign isntead of hmmalign')
    parser.add_argument('-s', action="store", default=None,dest='aligned_sequences', help='sequences that needs to be aligned')
    parser.add_argument('-a', action="store", default=None,dest='template_al', help='alignment template')
    parser.add_argument('-x', action="store", default=None,dest='taxonomy', help='taxonomy file path')
    parser.add_argument('-f', action='store_true', dest='force_rewrite', help='Set if you want to rewrite the file, even if it exists')
    parser.add_argument('-i', action="store", dest='fasta_input', default=None, help='input fasta sequences')
    parser.add_argument('-p', action="store", dest='protein_fasta_input', default=None, help='input fasta sequences, in protein format. Corresponding to the -i sequences')
    parser.add_argument('-w', action="store", dest='warning_file_check_input', default=None, help='for check_input there can be quite some warning messages. Use -w to save them to a file')
    parser.add_argument('-d', action="store", dest='database', default=None, help='file containing the database')
    parser.add_argument('-S', action="store", dest='intermediate_al', default=None, help='name of the file for the intermediate alignment')
    parser.add_argument('-C', action="store", dest='intermediate_cross_val', default=None, help='name of the file for the intermediate cross validation results')
    parser.add_argument('-m', action='store', type=int, default=None, dest='min_perc_state', help='Minimum number of mapping states, i.e. how many features of the classifier we cover')
    parser.add_argument('-l', action='store_true', dest='long_out', help='Print more columns for the classification pipeline')
    parser.add_argument('-r', action='store_true', dest='keep_all_genes', help='keep all genes when doing the classification of genomes')
    parser.add_argument('-D', action="store", dest='dir_input', default=None, help='input directory')
    parser.add_argument('-T', action="store", dest='file_thresholds', default=None, help='file with the thresholds for the genes in the genome classifier') # basically the minimum score required
    parser.add_argument('-e', action="store", default="l1", dest='penalty_logistic', help='penalty for the logistic regression',choices=['l1','l2','none'])
    parser.add_argument('-E', action="store", default="liblinear", dest='solver_logistic', help='solver for the logistic regression',choices=['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'])

    parser.add_argument('--version', action='version', version='%(prog)s {0} on python {1}'.format(tool_version, sys.version.split()[0]))

    args = parser.parse_args()

    # --------------------------------------------------------------------------
    # TEST routine
    # --------------------------------------------------------------------------
    if args.command == 'test':
        popenCMD = "stag_test"
        child = subprocess.Popen(popenCMD)
        child.communicate()
        rc = child.wait()
        return(rc)

    # --------------------------------------------------------------------------
    # set defaults for the parameters
    # --------------------------------------------------------------------------
    if args.command == 'correct_seq':
        if (args.min_perc_state is None): args.min_perc_state = 5
    else:
        if (args.min_perc_state is None): args.min_perc_state = 0

    if args.threads < 1:
        handle_error("number of threads (-t) is less than 1", None)
    if args.min_perc_state < 0 or args.min_perc_state > 100:
        handle_error("-m should be between 0 and 100. It represents the percentage\n"
                     "of internal states covered by the sequence (i.e. the number of features).", None)


    # --------------------------------------------------------------------------
    # ALIGN routine
    # --------------------------------------------------------------------------
    error = ""
    if args.command == 'align':
        # check that '-i' and '-a' have been provided
        if not args.fasta_input:
            error = "missing <seqfile> (-i)"
        elif not args.template_al:
            error = "missing <hmmfile>/<cmfile> (-a)"

        if error:
            handle_error(error, print_menu_align)

        # check that '-s' and '-a' are files
        check_file_exists(args.fasta_input, isfasta = True)
        check_file_exists(args.template_al, isfasta = False)

        # if -p is provided, then check that it is a fasta file
        if args.protein_fasta_input:
            check_file_exists(args.protein_fasta_input, isfasta = True)

        # call the function
        if args.output is None:
            for ali in align.align_generator(args.fasta_input, args.protein_fasta_input, args.template_al,
                                             args.use_cm_align, args.threads, args.verbose, False, args.min_perc_state):
                print(ali)
        else:
            align.align_file(args.fasta_input, args.protein_fasta_input, args.template_al, args.use_cm_align, args.threads,
                             args.verbose, args.output, args.min_perc_state)

    # --------------------------------------------------------------------------
    # CREATE_DB routine
    # --------------------------------------------------------------------------
    elif args.command == 'create_db':

        if not args.aligned_sequences:
            # check that '-s' has been provided (alignment produced by stag align)
            error = "missing <aligned_file> (-s)"
        elif not args.taxonomy:
            # check that '-x' has been provided (taxonomy file)
            error = "missing <taxonomy_file> (-x)"
        elif not args.template_al:
            # check that the hmm file is provided
            error = "missing <hmmfile>/<cmfile> (-a)"
        elif not args.output:
            # check that output is set
            error = "missing <output_DB> (-o)"

        if error:
            handle_error(error, print_menu_create_db)

        # check that '-s' and '-a' are files
        check_file_exists(args.aligned_sequences, isfasta = False)
        check_file_exists(args.taxonomy, isfasta = False)
        check_file_exists(args.template_al, isfasta = False)

        if not args.force_rewrite:
            check_file_doesnt_exists(args.output)

        # call the function to create the database
        create_db.create_db(args.aligned_sequences, args.taxonomy, args.verbose, args.output,
                            args.use_cm_align, args.intermediate_cross_val, tool_version,
                            args.penalty_logistic, args.solver_logistic,
                            hmm_file_path=args.template_al, protein_fasta_input=args.protein_fasta_input, procs=args.threads)

    # --------------------------------------------------------------------------
    # TRAIN routine
    # --------------------------------------------------------------------------
    elif args.command == 'train':
        # check that '-i' and '-a' have been provided
        if not args.fasta_input:
            error = "missing <seqfile> (-i)"
        elif not args.template_al:
            error = "missing <hmmfile>/<cmfile> (-a)"
        elif not args.taxonomy:
            # check that '-x' has been provided (taxonomy file)
            error = "missing <taxonomy_file> (-x)"
        elif not args.output:
            # check that output is set
            error = "missing <output_DB> (-o)"

        if error:
            handle_error(error, print_menu_train)

        # check that '-s' and '-a' are files
        check_file_exists(args.fasta_input, isfasta = True)
        check_file_exists(args.template_al, isfasta = False)

        # if -p is provided, then check that it is a fasta file
        if args.protein_fasta_input:
            check_file_exists(args.protein_fasta_input, isfasta = True)

        check_file_exists(args.taxonomy, isfasta = False)

        if not args.force_rewrite:
            check_file_doesnt_exists(args.output)

        # FIRST: ALIGN ---------------------------------------------------------
        # we create a temporary file that will contain the alignments
        al_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
        os.chmod(al_file.name, 0o644)
        # call the function
        align.align_file(args.fasta_input, args.protein_fasta_input, args.template_al, args.use_cm_align,
                         args.threads, args.verbose, al_file.name, args.min_perc_state)

        # SECOND: CREATE_DB ----------------------------------------------------
        # call the function to create the database
        create_db.create_db(al_file.name, args.taxonomy, args.verbose, args.output, args.use_cm_align,
                            args.intermediate_cross_val, tool_version,
                            args.penalty_logistic, args.solver_logistic,
                            hmm_file_path=args.template_al, protein_fasta_input=args.protein_fasta_input, procs=args.threads)

        # what to do with intermediate alignment -------------------------------
        if not args.intermediate_al:
            # remove it
            os.remove(al_file.name)
        else:
            # save it
            shutil.move(al_file.name, args.intermediate_al)


    # --------------------------------------------------------------------------
    # CLASSIFY routine
    # --------------------------------------------------------------------------
    elif args.command == 'classify':
        # check that '-i' has been provided (alignment produced by stag align)
        if not args.fasta_input and not args.aligned_sequences:
            error = "missing <fasta_seqs> (-i) or <aligned_seq> (-s)"
        elif not args.database:
            # check that '-d' has been provided (taxonomy file)
            error = "missing <database> (-d)"

        if error:
            handle_error(error, print_menu_classify)

        # check that they are files
        if args.fasta_input:
            check_file_exists(args.fasta_input, isfasta = True)
        check_file_exists(args.database, isfasta = False)
        # if -p is provided, then check that it is a fasta file
        if args.protein_fasta_input:
            check_file_exists(args.protein_fasta_input, isfasta = True)
        # if -S is provided, we remove the file if it exists, since in the
        # function it appends only
        if args.intermediate_al:
            if os.path.isfile(args.intermediate_al):
                os.remove(args.intermediate_al)


        # call the function
        classify.classify(args.database, fasta_input=args.fasta_input, protein_fasta_input=args.protein_fasta_input, 
                          verbose=args.verbose, threads=args.threads, output=args.output, long_out=args.long_out, 
                          current_tool_version=tool_version, aligned_sequences=args.aligned_sequences,
                          save_ali_to_file=args.intermediate_al, min_perc_state=args.min_perc_state)

    # --------------------------------------------------------------------------
    # CHECK_INPUT routine
    # --------------------------------------------------------------------------
    elif args.command == 'check_input':
        if not args.fasta_input:
            error = "missing <fasta_seqs> (-i)"
        elif not args.taxonomy:
            error = "missing <taxonomy_file> (-x)"
        elif not args.template_al:
            error = "missing <hmmfile>/<cmfile> (-a)"

        if error:
            handle_error(error, print_menu_check_input)

        check_create_db_input_files.check_input_files(args.fasta_input, args.protein_fasta_input, args.taxonomy,
                                                      args.template_al, args.use_cm_align, args.warning_file_check_input)

    # --------------------------------------------------------------------------
    # CORRECT_SEQ routine
    # --------------------------------------------------------------------------
    # check if the sequences are in correct orientation, if they are not, then
    # take reverse complement. Save to -o all the seqeunces is correct order
    elif args.command == 'correct_seq':
        # check that '-i' and '-a' have been provided
        if not args.fasta_input:
            error = "missing <seqfile> (-i)"
        elif not args.template_al:
            error = "missing <hmmfile>/<cmfile> (-a)"

        if error:
            handle_error(error, print_menu_correct_seq)

        # check that '-s' and '-a' are files
        check_file_exists(args.fasta_input, isfasta = True)
        check_file_exists(args.template_al, isfasta = False)

        # call the function
        correct_seq.correct_seq(args.fasta_input, args.template_al, args.use_cm_align, args.threads, args.verbose,
                                args.min_perc_state, args.output)

    # --------------------------------------------------------------------------
    # CONVERT_ALI routine
    # --------------------------------------------------------------------------
    elif args.command == 'convert_ali':
        # check that '-i' and '-o' have been provided
        if not args.fasta_input:
            error = "missing <file_in> (-i)"
        elif not args.output:
            error = "missing <file_out> (-o)"

        if error:
            handle_error(error, print_menu_convert_ali)

        # check that '-i' is a file
        check_file_exists(args.fasta_input, isfasta = False)

        # call function
        convert_ali.convert_ali(args.fasta_input, args.output, args.verbose)


    # --------------------------------------------------------------------------
    # UNZIP_db routine
    # --------------------------------------------------------------------------
    elif args.command == 'unzip_db':
        # check that '-d' and '-o' have been provided
        if not args.database:
            error = "missing <database> (-d)"
        elif not args.output:
            error = "missing <dir_out> (-o)" 

        if error:
            handle_error(error, print_menu_unzip_db)

        # check that '-d' is a file
        check_file_exists(args.database, isfasta = False)

        # call function
        unzip_db.unzip_db(args.database, args.verbose, args.output)

    # --------------------------------------------------------------------------
    # TRAIN_GENOME routine
    # --------------------------------------------------------------------------
    # We want to have a database for classifying genomes. It will contains the
    # calssifiers for the seingle genes
    # Input: single gene databases
    # Output: a database file (hdf5) that can be used by stag classify_genome
    elif args.command == 'train_genome':
        # check that parameters are set ----------------------------------------
        if not args.output:
            error = "missing <output_DB> (-o)"
        elif not args.fasta_input:
            error = "missing <list_gene_DBs> (-i)"
        elif not args.file_thresholds:
            error = "missing <gene_thresholds> (-T)"
        elif not args.intermediate_cross_val:
            error = "missing <concat_genes_DB> (-C)"

        if error:
            handle_error(error, print_menu_train_genome)

        # call the function
        train_genome.train_genome(args.output, args.fasta_input, args.file_thresholds,
                                  args.threads, args.verbose, args.intermediate_cross_val)

    # --------------------------------------------------------------------------
    # CLASSIFY_GENOME routine
    # --------------------------------------------------------------------------
    if args.command == 'classify_genome':
        # check input
        if not args.database:
            error = "missing <database> (-d)"
        elif not args.fasta_input and not args.dir_input:
            error = "you need to provide at least -i or -D."
        elif args.fasta_input and args.dir_input:
            error = "you need to provide -i or -D, not both."
        elif args.dir_input and not os.path.isdir(args.dir_input):
            error = "-D is not a directory."
        elif not args.output:
            # check that output dir is defined
            error = "missing output directory (-o)"

        if error:
            handle_error(error, print_menu_classify_genome)

        # find files to classify
        list_files = list()
        if args.fasta_input:
            check_file_exists(args.fasta_input, isfasta = True)
            list_files.append(args.fasta_input)
        else:
            for f in os.listdir(args.dir_input):
                f = os.path.join(args.dir_input, f)
                try:
                    if os.path.isfile(f) and open(f).read(1)[0] == ">":
                        list_files.append(f)
                except Exception as e:
                    if args.verbose > 1:
                        sys.stderr.write("[W::main] Warning: ")
                        sys.stderr.write("Cannot open file: {}\n".format(f))

            if not list_files:
                handle_error("no fasta files found in the directory.", None)
            sys.stderr.write(" Found "+str(len(list_files))+" fasta files\n")

        if os.path.isdir(args.output):
            if args.force_rewrite:
                shutil.rmtree(args.output)
            else:
                handle_error("output directory (-o) exists already.", None)

        # create output dir
        try:
            pathlib.Path(args.output).mkdir(exist_ok=True, parents=True)
        except:
            handle_error("creating the output directory (-o).", None)

        if list_files:
            from stag.classify_genome import validate_genome_files
            validate_genome_files(list_files)

        # call the function
        classify_genome.classify_genome(args.database, list_files, args.verbose, args.threads, args.output,
                                        args.long_out, tool_version, args.keep_all_genes)


    return None        # success

#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
    main()
