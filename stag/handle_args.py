# pylint: disable=C0301

import argparse
import sys

from . import __version__ as tool_version
from .helpers import bco, print_error


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


class CapitalisedHelpFormatter(argparse.HelpFormatter):
    def add_usage(self, usage, actions, groups, prefix=None):
        if prefix is None:
            prefix = ''
        return super(CapitalisedHelpFormatter, self).add_usage(usage, actions, groups, prefix)


def get_args():
    parser = argparse.ArgumentParser(usage=msg(), formatter_class=CapitalisedHelpFormatter, add_help=False)
    parser.add_argument(
        'command', action="store", default=None, help='mode to run stag',
        choices=[
            'align', 'train', 'classify', 'create_db', 'check_input', 'correct_seq', 'train_genome',
            'classify_genome', 'test', 'convert_ali', 'unzip_db'
        ]
    )
    parser.add_argument('-o', action="store", dest='output', default=None, help='name of output file')
    parser.add_argument('-t', type=int, action="store", dest='threads', default=1, help='Number of threads to be used.')
    parser.add_argument('-v', action='store', type=int, default=3, dest='verbose', help='Verbose levels', choices=list(range(1, 5)))
    parser.add_argument('-c', action='store_true', dest='use_cm_align', help='Set if you want to use cmalign isntead of hmmalign')
    parser.add_argument('-s', action="store", default=None, dest='aligned_sequences', help='sequences that needs to be aligned')
    parser.add_argument('-a', action="store", default=None, dest='template_al', help='alignment template')
    parser.add_argument('-x', action="store", default=None, dest='taxonomy', help='taxonomy file path')
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
    parser.add_argument('-T', action="store", dest='file_thresholds', default=None, help='file with the thresholds for the genes in the genome classifier')  # basically the minimum score required
    parser.add_argument('-e', action="store", default="l1", dest='penalty_logistic', help='penalty for the logistic regression', choices=['l1', 'l2', 'none'])
    parser.add_argument('-E', action="store", default="liblinear", dest='solver_logistic', help='solver for the logistic regression', choices=['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'])
    parser.add_argument('-G', action="store", dest="marker_genes", default=None, help="Set of identified marker genes in lieu of a genomic sequence")
    parser.add_argument('-N', action="store", dest="solver_iterations", default=5000, type=int, help="Increase this parameter if the output displays a ConvergenceWarning.")

    parser.add_argument('--version', action='version', version='%(prog)s {0} on python {1}'.format(tool_version, sys.version.split()[0]))

    return parser.parse_args()


def handle_error(error, help_f=None):
    if help_f:
        help_f()
    print_error()
    print(error, file=sys.stderr)
    sys.exit(1)


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
    sys.stderr.write(f"  {bco.LightBlue}-E{bco.ResetAll}  STR   solver for the logistic regression {bco.LightMagenta}[\"liblinear\"]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-N{bco.ResetAll}  INT   solver_iterations: increase this parameter if the output displays a ConvergenceWarning. (default=5000)\n")
    sys.stderr.write(f"  {bco.LightBlue}-t{bco.ResetAll}  INT   number of threads {bco.LightMagenta}[1]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")


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
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-e{bco.ResetAll}  STR   penalty for the logistic regression {bco.LightMagenta}[\"l1\"]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-E{bco.ResetAll}  STR   solver for the logistic regression {bco.LightMagenta}[\"liblinear\"]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-N{bco.ResetAll}  INT   solver_iterations: increase this parameter if the output displays a ConvergenceWarning. (default=5000)\n\n")
    sys.stderr.write(f"{bco.Cyan}Note:{bco.ResetAll} if -p is provided, then the alignment will be done at the level\nof the proteins and then converted to gene alignment (from -i input).\nThe order of the sequences in -i and -p should be the same.\n\n")


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


def print_menu_train_genome():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} train_genome {bco.LightBlue}-i{bco.ResetAll} <list_gene_DBs> {bco.LightBlue}-T{bco.ResetAll} <gene_thresholds>\n                         {bco.LightBlue}-o{bco.ResetAll} <output_DB> {bco.LightBlue}-C{bco.ResetAll} <concat_genes_DB> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  LIST  list of single gene databases to use (comma separated) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-T{bco.ResetAll}  FILE  hmm treshold for selecting the genes {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-C{bco.ResetAll}  FILE  stag database for the concatenated genes{bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  FILE  output file name (HDF5 format) {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-t{bco.ResetAll}  INT   number of threads {bco.LightMagenta}[1]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}")
    sys.stderr.write(f"  {bco.LightBlue}-N{bco.ResetAll}  INT   solver_iterations: increase this parameter if the output displays a ConvergenceWarning. (default=5000)\n\n")


def print_menu_classify_genome():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} classify_genome {bco.LightBlue}-d{bco.ResetAll} <genome_database> {bco.LightBlue}-o{bco.ResetAll} res_dir\n")
    sys.stderr.write(f"                            [{bco.LightBlue}-i{bco.ResetAll} <fasta_seq>/{bco.LightBlue}-D{bco.ResetAll} <directory>/{bco.LightBlue}-G{bco.ResetAll} <markers.json>] [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-d{bco.ResetAll}  FILE   database created with train_genome {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  FILE   genome fasta file\n")
    sys.stderr.write(f"  {bco.LightBlue}-D{bco.ResetAll}  DIR    directory containing genome fasta files (only fasta\n             files will be used)\n")
    sys.stderr.write(f"  {bco.LightBlue}-G{bco.ResetAll}  FILE   json file pointing at a marker gene set (in lieu of a full genome)\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  DIR    output directory {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-l{bco.ResetAll}         long output (with more information about the classification) {bco.LightMagenta}[False]\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT    verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-r{bco.ResetAll}         use all genes above the filter {bco.LightMagenta}[False]{bco.ResetAll}\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-t{bco.ResetAll}  INT   number of threads {bco.LightMagenta}[1]{bco.ResetAll}\n")


def print_menu_convert_ali():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} convert_ali {bco.LightBlue}-i{bco.ResetAll} <file_in> {bco.LightBlue}-o{bco.ResetAll} <file_out> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-i{bco.ResetAll}  FILE   Input file, either a 1-hot-encoding created by stag align,\n")
    sys.stderr.write(f"             or a fasta file of aligned sequences created by hmmalign.\n")
    sys.stderr.write(f"             The input type is detected automatically {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  FILE   A 1-hot-encoding if the input was fasta,\n")
    sys.stderr.write(f"             or a fasta file if the input was 1-hot-encoding {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT    verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")


def print_menu_unzip_db():
    sys.stderr.write("\n")
    sys.stderr.write(f"{bco.Cyan}Usage:{bco.ResetAll} {bco.Green}stag{bco.ResetAll} unzip_db {bco.LightBlue}-d{bco.ResetAll} <database> {bco.LightBlue}-o{bco.ResetAll} <dir_out> [options]\n\n")
    sys.stderr.write(f"  {bco.LightBlue}-d{bco.ResetAll}  FILE  database created with create_db or train {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-o{bco.ResetAll}  DIR   create a dir with the unzipped database {bco.LightMagenta}[required]{bco.ResetAll}\n")
    sys.stderr.write(f"  {bco.LightBlue}-v{bco.ResetAll}  INT   verbose level: 1=error, 2=warning, 3=message, 4+=debugging {bco.LightMagenta}[3]{bco.ResetAll}\n\n")
