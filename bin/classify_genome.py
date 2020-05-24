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

# ==============================================================================
# UNZIP THE DATABASE
# ==============================================================================
def load_genome_DB(database, tool_version, verbose):
    dirpath = tempfile.mkdtemp()
    shutil.unpack_archive(database, dirpath,"gztar")
    list_files = [f for f in os.listdir(dirpath) if os.path.isfile(os.path.join(dirpath, f))]
    return list_files,dirpath

#===============================================================================
#                                      MAIN
#===============================================================================
def classify_genome(database, fasta_input, verbose, threads, output, long_out, tool_version):
    # FIRST: unzip the database
    database_files, temp_dir = load_genome_DB(database, tool_version, verbose)

    # SECOND: run prodigal on the fasta genome

    # THIRD: find the marker genes from the predicted genes

    # FOURTH: classify the marker genes

    # we remove the temp dir
    shutil.rmtree(temp_dir)

    # FIFTH: join prediction
