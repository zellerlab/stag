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


#===============================================================================
#                                      MAIN
#===============================================================================
def classify_genome(database, fasta_input, verbose, threads, output, long_out, tool_version):
    # FIRST: unzip the database
    # SECOND: run prodigal on the fasta genome
    # THIRD: find the marker genes from the predicted genes
    # FOURTH: classify the marker genes
    # FIFTH: join prediction
