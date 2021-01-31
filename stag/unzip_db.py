"""
Scripts that saves the database to a directory
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

# Requirements:
# - numpy
# - h5py

import numpy as np
import sys
import time
import os
import h5py
import tempfile
import shutil

from stag.load_db import load_db

#===============================================================================
#                                      MAIN
#===============================================================================
def unzip_db(database, verbose, dir_output):
    # check if the output dir exists already
    if os.path.isdir(dir_output):
        sys.stderr.write("Error, output dir exists already\n")
        sys.exit(1)
    else:
        try:
            os.mkdir(dir_output)
        except OSError:
            sys.stderr.write("Error, failed to create the output directory\n")
            sys.exit(1)

    # load the database and save it
    load_db(database, dir_output=dir_output)
