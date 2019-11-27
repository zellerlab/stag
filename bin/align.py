"""
Scripts to align sequences and transoform them into 1-hot encoding
"""

# Author: Alessio Milanese <milanese.alessio@gmail.com>

import shutil
import time
import subprocess
import shlex
import os
import errno
import sys

#===============================================================================
#                                 FUNCTIONS
#===============================================================================

# ------------------------------------------------------------------------------
# function to check if a specific tool exists
def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True
