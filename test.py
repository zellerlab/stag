#!/usr/bin/env python

# ============================================================================ #
# test.py: test stag
#
# Author: Alessio Milanese (milanese@embl.de)
#
# ============================================================================ #

import os
import sys
import tempfile
import subprocess
import shlex
import errno

# position of the script -------------------------------------------------------
path_mOTUs = os.path.realpath(__file__)
path_array = path_mOTUs.split("/")
relative_path = "/".join(path_array[0:-1])
relative_path = relative_path + "/"


# ------------------------------------------------------------------------------
# function to check if a specific tool exists
# ------------------------------------------------------------------------------
def is_tool_and_return0(name):
    try:
        devnull = open(os.devnull)
        popenCMD = shlex.split(name)
        child = subprocess.Popen(popenCMD, stdout=devnull, stderr=devnull)
        streamdata = child.communicate()
        rc = child.wait()
        if rc == 0:
            return True
        else:
            return False
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

def is_tool(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True


# colors for the shell ---------------------------------------------------------
class bco:
    ResetAll = "\033[0m"
    Bold       = "\033[1m"
    Underlined = "\033[4m"
    Green        = "\033[32m"
    Yellow       = "\033[33m"
    Blue         = "\033[34m"
    Red          = "\033[31m"
    Magenta      = "\033[35m"
    Cyan         = "\033[36m"
    LightRed     = "\033[91m"
    LightGreen   = "\033[92m"
    LightYellow  = "\033[93m"
    LightBlue    = "\033[94m"
    LightMagenta = "\033[95m"
    LightCyan    = "\033[96m"

# ------------------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------------------
def main(argv=None):
    sys.stderr.write(f"{bco.Blue}{bco.Bold} ------------------------------------------------------------------------------{bco.ResetAll}\n")
    sys.stderr.write(f"{bco.Blue}{bco.Bold}|{bco.Green}                                    TEST STAG                                 {bco.Blue}|{bco.ResetAll}\n")
    sys.stderr.write(f"{bco.Blue}{bco.Bold} ------------------------------------------------------------------------------{bco.ResetAll}\n")

    error_found = False

    # CHECK TOOLS ==============================================================

    sys.stderr.write(f"{bco.Cyan}{bco.Bold}1-- Tools and versions:{bco.ResetAll}\n")
    # check python version -----------------------------------------------------
    sys.stderr.write("  ■ python:       ")
    python_version = sys.version_info
    if(python_version >= (3,0,0)):
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} WARNING: python2 is not supported{bco.ResetAll}\n\n")
        error_found = True

    # check hmmer --------------------------------------------------------------
    sys.stderr.write("  ■ hmmalign:     ")
    if is_tool("hmmalign"):
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} WARNING. HMMER3 is not in the path{bco.ResetAll}\n\n")
        error_found = True

    # check Easel --------------------------------------------------------------
    sys.stderr.write("  ■ esl-reformat: ")
    if is_tool("esl-reformat"):
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} WARNING. EASEL is not in the path{bco.ResetAll}\n\n")
        error_found = True

    # check Easel --------------------------------------------------------------
    sys.stderr.write("  ■ seqtk:        ")
    if is_tool("seqtk"):
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} WARNING. seqtk is not in the path{bco.ResetAll}\n\n")
        error_found = True

    # Python libraries:
    sys.stderr.write("  ■ (L)numpy:     ") #------------------------------------
    library_correct = True
    try:
        import numpy
    except ImportError as e:
        library_correct = False
    if library_correct:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} WARNING. numpy is missing{bco.ResetAll}\n\n")
        error_found = True

    sys.stderr.write("  ■ (L)pandas:    ") #------------------------------------
    library_correct = True
    try:
        import pandas
    except ImportError as e:
        library_correct = False
    if library_correct:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} WARNING. pandas is missing{bco.ResetAll}\n\n")
        error_found = True

    sys.stderr.write("  ■ (L)sklearn:   ") #------------------------------------
    library_correct = True
    try:
        import sklearn
    except ImportError as e:
        library_correct = False
    if library_correct:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} WARNING. sklearn is missing{bco.ResetAll}\n\n")
        error_found = True

    sys.stderr.write("  ■ (L)h5py:      ") #------------------------------------
    library_correct = True
    try:
        import h5py
    except ImportError as e:
        library_correct = False
    if library_correct:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} WARNING. h5py is missing{bco.ResetAll}\n\n")
        error_found = True


    if (error_found):
        return 1
    else:
        return 0        # success


#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
    status = main()
    sys.exit(status)
