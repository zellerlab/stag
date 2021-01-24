#!/usr/bin/env python

# ============================================================================ #
# stag_test Run stag tests
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
import pkg_resources

from .helpers import bco

TEST_DATA_PATH = pkg_resources.resource_filename("stag", "test")

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

    # TRY TO RUN STAG ==========================================================
    sys.stderr.write(f"{bco.Cyan}{bco.Bold}2-- Run stag:{bco.ResetAll}\n")

    sys.stderr.write("  ■ train:      ") #--------------------------------------
    sys.stderr.flush()
    seq_file = os.path.join(TEST_DATA_PATH, "sequences.fasta")
    tax_file = os.path.join(TEST_DATA_PATH, "sequences.taxonomy")
    hmm_file = os.path.join(TEST_DATA_PATH, "gene.hmm")
    temp_file_db = tempfile.NamedTemporaryFile(delete=False, mode="w")


    stag_command = "stag train -f -o "+temp_file_db.name+" -i "+seq_file+" -x "+tax_file+" -a "+hmm_file
    process = subprocess.run(stag_command.split())

    if process.returncode:
        sys.stderr.write(f"{bco.Red}{bco.Bold} Error{bco.ResetAll}\n")
        sys.exit(1)
    else:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")

    sys.stderr.write("  ■ classify:   ") #--------------------------------------
    sys.stderr.flush()
    temp_file_res = tempfile.NamedTemporaryFile(delete=False, mode="w")

    stag_command = "stag classify -v 1 -d "+temp_file_db.name+" -i "+seq_file+" -o "+temp_file_res.name
    process = subprocess.run(stag_command.split())

    if process.returncode:
        sys.stderr.write(f"{bco.Red}{bco.Bold} Error{bco.ResetAll}\n")
        sys.exit(1)
    else:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")

    # remove temp file
    os.remove(temp_file_db.name+".log")
    os.remove(temp_file_db.name)

    # CHECK THE RESULTING FILE =================================================
    sys.stderr.write(f"{bco.Cyan}{bco.Bold}3-- Check result of the classification:{bco.ResetAll}\n")
    sys.stderr.write("  ■ taxonomy of classified sequences: ")
    sys.stderr.flush()

    o = open(tax_file,"r")
    correct_tax = dict()
    for i in o:
        vals = i.rstrip().split("\t")
        correct_tax[vals[0]] = vals[1]
    o.close()

    o = open(temp_file_res.name,"r")
    o.readline() # remove header
    pred_tax = dict()
    for i in o:
        vals = i.rstrip().split("\t")
        if len(vals) < 2:
            sys.stderr.write(f"{bco.Red}{bco.Bold} Error: less than two values ("+str(vals)+f"){bco.ResetAll}\n")
            os.remove(temp_file_res.name)
            sys.exit(1)
        pred_tax[vals[0]] = vals[1]
    o.close()

    # remove temp file
    os.remove(temp_file_res.name)

    # let's check the values
    if not set(pred_tax.keys()) == set(correct_tax.keys()):
        sys.stderr.write(f"{bco.Red}{bco.Bold} Error: different number of predicted genes{bco.ResetAll}\n")
        sys.exit(1)
    # if we arrive here, we have the same set of predicted genes
    # let's check the predicted taxonomies
    error_flag = False
    for i in pred_tax:
        if pred_tax[i] != correct_tax[i]:
            error_flag = True
            sys.stderr.write(f"{bco.Red}{bco.Bold} Error: different taxonomy for "+i+"(correct:'"+correct_tax[i]+"', predicted:'"+pred_tax[i]+"')"+f"{bco.ResetAll}\n")

    if not error_flag:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")


    if error_found or error_flag:
        sys.exit(1)

    return None        # success


#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
    main()
