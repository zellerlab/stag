#!/usr/bin/env python

# ============================================================================ #
# stag_test Run stag tests
#
# Author: Alessio Milanese (milanese@embl.de)
#
# ============================================================================ #

import time
import os
import sys
import tempfile
import subprocess
import shlex
import errno
import pkg_resources
import urllib.request
import hashlib
from pathlib import Path
import shutil

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


# check md5
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

# download a file if it's not there
def download_file(url, filename):
    downloaded_correct = True
    try:
        urllib.request.urlretrieve(url, filename)
    except:
        downloaded_correct = False
    if downloaded_correct:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} error{bco.ResetAll}\n")


# function to download and decompress a tar.gz
def download_and_checkmd5_and_decompress(url, filename, md5_db, destination):
    # we remove a dir if it exist already
    shutil.rmtree(filename[0:-7], ignore_errors=True)
    # check if the file is already downloaded
    my_file = Path(filename)
    if my_file.is_file():
        sys.stderr.write("  ■ already downloaded\n")
    else:
        sys.stderr.write("  ■ download file:     ")
        download_file(url, filename)

    # check md5
    sys.stderr.write("  ■ check md5:            ")
    current_md5 = md5(filename)
    if current_md5 == md5_db:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} different{bco.ResetAll}\n")
        sys.stderr.write("  ■ Re-download file:     ")
        download_file(url, filename)
        sys.stderr.write("  ■ check md5:            ")
        current_md5 = md5(filename)
        if current_md5 == md5_db:
            sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")
        else:
            sys.stderr.write(f"{bco.Yellow}{bco.Bold} different{bco.ResetAll}\n")

    # decompress
    sys.stderr.write("  ■ Unzip file:           ")
    extract_cmd = "tar -zxvf "+filename+" -C "+destination
    try:
        FNULL = open(os.devnull, 'w')
        process = subprocess.Popen(extract_cmd.split(),stderr=FNULL,stdout=FNULL)
        output, error = process.communicate()
    except:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} error{bco.ResetAll}\n")
    if process.returncode:
        sys.stderr.write(f"{bco.Yellow}{bco.Bold} error{bco.ResetAll}\n")
    else:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")

    # we return the reulting directory, which should have the same name as
    # the file downloaded minus ".tar.gz"
    return filename[0:-7]+"/"



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

    t0 = time.time()
    stag_command = "stag train -f -o "+temp_file_db.name+" -i "+seq_file+" -x "+tax_file+" -a "+hmm_file
    process = subprocess.run(stag_command.split())
    runtime = time.time() - t0

    if process.returncode:
        sys.stderr.write(f"{bco.Red}{bco.Bold} Error{bco.ResetAll} ({runtime:.3f}s)\n")
        sys.exit(1)
    else:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll} ({runtime:.3f}s)\n")

    sys.stderr.write("  ■ classify:   ") #--------------------------------------
    sys.stderr.flush()
    temp_file_res = tempfile.NamedTemporaryFile(delete=False, mode="w")

    t0 = time.time()
    stag_command = "stag classify -v 1 -d "+temp_file_db.name+" -i "+seq_file+" -o "+temp_file_res.name
    process = subprocess.run(stag_command.split())
    runtime = time.time() - t0

    if process.returncode:
        sys.stderr.write(f"{bco.Red}{bco.Bold} Error{bco.ResetAll} ({runtime:.3f}s)\n")
        sys.exit(1)
    else:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll} ({runtime:.3f}s)\n")

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



    ############################################################################
    sys.stderr.write(f"{bco.Blue}{bco.Bold} ------------------------------------------------------------------------------{bco.ResetAll}\n")
    sys.stderr.write(f"{bco.Blue}{bco.Bold}|{bco.Green}                                    LONG STAG TEST                             {bco.Blue}|{bco.ResetAll}\n")
    sys.stderr.write(f"{bco.Blue}{bco.Bold} ------------------------------------------------------------------------------{bco.ResetAll}\n")

    # long test part 1: test building genome DB --------------------------------
    sys.stderr.write(f"{bco.Cyan}{bco.Bold}1-- Build genome database:{bco.ResetAll}\n")
    # prepare data
    link_db = "https://zenodo.org/record/4626959/files/train_genome_files.tar.gz"
    db_name = os.path.join(TEST_DATA_PATH, "train_genome_files.tar.gz")
    md5_db = "5ec5a527d25cc6d1f11a8ec50cd252a7"
    this_dir = download_and_checkmd5_and_decompress(link_db, db_name, md5_db, TEST_DATA_PATH)
    # train genome
    sys.stderr.write("  ■ train genome DB:      ")
    sys.stderr.flush()
    gene_files = this_dir + "COG0012," + this_dir + "COG0016," + this_dir + "COG0018"
    merged_db = this_dir + "genes_ali.stagDB"
    thresholds = this_dir + "gene_thresholds"
    result_genome_DB = this_dir + "TEST_DB.stagDB"

    t0 = time.time()
    stag_command = "stag train_genome -v 1 -o "+result_genome_DB+" -i "+gene_files+" -T "+thresholds+" -C "+merged_db
    process = subprocess.run(stag_command.split())
    runtime = time.time() - t0

    if process.returncode:
        sys.stderr.write(f"{bco.Red}{bco.Bold} Error{bco.ResetAll} ({runtime:.3f}s)\n")
        sys.exit(1)
    else:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll} ({runtime:.3f}s)\n")





    # long test part 2: test classify genome -----------------------------------
    sys.stderr.write(f"{bco.Cyan}{bco.Bold}2-- Test genome classification:{bco.ResetAll}\n")
    # prepare data
    link_db = "https://zenodo.org/record/4626959/files/classify_genome_files.tar.gz"
    db_name = os.path.join(TEST_DATA_PATH, "classify_genome_files.tar.gz")
    md5_db = "819cc77d463a797a330d8d1d9437feca"
    this_dir = download_and_checkmd5_and_decompress(link_db, db_name, md5_db, TEST_DATA_PATH)
    # train genome
    sys.stderr.write("  ■ classify 2 genomes:   ")
    sys.stderr.flush()
    genome_files = this_dir + "genomes"
    result = this_dir + "RESULT_TEMP"

    t0 = time.time()
    stag_command = "stag classify_genome -v 1 -o "+result+" -d "+result_genome_DB+" -D "+genome_files
    process = subprocess.run(stag_command.split())
    runtime = time.time() - t0

    if process.returncode:
        sys.stderr.write(f"{bco.Red}{bco.Bold} Error{bco.ResetAll} ({runtime:.3f}s)\n")
        sys.exit(1)
    else:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll} ({runtime:.3f}s)\n")

    # check result of the classification ---------------------------------------
    sys.stderr.write("  ■ check result:         ")
    # correct annotation
    o = open(this_dir + "CORRECT_ASSIGNMENT/genome_annotation")
    o.readline()
    correct_classification = dict()
    all_genomes = dict()
    for line in o:
        vals = line.rstrip().split("\t")
        correct_classification[vals[0]] = vals[1]
        all_genomes[vals[0]] = False
    o.close()
    # this annotation
    o = open(this_dir + "RESULT_TEMP/genome_annotation")
    o.readline()
    for line in o:
        vals = line.rstrip().split("\t")
        vals[0] = vals[0].split("/")[-1]
        all_genomes[vals[0]] = True
        if not vals[0] in correct_classification:
            sys.stderr.write(f"{bco.Red}{bco.Bold} Error, too many lines{bco.ResetAll}\n")
            sys.stderr.write(f"{bco.Red}{bco.Bold} Check "+ this_dir + f"RESULT_TEMP/genome_annotation {bco.ResetAll}\n")
            sys.exit(1)
        else:
            if correct_classification[vals[0]] != vals[1]:
                sys.stderr.write(f"{bco.Red}{bco.Bold} Error, wrong calssification{bco.ResetAll}\n")
                sys.stderr.write(f"{bco.Red}{bco.Bold} Check "+ this_dir + f"RESULT_TEMP/genome_annotation {bco.ResetAll}\n")
                sys.exit(1)
    o.close()
    # check that all the genomes were profiled
    for genome in all_genomes:
        if not all_genomes[genome]:
            sys.stderr.write(f"{bco.Red}{bco.Bold} Error, some genomes are missing{bco.ResetAll}\n")
            sys.stderr.write(f"{bco.Red}{bco.Bold} Check "+ this_dir + f"RESULT_TEMP/genome_annotation {bco.ResetAll}\n")
            sys.exit(1)
    # if we arrive till here, then it's correct
    sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")






    # long test part 3: test gene train and classification with real data ------
    sys.stderr.write(f"{bco.Cyan}{bco.Bold}3-- Test real genes:{bco.ResetAll}\n")
    # prepare data
    link_db = "https://zenodo.org/record/4626959/files/test_gene.tar.gz"
    db_name = os.path.join(TEST_DATA_PATH, "test_gene.tar.gz")
    md5_db = "bee91d9dc06fae153502af386c29ca5c"
    this_dir = download_and_checkmd5_and_decompress(link_db, db_name, md5_db, TEST_DATA_PATH)
    # train genome
    sys.stderr.write("  ■ train database:       ")
    sys.stderr.flush()
    seq_file = this_dir + "train.fna"
    protein_file = this_dir + "train.faa"
    tax_file = this_dir + "train.tax"
    hmm_file = this_dir + "COG0012.hmm"
    trained_db = this_dir + "TRAINED_TEMP"
    temp_file_db = tempfile.NamedTemporaryFile(delete=False, mode="w")

    t0 = time.time()
    stag_command = "stag train -f -o "+trained_db+" -i "+seq_file+" -p "+protein_file+" -x "+tax_file+" -a "+hmm_file
    process = subprocess.run(stag_command.split())
    runtime = time.time() - t0

    if process.returncode:
        sys.stderr.write(f"{bco.Red}{bco.Bold} Error{bco.ResetAll} ({runtime:.3f}s)\n")
        sys.exit(1)
    else:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll} ({runtime:.3f}s)\n")

    sys.stderr.write("  ■ classify:             ") #--------------------------------------
    sys.stderr.flush()
    res_classification = this_dir + "RES_TEMP"
    seq_file = this_dir + "test.fna"
    protein_file = this_dir + "test.faa"

    t0 = time.time()
    stag_command = "stag classify -v 1 -d "+trained_db+" -i "+seq_file+" -p "+protein_file+" -o "+res_classification
    process = subprocess.run(stag_command.split())
    runtime = time.time() - t0

    if process.returncode:
        sys.stderr.write(f"{bco.Red}{bco.Bold} Error{bco.ResetAll} ({runtime:.3f}s)\n")
        sys.exit(1)
    else:
        sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll} ({runtime:.3f}s)\n")


    # check result of the classification ---------------------------------------
    sys.stderr.write("  ■ check result:         ")
    # correct annotation
    o = open(this_dir + "test.CORRECT_ASSIGNMENT")
    o.readline()
    correct_classification = dict()
    all_genomes = dict()
    for line in o:
        vals = line.rstrip().split("\t")
        correct_classification[vals[0]] = vals[1]
        all_genomes[vals[0]] = False
    o.close()
    # this annotation
    o = open(res_classification)
    o.readline()
    for line in o:
        vals = line.rstrip().split("\t")
        all_genomes[vals[0]] = True
        if not vals[0] in correct_classification:
            sys.stderr.write(f"{bco.Red}{bco.Bold} Error, too many lines{bco.ResetAll}\n")
            sys.stderr.write(f"{bco.Red}{bco.Bold} Check "+ res_classification + f"{bco.ResetAll}\n")
            sys.exit(1)
        else:
            if correct_classification[vals[0]] != vals[1]:
                sys.stderr.write(f"\n{bco.Yellow} Corr: "+correct_classification[vals[0]]+f"{bco.ResetAll}\n")
                sys.stderr.write(f"{bco.Yellow} Pred: "+vals[1]+f"{bco.ResetAll}\n")
    o.close()
    # check that all the genomes were profiled
    for genome in all_genomes:
        if not all_genomes[genome]:
            sys.stderr.write(f"{bco.Red}{bco.Bold} Error, some genomes are missing{bco.ResetAll}\n")
            sys.stderr.write(f"{bco.Red}{bco.Bold} Check "+ res_classification + f"{bco.ResetAll}\n")
            sys.exit(1)
    # if we arrive till here, then it's correct
    sys.stderr.write(f"{bco.Green}{bco.Bold} correct{bco.ResetAll}\n")






    return None        # success


#-------------------------------- run main -------------------------------------
if __name__ == '__main__':
    main()
