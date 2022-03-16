import errno
import os
import sys
import shlex
import subprocess


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


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


def print_error():
    try:
        sys.stderr.write(f"\n{bco.Red}{bco.Bold}[E::main] Error: {bco.ResetAll}")
    except Exception:
        sys.stderr.write("[E::main] Error: ")


# function that checks if a file exists ----------------------------------------
def check_file_exists(file_name, isfasta=False):
    try:
        with open(file_name, "r") as _in:
            # if fasta file, then check that it starts with ">"
            if isfasta:
                if _in.readline()[0] != ">":
                    sys.stderr.write(f"{bco.Red}{bco.Bold}[E::main] Error: {bco.ResetAll}")
                    sys.stderr.write(f"Not a fasta file: {file_name}\n")
                    sys.stderr.write("          Fasta file is expected to start with '>'\n")
                    sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"{bco.Red}{bco.Bold}[E::main] Error: {bco.ResetAll}")
        sys.stderr.write(f"Cannot open file: {file_name}\n")
        sys.stderr.write(f"{e}\n")
        sys.exit(1)


# function that checks if a file exists already, and give an error -------------
def check_file_doesnt_exists(file_name):
    if os.path.exists(file_name):
        print_error()
        sys.stderr.write("Output file exists already: "+file_name+"\n")
        sys.exit(1)


def is_tool(name):
    with open(os.devnull) as devnull:
        try:
            with subprocess.Popen([name], stdout=devnull, stderr=devnull) as proc:
                proc.communicate()
        except OSError as e:
            return e.errno != errno.ENOENT

    return True


def is_tool_and_return0(name):
    with open(os.devnull) as devnull:
        try:
            with open(os.devnull) as devnull:
                with subprocess.Popen(shlex.split(name), stdout=devnull, stderr=devnull) as child:
                    _ = child.communicate()
                    rc = child.wait()
                    return rc == 0
        except OSError as e:
            return e.errno != errno.ENOENT

    return True


def read_fasta(fasta, head_start=0, is_binary=True):
    sid, seq = None, []
    fasta_stream = open(fasta) if isinstance(fasta, str) else fasta
    with fasta_stream:
        for line in fasta_stream:
            if is_binary:
                line = line.decode("utf-8")
            line = line.rstrip()
            if line.startswith(">"):
                if seq:
                    yield sid, "".join(seq)
                sid = line[head_start:]
                seq.clear()
            else:
                seq.append(line)
        if seq:
            yield sid, "".join(seq)
