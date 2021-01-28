import os
import sys
import shlex
import subprocess

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
	except Exception as e:
		sys.stderr.write("[E::main] Error: ")

## function that checks if a file exists ----------------------------------------
#def check_file_exists(file_name, isfasta = False):
#	try:
#		o = open(file_name,"r")
#		# if fasta file, then check that it starts with ">"
#		if isfasta:
#			if not(o.readline().startswith(">")):
#				print_error()
#				sys.stderr.write("Not a fasta file: "+file_name+"\n")
#				sys.stderr.write("          Fasta file is expected to start with '>'\n")
#				o.close()
#				sys.exit(1)
#		o.close()
#	except Exception as e:
#		print_error()
#		sys.stderr.write("Cannot open file: "+file_name+"\n")
#		sys.stderr.write(str(e)+"\n")
#		sys.exit(1)

# function that checks if a file exists ----------------------------------------
def check_file_exists(file_name, isfasta = False):
    try:
        o = open(file_name,"r")
        # if fasta file, then check that it starts with ">"
        if isfasta:
            if not(o.readline().startswith(">")):
                sys.stderr.write(f"{bco.Red}{bco.Bold}[E::main] Error: {bco.ResetAll}")
                sys.stderr.write("Not a fasta file: "+file_name+"\n")
                sys.stderr.write("          Fasta file is expected to start with '>'\n")
                o.close()
                sys.exit(1)
        o.close()
    except Exception as e:
        sys.stderr.write(f"{bco.Red}{bco.Bold}[E::main] Error: {bco.ResetAll}")
        sys.stderr.write("Cannot open file: "+file_name+"\n")
        sys.stderr.write(str(e)+"\n")
        sys.exit(1)






# function that checks if a file exists already, and give an error -------------
def check_file_doesnt_exists(file_name):
	if os.path.exists(file_name):
		print_error()
		sys.stderr.write("Output file exists already: "+file_name+"\n")
		sys.exit(1)

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

# ------------------------------------------------------------------------------
# function to convert a fasta file with multiple lines into a one line separated
# by a "\t"
# Example:
# >test_fasta_header
# ATTGCGATTTCT
# CGGTATCGGTAT
# CGGTTA
### TO:
# >test_fasta_header\tATTGCGATTTCTCGGTATCGGTATCGGTTA
def linearise_fasta(filein, head_start=0):
    # filein is a stream of data (from hmmalign)
    seq = ""
    for line_b in filein:
        line = line_b.decode("utf-8").rstrip()
        if line.startswith(">"):
            if seq:
                yield seq[head_start:] # we skip the 0 character, which is ">"
            seq = line + "\t"
        else:
            seq += line
    # give back the last sequence
    if seq:
        yield seq[head_start:] # we skip the 0 character, which is ">"

