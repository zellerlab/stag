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
