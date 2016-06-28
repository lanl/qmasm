###################################
# QASM global variables           #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import sys

# Name of this program
progname = sys.argv[0]

# Map from a symbol to a unique number
sym2num = {}

# One less than the next symbol number to assign
next_sym_num = -1
