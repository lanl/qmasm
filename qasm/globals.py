###################################
# QASM global variables           #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import qasm
import sys

# Name of this program
qasm.progname = sys.argv[0]

# Map from a symbol to a unique number
qasm.sym2num = {}

# One less than the next symbol number to assign
qasm.next_sym_num = -1
