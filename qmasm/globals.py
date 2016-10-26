###################################
# QMASM global variables          #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import qmasm
import sys

# Name of this program
qmasm.progname = sys.argv[0]

# Map from a symbol to a unique number
qmasm.sym2num = {}

# One less than the next symbol number to assign
qmasm.next_sym_num = -1

# List of Statement objects
qmasm.program = []

# Define our internal representation.
qmasm.chain_strength = 0    # Strength of chain couplers
qmasm.pin_strength = 0      # Strength of pin couplers
