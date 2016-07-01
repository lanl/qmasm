###################################
# QASM global variables           #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from collections import defaultdict
import qasm
import sys

# Name of this program
qasm.progname = sys.argv[0]

# Map from a symbol to a unique number
qasm.sym2num = {}

# One less than the next symbol number to assign
qasm.next_sym_num = -1

# List of Statement objects
qasm.program = []

# Define our internal representation.
qasm.weights = defaultdict(lambda: 0.0)    # Map from a spin to a point weight
qasm.strengths = defaultdict(lambda: 0.0)  # Map from a pair of spins to a coupler strength
qasm.chains = {}     # Subset of strengths keys that represents chains
qasm.pinned = []     # Pairs of {unique number, Boolean} to pin
qasm.chain_strength = 0    # Strength of chain couplers
qasm.pin_strength = 0      # Strength of pin couplers
