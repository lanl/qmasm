###################################
# QASM global variables           #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from collections import defaultdict
import sys

# Name of this program
progname = sys.argv[0]

# Map from a symbol to a unique number
sym2num = {}

# One less than the next symbol number to assign
next_sym_num = -1

# List of Statement objects
program = []

# Define our internal representation.
weights = defaultdict(lambda: 0.0)    # Map from a spin to a point weight
strengths = defaultdict(lambda: 0.0)  # Map from a pair of spins to a coupler strength
chains = {}     # Subset of strengths keys that represents chains
pinned = []     # Pairs of {unique number, Boolean} to pin

# Specify the minimum distinguishable difference between energy readings.
min_energy_delta = 0.005

# Strength of each chain and each pinned variable (to be assigned later)
chain_strength = -1.0
pin_strength = -1.0
