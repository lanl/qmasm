###################################
# QASM utility functions          #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from qasm.globals import *

# Define a function that maps from a symbol to a number, creating a
# new association if necessary.
def symbol_to_number(sym):
    try:
        return sym2num[sym]
    except KeyError:
        global next_sym_num
        next_sym_num += 1
        sym2num[sym] = next_sym_num
        return next_sym_num

# Define a function to abort the program on an error.
def abend(str):
    sys.stderr.write("%s: %s\n" % (progname, str))
    sys.exit(1)
