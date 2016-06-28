###################################
# QASM utility functions          #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from qasm.globals import *

# Define a function that maps from a symbol to a number, creating a
# new association if necessary.
def symbol_to_number(sym):
    global sym2num, next_sym_num
    try:
        return sym2num[sym]
    except KeyError:
        next_sym_num += 1
        sym2num[sym] = next_sym_num
        return next_sym_num

# Define a function to abort the program on an error.
def abend(str):
    global progname
    sys.stderr.write("%s: %s\n" % (progname, str))
    sys.exit(1)
