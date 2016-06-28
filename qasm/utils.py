###################################
# QASM utility functions          #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import qasm

# Define a function that maps from a symbol to a number, creating a
# new association if necessary.
def symbol_to_number(sym):
    try:
        return qasm.sym2num[sym]
    except KeyError:
        qasm.next_sym_num += 1
        qasm.sym2num[sym] = qasm.next_sym_num
        return qasm.next_sym_num

# Define a function to abort the program on an error.
def abend(str):
    sys.stderr.write("%s: %s\n" % (progname, str))
    sys.exit(1)
