###################################
# QASM utility functions          #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import qasm
import sys

# Define a function that maps from a symbol to a number, creating a
# new association if necessary.
def symbol_to_number(sym):
    global sym2num, next_sym_num
    try:
        return qasm.sym2num[sym]
    except KeyError:
        qasm.next_sym_num += 1
        qasm.sym2num[sym] = qasm.next_sym_num
        return qasm.next_sym_num

# Define a function to abort the program on an error.
def abend(str):
    sys.stderr.write("%s: %s\n" % (qasm.progname, str))
    sys.exit(1)

# Define a function that returns the topology of the chimera graph associated
# with a given solver.
def chimera_topology(solver):
    try:
        nominal_qubits = solver.properties["num_qubits"]
    except KeyError:
        # The Ising heuristic solver is an example of a solver that lacks a
        # fixed hardware representation.  We therefore claim an arbitrary
        # topology (the D-Wave 2X's 12x12x4x2 topology).
        return 4, 12, 12
    couplers = solver.properties["couplers"]
    deltas = [abs(c1 - c2) for c1, c2 in couplers]
    delta_tallies = {d: 0 for d in deltas}
    for d in deltas:
        delta_tallies[d] += 1
        sorted_tallies = sorted(delta_tallies.items(), key=lambda dt: dt[1], reverse=True)
    L = sorted_tallies[0][0]
    M = 1
    for d, t in sorted_tallies[1:]:
        if d > 2*L:
            M = d // (2*L)
            break
    N = (nominal_qubits + 2*L*M - 1) // (2*L*M)
    return L, M, N
