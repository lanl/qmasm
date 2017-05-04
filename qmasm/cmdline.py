###################################
# Parse the QMASM command line    #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import argparse
import qmasm
import sys

def parse_command_line():
    "Parse the QMASM command line.  Return an argparse.Namespace."
    cl_parser = argparse.ArgumentParser(description="Assemble a symbolic Hamiltonian into a numeric one")
    cl_parser.add_argument("input", nargs="*",
                           help="file from which to read a symbolic Hamiltonian")
    cl_parser.add_argument("-v", "--verbose", action="count", default=0,
                           help="increase output verbosity (can be specified repeatedly)")
    cl_parser.add_argument("-r", "--run", action="store_true",
                           help="run the program on the current solver")
    cl_parser.add_argument("-o", "--output", metavar="FILE", default="<stdout>",
                           help="file to which to write weights and strengths")
    cl_parser.add_argument("-f", "--format", choices=["qubist", "dw", "qbsolv", "qmasm", "minizinc"], default="qubist",
                           help="output-file format")
    cl_parser.add_argument("-O", action="store_true",
                           help="optimize the layout (i.e., use fewer unit cells)")
    cl_parser.add_argument("-p", "--pin", action="append",
                           help="pin a set of qubits to a set of true or false values")
    cl_parser.add_argument("-d", "--discard", choices=["yes", "no", "maybe"], default="yes",
                           help="always, never, or if otherwise no solutions, discard solutions with broken chains or broken pins (default: yes)")
    cl_parser.add_argument("--values", choices=["bools", "ints"], default="bools",
                           help="output solution values as Booleans or integers (default: bools)")
    cl_parser.add_argument("-a", "--all-solns", action="store_true",
                           help='output all solutions, not just those at the minimal energy level (implied by "-v -v"')
    cl_parser.add_argument("-C", "--chain-strength", metavar="NEG_NUM", type=float,
                           help="negative-valued chain strength (default: automatic)")
    cl_parser.add_argument("-P", "--pin-strength", metavar="NEG_NUM", type=float,
                           help="negative-valued pin strength (default: automatic)")
    cl_parser.add_argument("-q", "--qubo", action="store_true",
                           help="treat inputs as QUBOs rather than Ising systems")
    cl_parser.add_argument("-s", "--samples", metavar="POS_INT", type=int, default=1000,
                           help="specify the number of samples to take (default: 1000)")
    cl_parser.add_argument("--anneal-time", metavar="POS_INT", type=int, default=20,
                           help="specify the annealing time in microseconds (default: 20)")
    cl_parser.add_argument("--spin-revs", metavar="POS_INT", type=int, default=0,
                           help="specify the number of spin-reversal transforms to perform (default: 0)")
    cl_parser.add_argument("--extra-args", default="",
                           help="specify extra arguments to pass to a solver command (default: none)")
    cl_parser.add_argument("--postproc", choices=["", "sample", "opt"],
                           default="",
                           help="specify the type of postprocessing to perform (default: none)")
    cl_args = cl_parser.parse_args()
    if cl_args.chain_strength != None and cl_args.chain_strength >= 0.0:
        sys.stderr.write("%s: Warning: A non-negative chain strength (%.20g) was specified\n" % (qmasm.progname, cl_args.chain_strength))
    if cl_args.pin_strength != None and cl_args.pin_strength >= 0.0:
        sys.stderr.write("%s: Warning: A non-negative pin strength (%.20g) was specified\n" % (qmasm.progname, cl_args.pin_strength))
    return cl_args
