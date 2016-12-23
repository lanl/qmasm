#! /usr/bin/env python3

###################################
# Convert Qubist to QMASM         #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import argparse
import sys

# Parse the command line.
cl_parser = argparse.ArgumentParser(description="Convert Qubist input to QMASM input")
cl_parser.add_argument("input", nargs="?", metavar="FILE", default="-",
                       help="Qubist-format input file (default: standard input)")
cl_parser.add_argument("-o", "--output", metavar="FILE", default="-",
                           help="file to which to write QMASM code (default: stdandard output)")
cl_parser.add_argument("-f", "--format", metavar="FORMAT", default="%d",
                           help='printf-style format string for formatting qubit numbers (default: "%%d")')
cl_parser.add_argument("-r", "--renumber-from", metavar="INT", type=int,
                           help="starting number from which to renumber qubits")
cl_args = cl_parser.parse_args()

# Open the input file.
if cl_args.input == "-":
    infile = sys.stdin
else:
    try:
        infile = open(cl_args.input, "r")
    except IOError:
        sys.stderr.write("%s: Failed to open %s for input\n" % (sys.argv[0], cl_args.input))
        sys.exit(1)

# Open the output file.
if cl_args.output == "-":
    outfile = sys.stdout
else:
    try:
        outfile = open(cl_args.output, "w")
    except IOError:
        sys.stderr.write("%s: Failed to open %s for output\n" % cl_args.output)
        sys.exit(1)

# Read the input file into memory, keeping track of all qubit numbers seen.
qubist = []
qnums = set()
for line in infile:
    fields = line.split()
    if len(fields) != 3:
        continue
    q1, q2, val = int(fields[0]), int(fields[1]), fields[2]
    qnums.add(q1)
    qnums.add(q2)
    qubist.append((q1, q2, val))
qnums = sorted(qnums)

# Map old qubit numbers to new qubit numbers.
if cl_args.renumber_from != None:
    newq = dict(zip(qnums, range(cl_args.renumber_from, cl_args.renumber_from + len(qnums))))
else:
    newq = {q: q for q in qnums}

# Convert each line in turn.
for q1, q2, val in qubist:
    if q1 == q2:
        # Point weight
        fmt = cl_args.format + " %s\n"
        outfile.write(fmt % (newq[q1], val))
    else:
        # Coupler strength
        fmt = cl_args.format + " " + cl_args.format + " %s\n"
        outfile.write(fmt % (newq[q1], newq[q2], val))

# Wrap up.
if cl_args.input != "-":
    infile.close()
if cl_args.output != "-":
    outfile.close()
