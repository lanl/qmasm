#!/usr/bin/env python

###################################
# Wrap qbsolv when processing     #
# a QMASM-generated input file    #
#                                 #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import os
import qmasm
import re
import subprocess
import sys

# Parse the command line.
infile = None
nargs = len(sys.argv)
style = "bools"
verbosity = 0
pin_weight = -1
chain_strength = -1
qbsolv_args = []   # Subset of sys.argv to pass to qbsolv
for i in range(1, nargs):
    arg = sys.argv[i]
    qbsolv_args.append(arg)
    if arg == "-i" and i < nargs - 1:
        infile = sys.argv[i + 1]
    elif arg == "-o":
        sys.stderr.write("%s: qbsolv's -o option is not supported.  Aborting.\n" % sys.argv[0])
        sys.exit(1)
    elif arg[:9] == "--values=":
        # Accept qmasm's --values argument.
        style = arg[9:]
        qbsolv_args.pop()
    elif arg[:12] == "--verbosity=":
        # Accept a variation of qmasm's --verbose argument.
        verbosity = int(arg[12:])
        qbsolv_args.pop()
    elif arg[:13] == "--pin-weight=":
        # Accept qmasm's --pin-weight argument.
        pin_weight = float(arg[13:])
        qbsolv_args.pop()
    elif arg[:17] == "--chain-strength=":
        # Accept qmasm's --chain-strength argument.
        chain_strength = float(arg[17:])
        qbsolv_args.pop()
if infile == None:
    # No input file: Let qbsolv issue the error message.
    try:
        proc = subprocess.Popen(["qbsolv"] + qbsolv_args, stderr=sys.stderr, stdout=sys.stdout)
    except OSError as e:
        sys.stderr.write("qbsolv: %s\n" % str(e))
        sys.exit(1)
    retcode = proc.wait()
    if retcode < 0:
        os.kill(os.getpid(), -retcode)
    else:
        sys.exit(retcode)

# Parse the comments inserted by QMASM in the input file.
name2qubit = {}     # Map from a QMASM name to a qubit number
with open(infile) as f:
    for line in f:
        fields = line.split()
        if len(fields) >= 4 and fields[0] == "c" and fields[2] == "-->":
            name2qubit[fields[1]] = int(fields[3])
if len(name2qubit) == 0:
    sys.stderr.write("%s: No QMASM comments found in %s.  Aborting.\n" % (sys.argv[0], infile))
    sys.exit(1)

# Run qbsolv and store the solution bits and solution energy it outputs.
try:
    proc = subprocess.Popen(["qbsolv"] + qbsolv_args, stdout=subprocess.PIPE, stderr=sys.stderr)
except OSError as e:
    sys.stderr.write("qbsolv: %s\n" % str(e))
    sys.exit(1)
bits = []
bits_re = re.compile(r'^[01]+$')
energy = "?"
while True:
    line = proc.stdout.readline()
    if line == "":
        break
    line = line.rstrip()
    sys.stderr.write("# %s\n" % line)
    if bits_re.match(line):
        bits = [int(b) for b in list(line)]
    elif b"Energy of solution" in line:
        energy = float(line.split()[0])
retcode = proc.wait()
if retcode < 0:
    os.kill(os.getpid(), -retcode)
elif retcode > 0:
    # Some qbsolv errors go to stdout, not stderr.
    for line in proc.stdout:
        sys.stderr.write("qbsolv: %s\n" % line)
    sys.exit(retcode)
sys.stderr.write("\n")

# Fake various QMASM objects.
qmasm.pin_weight = pin_weight
qmasm.chain_strength = chain_strength
for n, q in sorted(name2qubit.items(), key=lambda k: k[1]):
    qmasm.sym_map.new_symbol(n)
answer = {"num_occurrences": [1],
          "energies": [energy],
          "solutions": [[2*b - 1 for b in bits]]}
problem = qmasm.Problem(False)
problem.embedding = [[i] for i in range(len(bits))]
problem.embedder_chains = set()
solutions = qmasm.Solutions(answer, problem, verbosity >= 2)

# Output the solution.  For now, we hard-wire show_asserts to False.
qmasm.output_solution(solutions, style, verbosity, False)
