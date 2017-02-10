#! /usr/bin/env python

###################################
# Wrap qbsolv when processing     #
# a QMASM-generated input file    #
#                                 #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import os
import subprocess
import sys

# Extract the name of the input file from the command line.
infile = None
nargs = len(sys.argv)
for i in range(1, nargs):
    if sys.argv[i] == "-i" and i < nargs - 1:
        infile = sys.argv[i + 1]
    if sys.argv[i] == "-o":
        sys.stderr.write("%s: qbsolv's -o option is not supported.  Aborting.\n" % sys.argv[0])
        sys.exit(1)
if infile == None:
    # No input file: Let qbsolv issue the error message.
    try:
        proc = subprocess.Popen(["qbsolv"] + sys.argv[1:])
    except OSError, e:
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
        if len(fields) == 4 and fields[0] == "c" and fields[2] == "-->":
            name2qubit[fields[1]] = int(fields[3])
if len(name2qubit) == 0:
    sys.stderr.write("%s: No QMASM comments found in %s.  Aborting.\n" % (sys.argv[0], infile))
    sys.exit(1)

# Run qbsolv and store its output bits.
try:
    proc = subprocess.Popen(["qbsolv"] + sys.argv[1:], stdout=subprocess.PIPE)
except OSError, e:
    sys.stderr.write("qbsolv: %s\n" % str(e))
    sys.exit(1)
pout, perr = proc.communicate()
retcode = proc.wait()
if retcode < 0:
    os.kill(os.getpid(), -retcode)
elif retcode > 0:
    sys.exit(retcode)
if perr != None:
    sys.stderr.write(perr)
output = pout.split("\n")
output = output[:-1]   # Skip empty final line.
bits = []
next_is_bits = False
energy = "?"
for line in output:
    sys.stderr.write("# %s\n" % line)
    if "Number of bits in solution" in line:
        next_is_bits = True
    elif next_is_bits:
        bits = [int(b) for b in list(line)]
        next_is_bits = False
    elif "Energy of solution" in line:
        energy = float(line.split()[0])
sys.stderr.write("\n")

# Report the output bits symbolically in QMASM style.
max_name_width = max([len(nm) for nm in name2qubit.keys()])
max_name_width = max(max_name_width, 7)
print("Solution #1 (energy = %.2f, tally = 1):\n" % energy)
print("    %-*s  Spin  Boolean" % (max_name_width, "Name(s)"))
print("    %s  ----  -------" % ("-" * max_name_width))
for name in sorted(name2qubit.keys()):
    num = name2qubit[name]
    spin = bits[num]*2 - 1
    print("    %-*s  %+4d  %s" % (max_name_width, name, spin, str(spin == +1)))
print("")
