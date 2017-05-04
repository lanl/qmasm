##############################################
# Run various classical (non-D-Wave) solvers #
# By Scott Pakin <pakin@lanl.gov>            #
##############################################

import os
import qmasm
import re
import subprocess
import sys
import tempfile

def run_qbsolv(logical_ising, oname, extra_args, verbosity):
    "Run qmasm-qbsolv on the problem and report the result."
    # Use the specified file name if provided.  Otherwise, write to a temporary
    # file."
    if oname == "<stdout>":
        tfile = tempfile.NamedTemporaryFile(suffix=".qubo", prefix="qmasm-", mode="w+")
        qubo_fname = tfile.name
        tfile.close()
    else:
        qubo_fname = oname

    # Create the .qubo file.
    qmasm.write_output(logical_ising, qubo_fname, "qbsolv", False)

    # Run qmasm-qbsolv on the .qubo file.
    args = ["qmasm-qbsolv", "-i", qubo_fname] + extra_args
    if verbosity >= 1:
        sys.stderr.write("Submitting the problem to qbsolv via qmasm-qbsolv.\n\n")
        if verbosity >= 2:
            sys.stderr.write("    Command line: %s\n\n" % str(" ".join(args)))
    try:
        subprocess.call(args, stdout=sys.stdout, stderr=sys.stderr)
    except OSError as e:
        qmasm.abend("Failed to run qmasm-qbsolv (%s)" % str(e.args[1]))

    # Delete the .qubo file if it's considered temporary.
    if oname == "<stdout>":
        os.remove(qubo_fname)

def run_minizinc(logical_ising, oname, extra_args, verbosity):
    "Run mzn-chuffed on the problem and report the result."
    # If the user specified a file, create that first.
    if oname != "<stdout>":
        qmasm.write_output(logical_ising, oname, "minizinc", False)

    # Always create a temporary file because we need something we can overwrite.
    tfile = tempfile.NamedTemporaryFile(suffix=".mzn", prefix="qmasm-", mode="w+")
    mzn_fname = tfile.name
    tfile.close()
    qmasm.write_output(logical_ising, mzn_fname, "minizinc", False)

    # Run mzn-chuffed on the .mzn file.
    args = ["mzn-chuffed"] + extra_args + [mzn_fname]
    if verbosity >= 1:
        sys.stderr.write("Submitting the problem to MiniZinc (specifically, mzn-chuffed).\n\n")
        if verbosity >= 2:
            sys.stderr.write("    Command line: %s\n" % str(" ".join(args)))
    try:
        mzn_output = subprocess.check_output(args, stderr=sys.stderr, universal_newlines=True)
    except Exception as e:
        qmasm.abend("Failed to run mzn-chuffed (%s)" % str(e))

    # Extract the energy value.
    match = re.search(r'\benergy = (-?\d+),', mzn_output)
    if match == None:
        qmasm.abend("Failed to find an energy level in MiniZinc output")
    energy = int(match.group(1))

    # Regenerate the MiniZinc input as a satisfaction problem instead of a
    # minimization problem.
    tfile = open(mzn_fname, "w")
    qmasm.output_minizinc(tfile, logical_ising, energy)
    tfile.close()

    # Run mzn-chuffed on the .mzn file.
    args = ["mzn-chuffed"]
    args += ["--all-solutions", "--solution-separator= ", "--search-complete-msg="]
    args += extra_args + [mzn_fname]
    if verbosity >= 2:
        sys.stderr.write("    Command line: %s\n" % str(" ".join(args)))
    try:
        mzn_output = subprocess.check_output(args, stderr=sys.stderr, universal_newlines=True)
    except Exception as e:
        qmasm.abend("Failed to run mzn-chuffed (%s)" % str(e))
    if verbosity >= 2:
        sys.stderr.write("\n")

    # Renumber the solutions.  Our MiniZinc code outputs them all as "Solution
    # #1" because I don't know how to get MiniZinc to number outputs itself.
    soln_re = re.compile(r'^Solution #(\d+)', re.MULTILINE)
    sn = 1
    sn_str = "Solution #%d" % sn
    for line in mzn_output.strip("\n").split("\n"):
        line, nsubs = soln_re.subn(sn_str, line)
        if nsubs > 0:
            sn += 1
            sn_str = "Solution #%d" % sn
        if line == " ":
            line = ""
        print(line)

    # Delete the .mzn file.
    os.remove(mzn_fname)
