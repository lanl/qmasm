##############################################
# Run various classical (non-D-Wave) solvers #
# By Scott Pakin <pakin@lanl.gov>            #
##############################################

import os
import qmasm
import subprocess
import sys
import tempfile

def run_qbsolv(logical_ising, oname, extra_args):
    "Run qmasm-qbsolv on the problem and report the result."
    # Use the specified file name if provided.  Otherwise, write to a temporary
    # file."
    qubo_fname = oname
    if qubo_fname == "<stdout>":
        tfile = tempfile.NamedTemporaryFile(suffix=".qubo", prefix="qmasm-", mode="w+")
        qubo_fname = tfile.name
        tfile.close()

    # Create the .qubo file.
    qmasm.write_output(logical_ising, qubo_fname, "qbsolv", False)

    # Run qmasm-qbsolv on the .qubo file.
    args = ["qmasm-qbsolv", "-i", qubo_fname] + extra_args
    try:
        subprocess.call(args, stdout=sys.stdout, stderr=sys.stderr)
    except OSError as e:
        qmasm.abend("Failed to run qmasm-qbsolv (%s)" % str(e.args[1]))

    # Delete the .qubo file if it's considered temporary.
    if oname == "<stdout>":
        os.remove(qubo_fname)
