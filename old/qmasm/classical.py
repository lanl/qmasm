##############################################
# Run various classical (non-D-Wave) solvers #
# By Scott Pakin <pakin@lanl.gov>            #
##############################################

import os
import qmasm
import re
import shlex
import subprocess
import sys
import tempfile

# Define a scale factor for converting floats to ints for MiniZinc's sake.
qmasm.minizinc_scale_factor = 10000.0

def run_qbsolv(ising, oname, extra_args, verbosity):
    "Run qmasm-qbsolv on the problem and report the result."
    # Use the specified file name if provided.  Otherwise, write to a temporary
    # file.
    if oname == "<stdout>":
        tfile = tempfile.NamedTemporaryFile(suffix=".qubo", prefix="qmasm-", mode="w+")
        qubo_fname = tfile.name
        tfile.close()
    else:
        qubo_fname = oname

    # Create the .qubo file.
    qmasm.write_output(ising, qubo_fname, "qbsolv", False)

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

def run_minizinc(ising, oname, extra_args, verbosity):
    "Run mzn-chuffed on the problem and report the result."
    # Always create a temporary file because we need something we can overwrite.
    tfile = tempfile.NamedTemporaryFile(suffix=".mzn", prefix="qmasm-", mode="w+")
    mzn_fname = tfile.name
    tfile.close()
    qmasm.write_output(ising, mzn_fname, "minizinc", False)

    # Patch the energy output to show the raw (integer) value.
    with open(mzn_fname) as tfile:
        code = re.sub(r'show\(energy.*?\)', "show(energy)", tfile.read())
    tfile = open(mzn_fname, "w")
    tfile.write(code)
    tfile.close()

    # Run MiniZinc on the .mzn file.
    args = ["minizinc", "--flatzinc-cmd=fzn-chuffed"] + extra_args + [mzn_fname]
    if verbosity >= 1:
        sys.stderr.write("Submitting the problem to MiniZinc.\n\n")
        if verbosity >= 2:
            sys.stderr.write("    Command line: %s\n" % str(" ".join(args)))
    try:
        mzn_output = subprocess.check_output(args, stderr=sys.stderr, universal_newlines=True)
    except Exception as e:
        qmasm.abend("Failed to run MiniZinc (%s)" % str(e))

    # Extract the energy value.
    match = re.search(r'\benergy = (-?\d+),', mzn_output)
    if match == None:
        qmasm.abend("Failed to find an energy level in MiniZinc output")
    energy = int(match.group(1))

    # Regenerate the MiniZinc input as a satisfaction problem instead of a
    # minimization problem.
    tfile = open(mzn_fname, "w")
    qmasm.output_minizinc(tfile, ising, energy)
    tfile.close()

    # If the user specified a file, write that, too.
    if oname != "<stdout>":
        tfile = open(oname, "w")
        qmasm.output_minizinc(tfile, ising, energy)
        tfile.close()

    # Run MiniZinc on the modified .mzn file.
    args = ["minizinc", "--flatzinc-cmd=fzn-chuffed"]
    args += ["--all-solutions", "--solution-separator= ", "--search-complete-msg="]
    args += extra_args + [mzn_fname]
    if verbosity >= 2:
        sys.stderr.write("    Command line: %s\n" % str(" ".join(args)))
    try:
        mzn_output = subprocess.check_output(args, stderr=sys.stderr, universal_newlines=True)
    except Exception as e:
        qmasm.abend("Failed to run MiniZinc (%s)" % str(e))
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

def process_classical(ising, cl_args):
    "Write a file for classical solution and optionally run it."
    # Extract variables we need from the command-line arguments,
    format = cl_args.format
    oname = cl_args.output
    run = cl_args.run
    extra_args = cl_args.extra_args
    as_qubo = cl_args.qubo
    verbosity = cl_args.verbose
    style = cl_args.values

    # The only formats we need to process at this point are qbsolv and minizinc.
    if format == "qbsolv":
        if run:
            # Pass both user-supplied and computed arguments to qbsolv.
            args = shlex.split(extra_args)
            args.append("--verbosity=%d" % verbosity)
            args.append("--values=%s" % style)
            args.append("--pin-weight=%f" % qmasm.pin_weight)
            args.append("--chain-strength=%f" % qmasm.chain_strength)
            qmasm.run_qbsolv(ising, oname, args, verbosity)
        else:
            qmasm.write_output(ising, oname, format, as_qubo)
    elif format == "minizinc":
        if run:
            qmasm.run_minizinc(ising, oname, shlex.split(extra_args), verbosity)
        else:
            qmasm.write_output(ising, oname, format, as_qubo)
