#! /usr/bin/env python

##################################################
# D-Wave quantum machine instruction "assembler" #
# By Scott Pakin <pakin@lanl.gov>                #
##################################################

import qmasm
import os
import os.path
import shlex
import string
import sys


# Parse the command line.
cl_args = qmasm.parse_command_line()

# Specify the minimum distinguishable difference between energy readings.
min_energy_delta = 0.005

# Parse the original input file(s) into an internal representation.
fparse = qmasm.FileParser()
fparse.parse_files(cl_args.input)

# Parse the variable pinnings specified on the command line.  Append these to
# the program.
if cl_args.pin != None:
    for pin in cl_args.pin:
        qmasm.program.extend(qmasm.process_pin("[command line]", 1, pin))

# Walk the statements in the program, processing each in turn.
logical_either = qmasm.Problem(cl_args.qubo)
for stmt in qmasm.program:
    stmt.update_qmi("", "<ERROR>", logical_either)

# Store all tallies for later reportage.
logical_stats = {
    "vars":      qmasm.next_sym_num + 1,
    "strengths": len(logical_either.strengths),
    "eqs":       len(logical_either.chains),
    "pins":      len(logical_either.pinned)
}

# Convert from QUBO to Ising in case the solver doesn't support QUBO problems.
if cl_args.qubo:
    logical_ising = logical_either.convert_to_ising()
else:
    logical_ising = logical_either

# Define a strength for each user-specified chain, and assign strengths
# to those chains.
qmasm.chain_strength = logical_ising.assign_chain_strength(cl_args.chain_strength)

# Define a strength for each user-specified pinned variable.
qmasm.pin_strength = logical_ising.assign_pin_strength(cl_args.pin_strength, qmasm.chain_strength)

# Output the chain and pin strengths.
if cl_args.verbose >= 1:
    sys.stderr.write("Computed the following strengths:\n\n")
    sys.stderr.write("    chain: %7.4f\n" % qmasm.chain_strength)
    sys.stderr.write("    pin:   %7.4f\n" % qmasm.pin_strength)
    sys.stderr.write("\n")

# Use a helper bit to help pin values to true or false.
logical_ising.pin_qubits(qmasm.pin_strength, qmasm.chain_strength)

# Convert chains to aliases where possible.
if cl_args.O:
    # Say what we're about to do
    if cl_args.verbose >= 2:
        sys.stderr.write("Replaced chains of equally weighted qubits with aliases:\n\n")
        sys.stderr.write("  %6d logical qubits before optimization\n" % (qmasm.next_sym_num + 1))

    # Replace chains with aliases wherever we can.
    logical_ising.convert_chains_to_aliases()

    # Summarize what we just did.
    if cl_args.verbose >= 2:
        sys.stderr.write("  %6d logical qubits after optimization\n\n" % (qmasm.next_sym_num + 1))

# This is a good time to update our logical statistics.
logical_stats["vars"] = qmasm.next_sym_num + 1

# Complain if we have no weights and no strengths.
if len(logical_ising.weights) == 0 and len(logical_ising.strengths) == 0:
    qmasm.abend("Nothing to do (no weights or strengths specified)")

# Complain if we have disconnected qubits.
discon_syms = logical_ising.find_disconnected_variables()
if len(discon_syms) > 0:
    qmasm.abend("Disconnected variables encountered: %s" % " ".join(sorted(discon_syms)))

# Output a normalized input file.
if cl_args.verbose >= 2:
    # Weights
    sys.stderr.write("Canonicalized the input file:\n\n")
    for q in sorted(logical_ising.weights.keys()):
        sys.stderr.write("    q%d %.20g\n" % (q + 1, logical_ising.weights[q]))
    sys.stderr.write("\n")

    # Chains
    if len(logical_ising.chains) > 0:
        for qs in sorted(logical_ising.chains.keys()):
            sys.stderr.write("    q%d = q%d\n" % (qs[0] + 1, qs[1] + 1))
        sys.stderr.write("\n")

    # Strengths (those that are not chains)
    for qs in sorted(logical_ising.strengths.keys()):
        if qs not in logical_ising.chains:
            sys.stderr.write("    q%d q%d %.20g\n" % (qs[0] + 1, qs[1] + 1, logical_ising.strengths[qs]))
    sys.stderr.write("\n")

    # Map each canonicalized name to one or more original symbols.
    canon2syms = [[] for _ in range(len(qmasm.sym2num))]
    max_sym_name_len = 8
    for s, n in list(qmasm.sym2num.items()):
        canon2syms[n].append(s)
        max_sym_name_len = max(max_sym_name_len, len(repr(canon2syms[n])) - 1)

    # Output the mapping we just computed.
    sys.stderr.write("Constructed a key to the above:\n\n")
    sys.stderr.write("    Canonical  %-*s\n" % (max_sym_name_len, "Original"))
    sys.stderr.write("    ---------  %s\n" % ("-" * max_sym_name_len))
    for i in range(len(canon2syms)):
        if canon2syms[i] == []:
            continue
        name_list = " ".join(sorted(canon2syms[i]))
        sys.stderr.write("    q%-8d  %s\n" % (i + 1, name_list))
    sys.stderr.write("\n")

# Establish a connection to the D-Wave, and use this to talk to a solver.  We
# rely on the qOp infrastructure to set the environment variables properly.
qmasm.connect_to_dwave()

# Output most or all solver properties.
if cl_args.verbose >= 1:
    # Introduce a few extra solver properties.
    ext_solver_properties = {}
    try:
        L, M, N = qmasm.chimera_topology(qmasm.solver)
        ext_solver_properties["chimera_toplogy_M_N_L"] = [M, N, L]
    except KeyError:
        pass
    ext_solver_properties["solver_name"] = qmasm.solver_name
    try:
        ext_solver_properties["connection_name"] = os.environ["DW_INTERNAL__CONNECTION"]
    except KeyError:
        pass
    ext_solver_properties.update(qmasm.solver.properties)

    # Determine the width of the widest key.
    max_key_len = len("Parameter")
    solver_props = list(ext_solver_properties.keys())
    solver_props.sort()
    for k in solver_props:
        max_key_len = max(max_key_len, len(k))

    # Output either "short" values (if verbose = 1) or all values (if
    # verbose > 1).
    short_value_len = 70 - max_key_len
    sys.stderr.write("Encountered the following solver properties:\n\n")
    sys.stderr.write("    %-*s  Value\n" % (max_key_len, "Parameter"))
    sys.stderr.write("    %s  -----\n" % ("-" * max_key_len))
    for k in solver_props:
        val_str = repr(ext_solver_properties[k])
        if cl_args.verbose >= 2 or len(val_str) <= short_value_len:
            sys.stderr.write("    %-*s  %s\n" % (max_key_len, k, val_str))
    sys.stderr.write("\n")

# As a special case, if the user specified --format=qbsolv we work with the
# pre-embedded version of the problem then exit.
if cl_args.format == "qbsolv":
    if cl_args.run:
        qmasm.run_qbsolv(logical_ising, cl_args.output,
                         shlex.split(cl_args.extra_args), cl_args.verbose)
    else:
        qmasm.write_output(logical_ising, cl_args.output, cl_args.format, cl_args.qubo)
    sys.exit(0)

# As a special case, if the user specified --format=minizinc we work with the
# pre-embedded version of the problem then exit.
if cl_args.format == "minizinc":
    if cl_args.run:
        qmasm.run_minizinc(logical_ising, cl_args.output, shlex.split(cl_args.extra_args), cl_args.verbose)
    else:
        qmasm.write_output(logical_ising, cl_args.output, cl_args.format, cl_args.qubo)
    sys.exit(0)

# As a special case, if the user requested QMASM output we output the
# pre-embedded version of the problem.  If we weren't asked to run, we can stop
# here.
if cl_args.format == "qmasm":
    qmasm.write_output(logical_ising, cl_args.output, cl_args.format, cl_args.qubo)
    if not cl_args.run:
        sys.exit(0)

# Embed the problem onto the D-Wave.
physical = qmasm.embed_problem_on_dwave(logical_ising, cl_args.O, cl_args.verbose)

# Set all chains to the user-specified strength then combine user-specified
# chains with embedder-created chains.
physical = qmasm.update_strengths_from_chains(physical)
if cl_args.verbose >= 2:
    sys.stderr.write("Introduced the following new chains:\n\n")
    if len(physical.chains) == 0:
        sys.stderr.write("    [none]\n")
    else:
        for c in physical.chains:
            num1, num2 = c
            if num1 > num2:
                num1, num2 = num2, num1
            sys.stderr.write("    %4d = %4d\n" % (num1, num2))
    sys.stderr.write("\n")

# Map each logical qubit to one or more symbols.
num2syms = [[] for _ in range(len(qmasm.sym2num))]
max_sym_name_len = 7
for s, n in list(qmasm.sym2num.items()):
    if cl_args.verbose >= 2 or "$" not in s:
        num2syms[n].append(s)
    max_sym_name_len = max(max_sym_name_len, len(repr(num2syms[n])) - 1)

# Output the embedding.
if cl_args.verbose >= 1:
    sys.stderr.write("Established a mapping from logical to physical qubits:\n\n")
    sys.stderr.write("    Logical  %-*s  Physical\n" % (max_sym_name_len, "Name(s)"))
    sys.stderr.write("    -------  %s  --------\n" % ("-" * max_sym_name_len))
    for i in range(len(physical.embedding)):
        if num2syms[i] == []:
            continue
        name_list = " ".join(sorted(num2syms[i]))
        phys_list = " ".join(["%4d" % e for e in sorted(physical.embedding[i])])
        sys.stderr.write("    %7d  %-*s  %s\n" % (i, max_sym_name_len, name_list, phys_list))
    sys.stderr.write("\n")
else:
    # Even at zero verbosity, we still note the logical-to-physical mapping.
    log2phys_comments = []
    for i in range(len(physical.embedding)):
        if num2syms[i] == []:
            continue
        name_list = " ".join(num2syms[i])
        phys_list = " ".join(["%d" % e for e in sorted(physical.embedding[i])])
        log2phys_comments.append("# %s --> %s" % (name_list, phys_list))
    log2phys_comments.sort()
    sys.stderr.write("\n".join(log2phys_comments) + "\n")

# Output some statistics about the embedding.
if cl_args.verbose >= 1:
    # Output a table.
    phys_wts = [elt for lst in physical.embedding for elt in lst]
    sys.stderr.write("Computed the following statistics of the logical-to-physical mapping:\n\n")
    sys.stderr.write("    Type      Metric          Value\n")
    sys.stderr.write("    --------  --------------  -----\n")
    sys.stderr.write("    Logical   Variables       %5d\n" % logical_stats["vars"])
    sys.stderr.write("    Logical   Strengths       %5d\n" % logical_stats["strengths"])
    sys.stderr.write("    Logical     Equivalences  %5d\n" % logical_stats["eqs"])
    sys.stderr.write("    Logical     Pins          %5d\n" % logical_stats["pins"])
    sys.stderr.write("    Physical  Qubits          %5d\n" % len(phys_wts))
    sys.stderr.write("    Physical  Couplers        %5d\n" % len(physical.strengths))
    sys.stderr.write("    Physical    Chains        %5d\n" % len(physical.chains))
    sys.stderr.write("\n")

    # Output some additional chain statistics.
    chain_lens = [len(c) for c in physical.embedding]
    max_chain_len = 0
    if chain_lens != []:
        max_chain_len = max(chain_lens)
    num_max_chains = len([l for l in chain_lens if l == max_chain_len])
    sys.stderr.write("    Maximum chain length = %d (occurrences = %d)\n\n" % (max_chain_len, num_max_chains))

# Manually scale the weights and strengths so Qubist doesn't complain.
physical = qmasm.scale_weights_strengths(physical, cl_args.verbose)

# Output a file in any of a variety of formats.  Note that we've already
# handled qbsolv output and QMASM output as special cases.
if cl_args.format != "qmasm" and (cl_args.output != "<stdout>" or not cl_args.run):
    qmasm.write_output(physical, cl_args.output, cl_args.format, cl_args.qubo)

# If we weren't told to run anything we can exit now.
if not cl_args.run:
    sys.exit(0)

# Submit the problem to the D-Wave.
if cl_args.verbose >= 1:
    sys.stderr.write("Submitting the problem to the %s solver.\n\n" % qmasm.solver_name)
dwave_response = qmasm.submit_dwave_problem(cl_args.verbose,
                                            physical,
                                            cl_args.samples,
                                            cl_args.anneal_time,
                                            cl_args.spin_revs,
                                            cl_args.postproc,
                                            cl_args.discard)
answer, final_answer, num_occurrences, num_not_broken = dwave_response

# Output solver timing information.
if cl_args.verbose >= 1:
    try:
        timing_info = list(answer["timing"].items())
        sys.stderr.write("Timing information:\n\n")
        sys.stderr.write("    %-30s %-10s\n" % ("Measurement", "Value (us)"))
        sys.stderr.write("    %s %s\n" % ("-" * 30, "-" * 10))
        for timing_value in timing_info:
            sys.stderr.write("    %-30s %10d\n" % timing_value)
        sys.stderr.write("\n")
    except KeyError:
        # Not all solvers provide timing information.
        pass

# Define a class to represent a valid solution.
class ValidSolution:
    "Represent a minimal state of a spin system."

    def __init__(self, soln, energy):
        self.solution = soln
        self.energy = energy
        self.names = []   # List of names for each named row
        self.spins = []   # Spin for each named row
        self.id = 0       # Map from spins to an int
        for q in range(len(soln)):
            if num2syms[q] == []:
                continue
            self.names.append(" ".join(num2syms[q]))
            self.spins.append(soln[q])
            self.id = self.id*2 + soln[q]

# Determine the set of solutions to output.
energies = answer["energies"]
n_low_energies = len([e for e in energies if abs(e - energies[0]) < min_energy_delta])
if cl_args.all_solns:
    n_solns_to_output = len(final_answer)
else:
    n_solns_to_output = min(n_low_energies, len(final_answer))
id2solution = {}   # Map from an int to a solution
for snum in range(n_solns_to_output):
    soln = ValidSolution(final_answer[snum], energies[snum])
    if soln.id not in id2solution:
        id2solution[soln.id] = soln

# Output information about the raw solutions.
if cl_args.verbose >= 1:
    sys.stderr.write("Number of solutions found:\n\n")
    sys.stderr.write("    %6d total\n" % len(energies))
    sys.stderr.write("    %6d with no broken chains or broken pins\n" % num_not_broken)
    sys.stderr.write("    %6d at minimal energy\n" % n_low_energies)
    sys.stderr.write("    %6d excluding duplicate variable assignments\n" % len(id2solution))
    sys.stderr.write("\n")

# Output energy tallies.  We first recompute these because some entries seem to
# be multiply listed.
if cl_args.verbose >= 2:
    try:
        tallies = answer["num_occurrences"]
    except KeyError:
        tallies = [1] * len(energies)
    new_energy_tallies = {}
    for i in range(len(energies)):
        e = float(energies[i])
        t = int(tallies[i])
        try:
            new_energy_tallies[e] += t
        except KeyError:
            new_energy_tallies[e] = t
    new_energies = list(new_energy_tallies.keys())
    new_energies.sort()
    min_energy_possible = -sum([abs(w) for w in physical.weights] + [abs(s) for s in list(physical.strengths.values())])
    sys.stderr.write("Energy histogram (theoretical minimum = %.4f):\n\n" % min_energy_possible)
    sys.stderr.write("    Energy      Tally\n")
    sys.stderr.write("    ----------  ------\n")
    for e in new_energies:
        sys.stderr.write("    %10.4f  %6d\n" % (e, new_energy_tallies[e]))
    sys.stderr.write("\n")

# Output the solution to the standard output device.
qmasm.output_solution(id2solution, num_occurrences, cl_args.values)
