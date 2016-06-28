#! /usr/bin/env python

##################################################
# D-Wave quantum machine instruction "assembler" #
# By Scott Pakin <pakin@lanl.gov>                #
##################################################

from qasm.cmdline import parse_command_line
from qasm.globals import *
from qasm.utils import *
from collections import defaultdict
from dwave_sapi2.core import solve_ising
from dwave_sapi2.embedding import find_embedding, embed_problem, unembed_answer
from dwave_sapi2.local import local_connection
from dwave_sapi2.remote import RemoteConnection
from dwave_sapi2.util import get_hardware_adjacency, qubo_to_ising, ising_to_qubo, linear_index_to_chimera
import os
import os.path
import math
import random
import re
import shlex
import string
import sys

# Parse the command line.
cl_args = parse_command_line()

# Define our internal representation.
weights = defaultdict(lambda: 0.0)    # Map from a spin to a point weight
strengths = defaultdict(lambda: 0.0)  # Map from a pair of spins to a coupler strength
chains = {}     # Subset of strengths keys that represents chains
pinned = []     # Pairs of {unique number, Boolean} to pin

# Define synonyms for "true" and "false".
str2bool = {s: True for s in ["1", "+1", "T", "TRUE"]}
str2bool.update({s: False for s in ["0", "-1", "F", "FALSE"]})

# Specify the minimum distinguishable difference between energy readings.
min_energy_delta = 0.005

# Define a function that creates a new internal symbol.
def new_internal_sym():
    while True:
        sym = "$"
        for i in range(5):
            sym += random.choice(string.lowercase)
        if not sym2num.has_key(sym):
            return sym

# Define a function that says if a string can be treated as a float.
def is_float(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

# Define a function that aborts the program, reporting an invalid
# input line as part of the error message.
filename = "<stdin>"
lineno = 0
def error_in_line(str):
    sys.stderr.write('%s:%d: error: %s\n' % (filename, lineno, str))
    sys.exit(1)

# Define a function that searches a list of directories for a file.
def find_file_in_path(pathnames, filename):
    for pname in pathnames:
        fname = os.path.join(pname, filename)
        if os.path.exists(fname):
            return fname
        fname_qasm = fname + ".qasm"
        if os.path.exists(fname_qasm):
            return fname_qasm
    return None

# Define a function that returns the topology of the chimera graph associated
# with a given solver.
def chimera_topology(solver):
    nominal_qubits = solver.properties["num_qubits"]
    couplers = solver.properties["couplers"]
    deltas = [abs(c1 - c2) for c1, c2 in couplers]
    delta_tallies = {d: 0 for d in deltas}
    for d in deltas:
        delta_tallies[d] += 1
        sorted_tallies = sorted(delta_tallies.items(), key=lambda dt: dt[1], reverse=True)
    L = sorted_tallies[0][0]
    M = sorted_tallies[1][0] // (2*L)
    N = (nominal_qubits + 2*L*M - 1) // (2*L*M)
    return L, M, N

# Define a function that maps a pair of qubits to a coupler number a la the dw
# command.
def coupler_number(M, N, L, q1, q2):
    qmin = min(q1, q2)
    qmax = max(q1, q2)
    [[imin, jmin, umin, kmin], [imax, jmax, umax, kmax]] = linear_index_to_chimera([qmin, qmax], M, N, L)
    cell_links = L*L
    if imin == imax and jmin == jmax and umin != umax:
        # Same unit cell
        return cell_links*(imin*M + jmin) + kmin*L + kmax
    total_intra = cell_links*M*N
    if imin == imax and jmin + 1 == jmax and umin == umax and kmin == kmax:
        # Horizontal (same cell row)
        return total_intra + L*(imin*(M - 1) + jmin) + kmin
    total_horiz = (M - 1)*N*L
    if imin + 1 == imax and jmin == jmax and umin == umax and kmin == kmax:
        # Vertical (same cell column)
        return total_intra + total_horiz + L*(imin*M + jmin) + kmin
    raise IndexError("No coupler exists between Q%04d and Q%04d" % (q1, q2))

class Statement(object):
    "One statement in a QASM source file."

    def __init__(self, lineno):
        self.lineno = lineno

    def error_in_line(self, msg):
        if self.lineno == None:
            abend(msg)
        else:
            sys.stderr.write('%s:%d: error: %s\n' % (cl_args.input, self.lineno, msg))
        sys.exit(1)

class Weight(Statement):
    "Represent a point weight on a qubit."
    def __init__(self, lineno, sym, weight):
        super(Weight, self).__init__(lineno)
        self.sym = sym
        self.weight = weight

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        num = symbol_to_number(prefix + self.sym)
        weights[num] += self.weight

class Chain(Statement):
    "Chain between qubits."
    def __init__(self, lineno, sym1, sym2):
        super(Chain, self).__init__(lineno)
        self.sym1 = sym1
        self.sym2 = sym2

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        num1 = symbol_to_number(prefix + self.sym1)
        num2 = symbol_to_number(prefix + self.sym2)
        if num1 == num2:
            self.error_in_line("A chain cannot connect a spin to itself")
        elif num1 > num2:
            num1, num2 = num2, num1
        chains[(num1, num2)] = None   # Value is a don't-care.

class Pin(Statement):
    "Pinning of a qubit to true or false."
    def __init__(self, lineno, sym, goal):
        super(Pin, self).__init__(lineno)
        self.sym = sym
        self.goal = goal

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        num = symbol_to_number(prefix + self.sym)
        pinned.append((num, self.goal))

class Alias(Statement):
    "Alias of one symbol to another."
    def __init__(self, lineno, sym1, sym2):
        super(Alias, self).__init__(lineno)
        self.sym1 = sym1
        self.sym2 = sym2

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        sym1 = prefix + self.sym1
        sym2 = prefix + self.sym2
        try:
            sym2num[sym1] = sym2num[sym2]
        except KeyError:
            self.error_in_line("Cannot make symbol %s an alias of undefined symbol %s" % (sym1, sym2))
        if sym1 == sym2:
            self.error_in_line("Fields cannot alias themselves")

class Strength(Statement):
    "Coupler strength between two qubits."
    def __init__(self, lineno, sym1, sym2, strength):
        super(Strength, self).__init__(lineno)
        self.sym1 = sym1
        self.sym2 = sym2
        self.strength = strength

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        num1 = symbol_to_number(prefix + self.sym1)
        num2 = symbol_to_number(prefix + self.sym2)
        if num1 == num2:
            self.error_in_line("A coupler cannot connect a spin to itself")
        elif num1 > num2:
            num1, num2 = num2, num1
        strengths[(num1, num2)] += self.strength

class MacroUse(Statement):
    "Instantiation of a macro definition."
    def __init__(self, lineno, name, body, prefix):
        super(MacroUse, self).__init__(lineno)
        self.name = name
        self.body = body
        self.prefix = prefix

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        for stmt in self.body:
            stmt.update_qmi(prefix + self.prefix, weights, strengths, chains, pinned)

# Define a function that parses an input file into an internal representation.
# This function can be called recursively (due to !include directives).
program = []       # List of Statement objects
macros = {}        # Map from a macro name to a list of Statement objects
current_macro = (None, [])   # Macro currently being defined (name and statements)
aliases = {}       # Map from a symbol to its textual expansion
target = program   # Reference to either the program or the current macro
def parse_file(infilename, infile):
    global program, macros, current_macro, aliases, target, filename, lineno
    filename = infilename
    for line in infile:
        # Split the line into fields and apply text aliases.
        lineno += 1
        if line.strip() == "":
            continue
        fields = shlex.split(line, True)
        for i in range(len(fields)):
            try:
                fields[i] = aliases[fields[i]]
            except KeyError:
                pass

        # Process the line.
        if len(fields) == 0:
            # Ignore empty lines.
            continue
        elif len(fields) >= 2 and fields[0] == "!include":
            # "!include" "<filename>" -- process a named auxiliary file.
            incname = string.join(fields[1:], " ")
            if len(incname) >= 2 and incname[0] == "<" and incname[-1] == ">":
                # Search QASMPATH for the filename.
                incname = incname[1:-1]
                try:
                    qasmpath = string.split(os.environ["QASMPATH"], ":")
                    qasmpath.append(".")
                except KeyError:
                    qasmpath = ["."]
                found_incname = find_file_in_path(qasmpath, incname)
                if found_incname != None:
                    incname = found_incname
            elif len(incname) >= 2:
                # Search only the current directory for the filename.
                found_incname = find_file_in_path(["."], incname)
                if found_incname != None:
                    incname = found_incname
            try:
                incfile = open(incname)
            except IOError:
                abend('Failed to open %s for input' % incname)
            parse_file(incname, incfile)
            incfile.close()
        elif len(fields) == 2:
            if fields[0] == "!begin_macro":
                # "!begin_macro" <name> -- begin a macro definition.
                name = fields[1]
                if macros.has_key(name):
                    error_in_line("Macro %s is multiply defined" % name)
                if current_macro[0] != None:
                    error_in_line("Nested macros are not supported")
                current_macro = (name, [])
                target = current_macro[1]
            elif fields[0] == "!end_macro":
                # "!end_macro" <name> -- end a macro definition.
                name = fields[1]
                if current_macro[0] == None:
                    error_in_line("Ended macro %s with no corresponding begin" % name)
                if current_macro[0] != name:
                    error_in_line("Ended macro %s after beginning macro %s" % (name, current_macro[0]))
                macros[name] = current_macro[1]
                target = program
                current_macro = (None, [])
            else:
                # <symbol> <weight> -- increment a symbol's point weight.
                try:
                    val = float(fields[1])
                except ValueError:
                    error_in_line('Failed to parse "%s %s" as a symbol followed by a numerical weight' % (fields[0], fields[1]))
                target.append(Weight(lineno, fields[0], val))
        elif len(fields) == 3:
            if fields[1] == "=":
                # <symbol_1> = <symbol_2> -- create a chain between <symbol_1> and
                # <symbol_2>.
                target.append(Chain(lineno, fields[0], fields[2]))
            elif fields[1] == ":=":
                # <symbol> := <value> -- force symbol <symbol> to have value
                # <value>.
                try:
                    goal = str2bool[fields[2].upper()]
                except KeyError:
                    error_in_line('Right-hand side ("%s") must be a Boolean value' % fields[2])
                target.append(Pin(lineno, fields[0], goal))
            elif fields[1] == "<->":
                # <symbol_1> <-> <symbol_2> -- make <symbol_1> an alias of
                # <symbol_2>.
                target.append(Alias(lineno, fields[0], fields[2]))
            elif fields[0] == "!use_macro":
                # "!use_macro" <macro_name> <instance_name> -- instantiate
                # a macro using <instance_name> as each variable's prefix.
                name = fields[1]
                try:
                    target.append(MacroUse(lineno, name, macros[name], fields[2] + "."))
                except KeyError:
                    error_in_line("Unknown macro %s" % name)
            elif fields[0] == "!alias":
                # "!alias" <symbol> <text> -- replace a field of <symbol> with
                # <text>.
                aliases[fields[1]] = fields[2]
            elif is_float(fields[2]):
                # <symbol_1> <symbol_2> <strength> -- increment a coupler strength.
                try:
                    strength = float(fields[2])
                except ValueError:
                    error_in_line('Failed to parse "%s" as a number' % fields[2])
                target.append(Strength(lineno, fields[0], fields[1], strength))
            else:
                # Three fields but none of the above cases
                error_in_line('Cannot parse "%s"' % line)
        else:
            # Neither two nor three fields
            error_in_line('Cannot parse "%s"' % line)

# Parse the original input file(s) into an internal representation.
if cl_args.input == []:
    # No files were specified: Read from standard input.
    parse_file("<stdin>", sys.stdin)
    if current_macro[0] != None:
        error_in_line("Unterminated definition of macro %s" % current_macro[0])
else:
    # Files were specified: Process each in turn.
    for infilename in cl_args.input:
        try:
            infile = open(infilename)
        except IOError:
            abend('Failed to open %s for input' % infilename)
        parse_file(infilename, infile)
        if current_macro[0] != None:
            error_in_line("Unterminated definition of macro %s" % current_macro[0])
        infile.close()

# Parse the variable pinnings specified on the command line.  Append these to
# the program.
if cl_args.pin != None:
    for pstr in cl_args.pin:
        lhs_rhs = pstr.split(":=")
        if len(lhs_rhs) != 2:
            abend('Failed to parse --pin="%s"' % pstr)
        lhs = lhs_rhs[0].split()
        rhs = []
        for r in lhs_rhs[1].upper().split():
            try:
                rhs.append(str2bool[r])
            except KeyError:
                for subr in r:
                    try:
                        rhs.append(str2bool[subr])
                    except KeyError:
                        abend('Failed to parse --pin="%s"' % pstr)
            if len(lhs) != len(rhs):
                abend('Different number of left- and right-hand-side values in --pin="%s" (%d vs. %d)' % (pstr, len(lhs), len(rhs)))
            for l, r in zip(lhs, rhs):
                program.append(Pin(None, l, r))

# Walk the statements in the program, processing each in turn.
for stmt in program:
    stmt.update_qmi("", weights, strengths, chains, pinned)

# Store all tallies for later reportage.
logical_stats = {
    "vars":      next_sym_num + 1,
    "strengths": len(strengths),
    "eqs":       len(chains),
    "pins":      len(pinned)
}

# Convert from QUBO to Ising in case the solver doesn't support QUBO problems.
if cl_args.qubo:
    qmatrix = {(q, q): w for q, w in weights.items()}
    qmatrix.update(strengths)
    hvals, strengths, _ = qubo_to_ising(qmatrix)
    strengths = defaultdict(lambda: 0.0, strengths)
    weights.update({i: hvals[i] for i in range(len(hvals))})

# Define a strength for each user-specified chain, and assign strengths
# to those chains.
chain_strength = cl_args.chain_strength
if chain_strength == None:
    # Chain strength defaults to the maximum strength in the data.
    try:
        chain_strength = -max([abs(w) for w in strengths.values()])
    except ValueError:
        # No strengths -- use weights instead.
        try:
            chain_strength = -max([abs(w) for w in weights.values()])
        except ValueError:
            # No weights or strengths -- arbitrarily choose -1.
            chain_strength = -1.0
elif cl_args.qubo:
    # With QUBO input we need to divide the chain strength by 4 for consistency
    # with the other coupler strengths.
    chain_strength /= 4.0
for c in chains.keys():
    strengths[c] += chain_strength

# Define a strength for each user-specified pinned variable.
pin_strength = cl_args.pin_strength
if pin_strength == None:
    # Pin strength defaults to half the chain strength.
    pin_strength = chain_strength/2.0
elif cl_args.qubo:
    # With QUBO input we need to divide the chain strength by 4 for consistency
    # with the other coupler strengths.
    pin_strength /= 4.0

# Output the chain and pin strengths.
if cl_args.verbose >= 1:
    sys.stderr.write("Computed the following strengths:\n\n")
    sys.stderr.write("    chain: %7.4f\n" % chain_strength)
    sys.stderr.write("    pin:   %7.4f\n" % pin_strength)
    sys.stderr.write("\n")

# Use a helper bit to help pin values to true or false.  A helper bit isn't
# strictly necessary, but I'm thinking it'll emphasize the pin more than a
# single point weight.  We recycle the coupler strength to set the pin
# strength.
for q_user, b in pinned:
    q_helper = symbol_to_number(new_internal_sym())
    q1, q2 = q_helper, q_user
    if q1 > q2:
        q1, q2 = q2, q1
    if b:
        weights[q_helper] += pin_strength/2.0
        weights[q_user] += pin_strength
        strengths[(q1, q2)] += -pin_strength/2.0
    else:
        weights[q_helper] += pin_strength/2.0
        weights[q_user] += -pin_strength
        strengths[(q1, q2)] += pin_strength/2.0

# Convert chains to aliases where possible.
if cl_args.O:
    # Say what we're about to do
    if cl_args.verbose >= 2:
        sys.stderr.write("Replaced chains of equally weighted qubits with aliases:\n\n")
        sys.stderr.write("  %6d logical qubits before optimization\n" % (next_sym_num + 1))

    # Identify all chains that can be converted to aliases.
    num2allsyms = [[] for _ in range(len(sym2num))]
    for s, n in sym2num.items():
        num2allsyms[n].append(s)
    make_aliases = []
    for q1, q2 in chains:
        if weights[q1] == weights[q2]:
            make_aliases.append((q1, q2))
    make_aliases.sort(reverse=True, key=lambda qs: (qs[1], qs[0]))

    # Replace each chain in make_aliases with an alias.  Work in reverse order
    # of qubit number and shift all greater qubit numbers downward.
    for q1, q2 in make_aliases:
        # Map q2's symbolic names to q1's.  Shift everything above q2 downwards.
        alias_sym2num = {}
        for s, sq in sym2num.items():
            if sq == q2:
                sq = q1
            elif sq > q2:
                sq -= 1
            alias_sym2num[s] = sq
        sym2num = alias_sym2num

        # Elide q2 from the list of weights.
        alias_weights = defaultdict(lambda: 0.0)
        for wq, wt in weights.items():
            if wq == q2:
                continue
            if wq > q2:
                wq -= 1
            if wt != 0.0:
                alias_weights[wq] = wt
        alias_weights[q1] += weights[q2]   # Conserve overall energy.
        weights = alias_weights

        # Replace q2 with q1 in all strengths.  Shift everything above q2
        # downwards.
        alias_strengths = defaultdict(lambda: 0.0)
        for (sq1, sq2), wt in strengths.items():
            if sq1 == q2:
                sq1 = q1
            if sq1 > q2:
                sq1 -= 1
            if sq2 == q2:
                sq2 = q1
            if sq2 > q2:
                sq2 -= 1
            if sq1 != sq2:
                alias_strengths[(sq1, sq2)] = wt
        strengths = alias_strengths

        # Replace q2 with q1 in all strengths.  Shift everything above q2
        # downwards.
        alias_chains = {}
        for cq1, cq2 in chains.keys():
            if cq1 == q2:
                cq1 = q1
            if cq1 > q2:
                cq1 -= 1
            if cq2 == q2:
                cq2 = q1
            if cq2 > q2:
                cq2 -= 1
            if cq1 != cq2:
                alias_chains[(cq1, cq2)] = None
        chains = alias_chains

        # Replace q2 with q1 in all pinned qubits.  Shift everything above q2
        # downwards.
        alias_pinned = []
        for pq, b in pinned:
            if pq == q2:
                pq = q1
            if pq > q2:
                pq -= 1
            alias_pinned.append((pq, b))
        pinned = alias_pinned

        # We now have one fewer symbol.
        next_sym_num -= 1

    # Summarize what we just did.
    if cl_args.verbose >= 2:
        sys.stderr.write("  %6d logical qubits after optimization\n\n" % (next_sym_num + 1))

# This is a good time to update our logical statistics.
logical_stats["vars"] = next_sym_num + 1

# Complain if we have no weights and no strengths.
if len(weights) == 0 and len(strengths) == 0:
    abend("Nothing to do (no weights or strengths specified)")

# Output a normalized input file.
if cl_args.verbose >= 2:
    # Weights
    sys.stderr.write("Canonicalized the input file:\n\n")
    for q in sorted(weights.keys()):
        sys.stderr.write("    q%d %.20g\n" % (q + 1, weights[q]))
    sys.stderr.write("\n")

    # Chains
    if len(chains) > 0:
        for qs in sorted(chains.keys()):
            sys.stderr.write("    q%d = q%d\n" % (qs[0], qs[1]))
        sys.stderr.write("\n")

    # Strengths (those that are not chains)
    for qs in sorted(strengths.keys()):
        if not chains.has_key(qs):
            sys.stderr.write("    q%d q%d %.20g\n" % (qs[0] + 1, qs[1] + 1, strengths[qs]))
    sys.stderr.write("\n")

    # Map each canonicalized name to one or more original symbols.
    canon2syms = [[] for _ in range(len(sym2num))]
    max_sym_name_len = 8
    for s, n in sym2num.items():
        canon2syms[n].append(s)
        max_sym_name_len = max(max_sym_name_len, len(repr(canon2syms[n])) - 1)

    # Output the mapping we just computed.
    sys.stderr.write("Constructed a key to the above:\n\n")
    sys.stderr.write("    Canonical  %-*s\n" % (max_sym_name_len, "Original"))
    sys.stderr.write("    ---------  %s\n" % ("-" * max_sym_name_len))
    for i in range(len(canon2syms)):
        if canon2syms[i] == []:
            continue
        name_list = string.join(sorted(canon2syms[i]))
        sys.stderr.write("    q%-8d  %s\n" % (i + 1, name_list))
    sys.stderr.write("\n")

# Establish a connection to the D-Wave, and use this to talk to a solver.  We
# rely on the qOp infrastructure to set the environment variables properly.
try:
    url = os.environ["DW_INTERNAL__HTTPLINK"]
    token = os.environ["DW_INTERNAL__TOKEN"]
    conn = RemoteConnection(url, token)
except KeyError:
    url = "<local>"
    token = "<N/A>"
    conn = local_connection
except IOError as e:
    abend("Failed to establish a remote connection (%s)" % e)
try:
    solver_name = os.environ["DW_INTERNAL__SOLVER"]
except:
    # Solver was not specified: Use the first available solver.
    solver_name = conn.solver_names()[0]
try:
    solver = conn.get_solver(solver_name)
except KeyError:
    abend("Failed to find solver %s on connection %s" % (solver_name, url))

# Output most or all solver properties.
if cl_args.verbose >= 1:
    # Determine the width of the widest key.
    max_key_len = len("Parameter")
    ext_solver_properties = {}
    try:
        L, M, N = chimera_topology(solver)
        ext_solver_properties["chimera_toplogy_M_N_L"] = [M, N, L]
    except KeyError:
        pass
    ext_solver_properties.update(solver.properties)
    solver_props = ext_solver_properties.keys()
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

# Find an embedding of the problem.
edges = strengths.keys()
edges.sort()
try:
    hw_adj = get_hardware_adjacency(solver)
except KeyError:
    # The Ising heuristic solver is an example of a solver that lacks a fixed
    # hardware representation.  We therefore assert that the hardware exactly
    # matches the problem'input graph.
    hw_adj = edges
if cl_args.verbose >= 2:
    sys.stderr.write("Embedding the logical adjacency within the physical topology.\n\n")
    stdout_fd = os.dup(sys.stdout.fileno())
    os.dup2(sys.stderr.fileno(), sys.stdout.fileno())
    embedding = find_embedding(edges, hw_adj, verbose=1)
    os.dup2(stdout_fd, sys.stdout.fileno())
    sys.stderr.write("\n")
else:
    embedding = find_embedding(edges, hw_adj, verbose=0)
if embedding == []:
    # We received an empty embedding.  I've seen this happen with the
    # ising-heuristic solver.  A workaround seems to be to fabricate a trivial
    # embedding in which logical qubit X maps to physical qubit X.
    embedding = [[q] for q in range(next_sym_num + 1)]

# Embed the problem using the embedding we found.
try:
    h_range = solver.properties["h_range"]
    j_range = solver.properties["j_range"]
except KeyError:
    h_range = [-1.0, 1.0]
    j_range = [-1.0, 1.0]
weight_list = [weights[q] for q in range(next_sym_num + 1)]
smearable = any([s != 0.0 for s in strengths.values()])
try:
    [new_weights, new_strengths, new_chains, new_embedding] = embed_problem(
        weight_list, strengths, embedding, hw_adj, True, smearable, h_range, j_range)
except ValueError as e:
    abend("Failed to embed the problem in the solver (%s)" % e)

# If told to optimize the layout, iteratively search for a better embedding.
if cl_args.O:
    try:
        # Say what we're about to do
        if cl_args.verbose >= 2:
            sys.stderr.write("Optimizing the embedding.\n\n")

        # Determine the edges of a rectangle of cells we want to use.
        L, M, N = chimera_topology(solver)
        L2 = 2*L
        ncells = (next_sym_num + 1) // L2
        edgey = max(int(math.sqrt(ncells)), 1)
        edgex = max((ncells + edgey - 1) // edgey, 1)

        # Repeatedly expand edgex and edgey until the embedding works.
        while edgex < M and edgey < N:
            # Retain adjacencies only within the rectangle.
            alt_hw_adj = []
            for q1, q2 in hw_adj:
                c1 = q1//L2
                if c1 % M >= edgex:
                    continue
                if c1 // M >= edgey:
                    continue
                c2 = q2//L2
                if c2 % M >= edgex:
                    continue
                if c2 // M >= edgey:
                    continue
                alt_hw_adj.append((q1, q2))
            alt_hw_adj = set(alt_hw_adj)

            # Embed again.
            try:
                if cl_args.verbose >= 2:
                    sys.stderr.write("  Trying a %dx%d unit-cell embedding ... " % (edgex, edgey))
                alt_embedding = find_embedding(edges, alt_hw_adj, verbose=0)
                [new_weights, new_strengths, new_chains, new_embedding] = embed_problem(
                    weight_list, strengths, alt_embedding, alt_hw_adj, True, True, h_range, j_range)
                if cl_args.verbose >= 2:
                    sys.stderr.write("succeeded\n")
                break  # Success!
            except ValueError as e:
                # Failure -- increase edgex or edgey and try again.
                if cl_args.verbose >= 2:
                    sys.stderr.write("failed\n")
                if edgex < edgey:
                    edgex += 1
                else:
                    edgey += 1
    except KeyError:
        if cl_args.verbose >= 2:
            sys.stderr.write("  - Failed to query the machine topology for embedding parameters\n")
    if cl_args.verbose >= 2:
        sys.stderr.write("\n")

# Set all chains to the user-specified strength then combine
# user-specified chains with embedder-created chains.
new_chains = {c: chain_strength for c in new_chains.keys()}
combined_strengths = new_strengths.copy()
combined_strengths.update(new_chains)
if cl_args.verbose >= 2:
    sys.stderr.write("Introduced the following new chains:\n\n")
    if len(new_chains) == 0:
        sys.stderr.write("    [none]\n")
    else:
        for c in new_chains:
            num1, num2 = c
            if num1 > num2:
                num1, num2 = num2, num1
            sys.stderr.write("    %4d = %4d\n" % (num1, num2))
    sys.stderr.write("\n")

# Map each logical qubit to one or more symbols.
num2syms = [[] for _ in range(len(sym2num))]
max_sym_name_len = 7
for s, n in sym2num.items():
    if cl_args.verbose >= 2 or "$" not in s:
        num2syms[n].append(s)
    max_sym_name_len = max(max_sym_name_len, len(repr(num2syms[n])) - 1)

# Output the embedding.
if cl_args.verbose >= 1:
    sys.stderr.write("Established a mapping from logical to physical qubits:\n\n")
    sys.stderr.write("    Logical  %-*s  Physical\n" % (max_sym_name_len, "Name(s)"))
    sys.stderr.write("    -------  %s  --------\n" % ("-" * max_sym_name_len))
    for i in range(len(new_embedding)):
        if num2syms[i] == []:
            continue
        name_list = string.join(sorted(num2syms[i]))
        phys_list = string.join(["%4d" % e for e in sorted(new_embedding[i])])
        sys.stderr.write("    %7d  %-*s  %s\n" % (i, max_sym_name_len, name_list, phys_list))
    sys.stderr.write("\n")
else:
    # Even at zero verbosity, we still note the logical-to-physical mapping.
    log2phys_comments = []
    for i in range(len(new_embedding)):
        if num2syms[i] == []:
            continue
        name_list = string.join(num2syms[i])
        phys_list = string.join(["%d" % e for e in sorted(new_embedding[i])])
        log2phys_comments.append("# %s --> %s" % (name_list, phys_list))
    log2phys_comments.sort()
    sys.stderr.write("\n".join(log2phys_comments) + "\n")

# Output some statistics about the embedding.
if cl_args.verbose >= 1:
    # Output a table.
    phys_wts = [elt for lst in new_embedding for elt in lst]
    sys.stderr.write("Computed the following statistics of the logical-to-physical mapping:\n\n")
    sys.stderr.write("    Type      Metric          Value\n")
    sys.stderr.write("    --------  --------------  -----\n")
    sys.stderr.write("    Logical   Variables       %5d\n" % logical_stats["vars"])
    sys.stderr.write("    Logical   Strengths       %5d\n" % logical_stats["strengths"])
    sys.stderr.write("    Logical     Equivalences  %5d\n" % logical_stats["eqs"])
    sys.stderr.write("    Logical     Pins          %5d\n" % logical_stats["pins"])
    sys.stderr.write("    Physical  Qubits          %5d\n" % len(phys_wts))
    sys.stderr.write("    Physical  Couplers        %5d\n" % len(combined_strengths))
    sys.stderr.write("    Physical    Chains        %5d\n" % len(new_chains))
    sys.stderr.write("\n")

    # Output some additional chain statistics.
    chain_lens = [len(c) for c in new_embedding]
    max_chain_len = 0
    if chain_lens != []:
        max_chain_len = max(chain_lens)
    num_max_chains = len([l for l in chain_lens if l == max_chain_len])
    sys.stderr.write("    Maximum chain length = %d (occurrences = %d)\n\n" % (max_chain_len, num_max_chains))

# Manually scale the weights and strengths so Qubist doesn't complain.
old_cap = max([abs(w) for w in new_weights + combined_strengths.values()])
new_cap = min(-h_range[0], h_range[1], -j_range[0], j_range[1])
if old_cap == 0.0:
    # Handle the obscure case of a zero old_cap.
    old_cap = new_cap
new_weights = [w*new_cap/old_cap for w in new_weights]
combined_strengths = {js: w*new_cap/old_cap for js, w in combined_strengths.items()}
if cl_args.verbose >= 1 and old_cap != new_cap:
    sys.stderr.write("Scaling weights and strengths from [%.10g, %.10g] to [%.10g, %.10g].\n\n" % (-old_cap, old_cap, -new_cap, new_cap))

# Output a file in any of a variety of formats.
if cl_args.output != "<stdout>" or not cl_args.run:
    # Open the output file.
    if cl_args.output == "<stdout>":
        outfile = sys.stdout
    else:
        try:
            outfile = open(cl_args.output, "w")
        except IOError:
            abend('Failed to open %s for output' % cl_args.output)

    # Convert from Ising back to QUBO if --qubo was specified or if we're
    # outputting in dw format, which expects QUBO coefficients.
    if cl_args.qubo or cl_args.format == "dw":
        qmatrix, _ = ising_to_qubo(new_weights, combined_strengths)
        output_weights = [0] * len(new_weights)
        for (q1, q2), wt in qmatrix.items():
            if q1 == q2:
                output_weights[q1] = wt
        output_strengths = {(q1, q2): wt for (q1, q2), wt in qmatrix.items() if q1 != q2}
    else:
        output_weights = new_weights
        output_strengths = combined_strengths

    # Output the weights and strengths in the specified format.
    if cl_args.format == "qubist":
        # Output the data in Qubist format.
        data = []
        for q in range(len(new_weights)):
            if output_weights[q] != 0.0:
                data.append("%d %d %.10g" % (q, q, output_weights[q]))
        for sp, str in output_strengths.items():
            if str != 0.0:
                data.append("%d %d %.10g" % (sp[0], sp[1], str))

        # Output the header and data in Qubist format.
        try:
            num_qubits = solver.properties["num_qubits"]
        except KeyError:
            # The Ising heuristic solver is an example of a solver that lacks a
            # fixed hardware representation.  We therefore assert that the number
            # of qubits is exactly the number of qubits we require.
            num_qubits = len(new_weights)
        outfile.write("%d %d\n" % (num_qubits, len(data)))
        for d in data:
            outfile.write("%s\n" % d)
    elif cl_args.format == "dw":
        # Output weights and strengths in dw format.
        try:
            L, M, N = chimera_topology(solver)
        except KeyError:
            abend("Failed to query the chimera topology")
        wdata = []
        for q in range(len(new_weights)):
            if output_weights[q] != 0.0:
                wdata.append("Q%0d <== %.25g" % (q, output_weights[q]))
        wdata.sort()
        sdata = []
        for sp, str in output_strengths.items():
            if str != 0.0:
                coupler = coupler_number(M, N, L, sp[0], sp[1])
                sdata.append("C%04d <== %.25g" % (coupler, str))
        sdata.sort()
        outfile.write("\n".join(wdata + sdata) + "\n")
    elif cl_args.format == "qbsolv":
        # Output weights and strengths in qbsolv format.
        nonzero_strengths = [s for s in output_strengths.values() if s != 0.0]
        outfile.write("p qubo 0 %d %d %d\n" % (len(output_weights), len(output_weights), len(nonzero_strengths)))
        for q in range(len(output_weights)):
            outfile.write("%d %d %.10g\n" % (q, q, output_weights[q]))
        for qs in sorted(output_strengths.keys()):
            s = output_strengths[qs]
            if s != 0.0:
                outfile.write("%d %d %.10g\n" % (qs[0], qs[1], s))

    # Close the output file.
    if cl_args.output != "<stdout>":
        outfile.close()

# If we weren't told to run anything we can exit now.
if not cl_args.run:
    sys.exit(0)

# Submit the problem to the D-Wave.
if cl_args.verbose >= 1:
    sys.stderr.write("Submitting the problem to the %s solver.\n\n" % solver_name)
solver_params = dict(chains=new_embedding, num_reads=cl_args.samples, annealing_time=cl_args.anneal_time)
while True:
    # Repeatedly remove parameters the particular solver doesn't like until it
    # actually works -- or fails for a different reason.
    try:
        answer = solve_ising(solver, new_weights, combined_strengths, **solver_params)
        break
    except ValueError as e:
        # Is there a better way to extract the failing symbol than a regular
        # expression match?
        bad_name = re.match(r'"(.*?)"', str(e))
        if bad_name == None:
            raise e
        del solver_params[bad_name.group(1)]
    except RuntimeError as e:
        abend(e)
final_answer = unembed_answer(answer["solutions"], new_embedding, broken_chains="discard")

# Output solver timing information.
if cl_args.verbose >= 1:
    try:
        timing_info = answer["timing"].items()
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
    "Represent a valid ground state of a spin system."

    def __init__(self, soln, energy):
        self.solution = soln
        self.energy = energy
        self.names = []   # List of names for each named row
        self.spins = []   # Spin for each named row
        self.id = 0       # Map from spins to an int
        for q in range(len(soln)):
            if num2syms[q] == []:
                continue
            self.names.append(string.join(num2syms[q]))
            self.spins.append(soln[q])
            self.id = self.id*2 + soln[q]

# Define a function that says whether an answer contains no broken pins and no
# broken (user-specified) chains.
def answer_is_intact(answer):
    # Reject broken pins.
    bool2spin = [-1, +1]
    for pnum, pin in pinned:
        if answer[pnum] != bool2spin[pin]:
            return False

    # Reject broken chains.
    for q1, q2 in chains.keys():
        if answer[q1] != answer[q2]:
            return False

    # The answer looks good!
    return True

# Determine the set of solutions to output.
final_answer = [a for a in final_answer if answer_is_intact(a)]
if len(final_answer) == 0:
    print "No valid solutions found,"
    sys.exit(0)
energies = answer["energies"]
n_low_energies = len([e for e in energies if abs(e - energies[0]) < min_energy_delta])
if cl_args.verbose >= 2:
    n_solns_to_output = len(final_answer)
else:
    n_solns_to_output = min(n_low_energies, len(final_answer))
id2solution = {}   # Map from an int to a solution
for snum in range(n_solns_to_output):
    soln = ValidSolution(final_answer[snum], energies[snum])
    if not id2solution.has_key(soln.id):
        id2solution[soln.id] = soln

# Output information about the raw solutions.
if cl_args.verbose >= 1:
    sys.stderr.write("Number of solutions found:\n\n")
    sys.stderr.write("    %6d total\n" % len(energies))
    sys.stderr.write("    %6d with no broken chains or broken pins\n" % len(final_answer))
    sys.stderr.write("    %6d in the ground state\n" % n_low_energies)
    sys.stderr.write("    %6d excluding duplicate variable assignments\n" % len(id2solution))
    sys.stderr.write("\n")

# Output energy tallies.  We first recompute these because some entries seem to
# be multiply listed.
if cl_args.verbose >= 2:
    tallies = answer["num_occurrences"]
    new_energy_tallies = {}
    for i in range(len(energies)):
        e = float(energies[i])
        t = int(tallies[i])
        try:
            new_energy_tallies[e] += t
        except KeyError:
            new_energy_tallies[e] = t
    new_energies = new_energy_tallies.keys()
    new_energies.sort()
    sys.stderr.write("Energy histogram:\n\n")
    sys.stderr.write("    Energy      Tally\n")
    sys.stderr.write("    ----------  ------\n")
    for e in new_energies:
        sys.stderr.write("    %10.4f  %6d\n" % (e, new_energy_tallies[e]))
    sys.stderr.write("\n")

# Output to standard output (regardless of --output) the solutions with no
# broken chains.
bool_str = {-1: "False", +1: "True"}
sorted_solns = [id2solution[s] for s in sorted(id2solution.keys())]
for snum in range(len(sorted_solns)):
    print "Solution #%d (energy = %.2f):\n" % (snum + 1, sorted_solns[snum].energy)
    print "    %-*s  Spin  Boolean" % (max_sym_name_len, "Name(s)")
    print "    %s  ----  -------" % ("-" * max_sym_name_len)
    soln = sorted_solns[snum]
    output_lines = []
    for q in range(len(soln.spins)):
        names = soln.names[q]
        spin = soln.spins[q]
        output_lines.append("    %-*s  %+4d  %-7s" % (max_sym_name_len, names, spin, bool_str[spin]))
    output_lines.sort()
    print "\n".join(output_lines), "\n"
