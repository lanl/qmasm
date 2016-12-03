###################################
# Output QUBOs in various formats #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import qmasm
import re
import sys
from dwave_sapi2.util import ising_to_qubo, linear_index_to_chimera

def open_output_file(oname):
    "Open a file or standard output."
    if oname == "<stdout>":
        outfile = sys.stdout
    else:
        try:
            outfile = open(oname, "w")
        except IOError:
            qmasm.abend('Failed to open %s for output' % oname)
    return outfile

def coupler_number(M, N, L, q1, q2):
    "Map a pair of qubits to a coupler number a la the dw command."
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

def output_qubist(outfile, as_qubo, problem):
    "Output weights and strengths in Qubist format, either Ising or QUBO."
    if as_qubo and not problem.qubo:
        qprob = problem.convert_to_qubo()
        output_weights, output_strengths = qprob.weights, qprob.strengths
    elif not as_qubo and problem.qubo:
        iprob = problem.convert_to_ising()
        output_weights, output_strengths = iprob.weights, iprob.strengths
    else:
        output_weights = problem.weights
        output_strengths = problem.strengths
    data = []
    for q, wt in sorted(output_weights.items()):
        if wt != 0.0:
            data.append("%d %d %.10g" % (q, q, wt))
    for sp, str in sorted(output_strengths.items()):
        if str != 0.0:
            data.append("%d %d %.10g" % (sp[0], sp[1], str))

    # Output the header and data in Qubist format.
    try:
        num_qubits = qmasm.solver.properties["num_qubits"]
    except KeyError:
        # The Ising heuristic solver is an example of a solver that lacks a
        # fixed hardware representation.  We therefore assert that the number
        # of qubits is exactly the number of qubits we require.
        num_qubits = len(output_weights)
    outfile.write("%d %d\n" % (num_qubits, len(data)))
    for d in data:
        outfile.write("%s\n" % d)

def output_dw(outfile, problem):
    "Output weights and strengths in dw format."
    if not problem.qubo:
        qprob = problem.convert_to_qubo()
        output_weights, output_strengths = qprob.weights, qprob.strengths
    else:
        output_weights = problem.weights
        output_strengths = problem.strengths
    try:
        L, M, N = qmasm.chimera_topology(qmasm.solver)
    except KeyError:
        qmasm.abend("Failed to query the chimera topology")
    wdata = []
    for q in range(len(output_weights)):
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

def output_qbsolv(outfile, problem):
    "Output weights and strengths in qbsolv format."
    key_width = 0
    val_width = 0
    items = []
    for s, n in qmasm.sym2num.items():
        if "$" in s:
            continue
        if len(s) > key_width:
            key_width = len(s)
        nstr = str(n)
        if len(nstr) > val_width:
            val_width = len(nstr)
        items.append((s, nstr))
    items.sort()
    for s, nstr in items:
        outfile.write("c %-*s --> %-*s\n" % (key_width, s, val_width, nstr))
    if not problem.qubo:
        qprob = problem.convert_to_qubo()
        output_weights, output_strengths = qprob.weights, qprob.strengths
    else:
        output_weights = problem.weights
        output_strengths = problem.strengths
    num_output_weights = len(output_weights)
    for q1, q2 in output_strengths.keys():
        # A large-numbered qubit might have zero weight but nonzero strength.
        num_output_weights = max(num_output_weights, q1 + 1, q2 + 1)
    nonzero_strengths = [s for s in output_strengths.values() if s != 0.0]
    outfile.write("p qubo 0 %d %d %d\n" % (num_output_weights, num_output_weights, len(nonzero_strengths)))
    for q in range(num_output_weights):
        outfile.write("%d %d %.10g\n" % (q, q, output_weights[q]))
    for qs in sorted(output_strengths.keys()):
        s = output_strengths[qs]
        if s != 0.0:
            outfile.write("%d %d %.10g\n" % (qs[0], qs[1], s))

def output_qmasm(outfile):
    "Output weights and strengths as a flattened QMASM source file."
    for p in qmasm.program:
        outfile.write("%s\n" % p.as_str())

# quote was adapted from Python 3's shlex module because the quote method isn't
# included in Python 2's shlex.
_find_unsafe = re.compile(r'[^\w@%+=:,./-]').search
def quote(s):
    """Return a shell-escaped version of the string *s*."""
    if not s:
        return "''"
    if _find_unsafe(s) is None:
        return s

    # Use single quotes, and put single quotes into double quotes
    # the string $'b is then quoted as '$'"'"'b'.
    return "'" + s.replace("'", "'\"'\"'") + "'"

def output_minizinc(outfile, problem):
    "Output weights and strengths as a MiniZinc constraint problem."
    # Write some header information.
    outfile.write("""% Use MiniZinc to minimize a given Hamiltonian.
%
% Producer:     QMASM (https://github.com/losalamos/qmasm/)
% Author:       Scott Pakin (pakin@lanl.gov)
""")
    outfile.write("%% Command line: %s\n\n" % " ".join([quote(a) for a in sys.argv]))

    # Output all QMASM variables as MiniZinc variables.
    all_weights = set(problem.weights.keys())
    all_weights.update([qs[0] for qs in problem.strengths.keys()])
    all_weights.update([qs[1] for qs in problem.strengths.keys()])
    for q in sorted(all_weights):
        outfile.write("var bool: q%d;\n" % q)
    outfile.write("\n")

    # We'll dynamically convert Booleans to spins.
    if not problem.qubo:
        outfile.write("function var int: to_spin(var bool: q) = 2*q - 1;\n\n")

    # Express energy as one, big Hamiltonian.
    scale_to_int = lambda f: int(round(10000.0*f))
    outfile.write("var int: energy =\n")
    if problem.qubo:
        weight_terms = ["%8d * q%d" % (scale_to_int(wt), q) for q, wt in sorted(problem.weights.items())]
        strength_terms = ["%8d * q%d * q%d" % (scale_to_int(s), qs[0], qs[1]) for qs, s in sorted(problem.strengths.items())]
    else:
        weight_terms = ["%8d * to_spin(q%d)" % (scale_to_int(wt), q) for q, wt in sorted(problem.weights.items())]
        strength_terms = ["%8d * to_spin(q%d) * to_spin(q%d)" % (scale_to_int(s), qs[0], qs[1]) for qs, s in sorted(problem.strengths.items())]
    all_terms = weight_terms + strength_terms
    outfile.write("  %s;\n" % " +\n  ".join(all_terms))

    # Because we can't both minimize and enumerate all solutions, we do only
    # the former with instructions for the user on how to switch to the latter.
    outfile.write("""
% First pass: Compute the minimum energy.
solve minimize energy;

% Second pass: Find all minimum-energy solutions.
%
% Once you've solved for minimum energy, comment out the "solve minimize
% energy" line, plug the minimal energy value into the following line,
% uncomment it and the "solve satisfy" line, and re-run MiniZinc, requesting
% all solutions this time.
%constraint energy = -12345;
%solve satisfy;

""")

    # Map each logical qubit to one or more symbols.
    num2syms = [[] for _ in range(len(qmasm.sym2num))]
    max_sym_name_len = 7
    for s, n in qmasm.sym2num.items():
        num2syms[n].append(s)
        max_sym_name_len = max(max_sym_name_len, len(repr(num2syms[n])) - 1)
    def compare_syms(a, b):
        "Compare with internal variables appearing after external variables."
        if "$" in a and "$" in b:
            return cmp(a, b)   # Both contain "$"
        elif "$" in a:
            return +1          # Only a contains "$"
        elif "$" in b:
            return -1          # Only b contains "$"
        else:
            return cmp(a, b)   # Neither contains "$"
    for n in range(len(num2syms)):
        num2syms[n].sort(cmp=compare_syms)

    # Output code to show the results symbolically.
    outfile.write("output [\n")
    outfile.write('  "Energy = ", show(energy), "\\n\\n",\n')
    outfile.write('  "%-*s  Spin  Boolean\\n",\n' % (max_sym_name_len, "Name(s)"))
    outfile.write('  "%s  ----  -------\\n",\n' % ("-" * max_sym_name_len))
    outlist = []
    for n in range(len(num2syms)):
        try:
            phys = problem.embedding[n][0]
        except IndexError:
            continue
        syms = " ".join(num2syms[n])
        line = ""
        line += '"%-*s  ", ' % (max_sym_name_len, syms)
        if problem.qubo:
            line += 'show_int(4, q%d), ' % phys
        else:
            line += 'show_int(4, 2*q%d - 1), ' % phys
        line += '"  ", show(if q%d then "True" else "False" endif), ' % phys
        line += '"\\n"'
        outlist.append(line)
    outlist.sort()
    outfile.write("  %s\n];\n" % ",\n  ".join(outlist))

def write_output(problem, oname, oformat, as_qubo):
    "Write an output file in one of a variety of formats."

    # Open the output file.
    outfile = open_output_file(oname)

    # Output the weights and strengths in the specified format.
    if oformat == "qubist":
        output_qubist(outfile, as_qubo, problem)
    elif oformat == "dw":
        output_dw(outfile, problem)
    elif oformat == "qbsolv":
        output_qbsolv(outfile, problem)
    elif oformat == "qmasm":
        output_qmasm(outfile)
    elif oformat == "minizinc":
        output_minizinc(outfile, problem)

    # Close the output file.
    if oname != "<stdout>":
        outfile.close()
