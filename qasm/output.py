###################################
# Output QUBOs in various formats #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import qasm
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
            abend('Failed to open %s for output' % oname)
    return outfile

def convert_ising_to_qubo(new_weights, combined_strengths):
    "Convert weights and strengths from Ising to QUBO mode."
    qmatrix, _ = ising_to_qubo(new_weights, combined_strengths)
    output_weights = [0] * len(new_weights)
    for (q1, q2), wt in qmatrix.items():
        if q1 == q2:
            output_weights[q1] = wt
    output_strengths = {(q1, q2): wt for (q1, q2), wt in qmatrix.items() if q1 != q2}
    return output_weights, output_strengths

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

def output_qubist(outfile, as_qubo, new_weights, combined_strengths):
    "Output weights and strengths in Qubist format."
    if as_qubo:
        output_weights, output_strengths = convert_ising_to_qubo(new_weights, combined_strengths)
    else:
        output_weights = new_weights
        output_strengths = combined_strengths
    data = []
    for q in range(len(new_weights)):
        if output_weights[q] != 0.0:
            data.append("%d %d %.10g" % (q, q, output_weights[q]))
    for sp, str in output_strengths.items():
        if str != 0.0:
            data.append("%d %d %.10g" % (sp[0], sp[1], str))

    # Output the header and data in Qubist format.
    try:
        num_qubits = qasm.solver.properties["num_qubits"]
    except KeyError:
        # The Ising heuristic solver is an example of a solver that lacks a
        # fixed hardware representation.  We therefore assert that the number
        # of qubits is exactly the number of qubits we require.
        num_qubits = len(new_weights)
    outfile.write("%d %d\n" % (num_qubits, len(data)))
    for d in data:
        outfile.write("%s\n" % d)

def output_dw(outfile, new_weights, combined_strengths):
    # Output weights and strengths in dw format.
    try:
        L, M, N = qasm.chimera_topology(qasm.solver)
    except KeyError:
        abend("Failed to query the chimera topology")
    output_weights, output_strengths = convert_ising_to_qubo(new_weights, combined_strengths)
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

def output_qbsolv(outfile, new_weights, combined_strengths):
    # Output weights and strengths in qbsolv format.
    output_weights, output_strengths = convert_ising_to_qubo(new_weights, combined_strengths)
    nonzero_strengths = [s for s in output_strengths.values() if s != 0.0]
    outfile.write("p qubo 0 %d %d %d\n" % (len(output_weights), len(output_weights), len(nonzero_strengths)))
    for q in range(len(output_weights)):
        outfile.write("%d %d %.10g\n" % (q, q, output_weights[q]))
    for qs in sorted(output_strengths.keys()):
        s = output_strengths[qs]
        if s != 0.0:
            outfile.write("%d %d %.10g\n" % (qs[0], qs[1], s))

def write_output(problem, oname, oformat, as_qubo):
    "Write an output file in one of a variety of formats."

    # Open the output file.
    outfile = open_output_file(oname)

    # Convert from Ising back to QUBO if --qubo was specified or if we're
    # outputting in dw format, which expects QUBO coefficients.
    if as_qubo or oformat == "dw":
        output_weights, output_strengths = convert_ising_to_qubo(problem.weights, problem.strengths)
    else:
        output_weights = problem.weights
        output_strengths = problem.strengths

    # Output the weights and strengths in the specified format.
    if oformat == "qubist":
        output_qubist(outfile, as_qubo, problem.weights, problem.strengths)
    elif oformat == "dw":
        output_dw(outfile, problem.weights, problem.strengths)
    elif oformat == "qbsolv":
        output_qbsolv(outfile, problem.weights, problem.strengths)

    # Close the output file.
    if oname != "<stdout>":
        outfile.close()
