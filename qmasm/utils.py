###################################
# QMASM utility functions         #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from collections import defaultdict
import math
import qmasm
import sys

class RemainingNextException(Exception):
    'This exception is thrown if a "!next." directive can\'t be replaced.'
    pass

def symbol_to_number(sym, prefix=None, next_prefix=None):
    "Map from a symbol to a number, creating a new association if necessary."
    global sym2num, next_sym_num

    # Replace "!next." by substituting prefixes in the name.
    if "!next." in sym:
        if prefix == None or next_prefix == None:
            raise RemainingNextException
        sym = sym.replace(prefix + "!next.", next_prefix)

    # Return the symbol's logical qubit number or allocate a new one.
    try:
        return qmasm.sym2num[sym]
    except KeyError:
        qmasm.next_sym_num += 1
        qmasm.sym2num[sym] = qmasm.next_sym_num
        return qmasm.next_sym_num

def abend(str):
    "Abort the program on an error."
    sys.stderr.write("%s: %s\n" % (qmasm.progname, str))
    sys.exit(1)

def dict_to_list(d):
    "Convert a dictionary to a list."
    if len(d) == 0:
        return []
    llen = max(d.keys()) + 1
    lst = [0] * llen
    for k, v in d.items():
        lst[k] = v
    return lst

def list_to_dict(lst):
    "Convert a list to a dictionary."
    return defaultdict(lambda: 0.0,
                       {k: lst[k]
                        for k in range(len(lst))
                        if lst[k] != 0.0})

def chimera_topology(solver):
    "Return the topology of the chimera graph associated with a given solver."
    try:
        nominal_qubits = solver.properties["num_qubits"]
    except KeyError:
        # The Ising heuristic solver is an example of a solver that lacks a
        # fixed hardware representation.  We therefore claim an arbitrary
        # topology (the D-Wave 2X's 12x12x4x2 topology).
        return 4, 12, 12
    couplers = solver.properties["couplers"]
    deltas = [abs(c1 - c2) for c1, c2 in couplers]
    delta_tallies = {d: 0 for d in deltas}
    for d in deltas:
        delta_tallies[d] += 1
        sorted_tallies = sorted(delta_tallies.items(), key=lambda dt: dt[1], reverse=True)
    L = sorted_tallies[0][0]
    M = 1
    for d, t in sorted_tallies[1:]:
        if d > 2*L:
            M = d // (2*L)
            break
    N = (nominal_qubits + 2*L*M - 1) // (2*L*M)
    return L, M, N

def edges_to_neighbor_list(pairs):
    "Return a mapping from each node to every node it touches."
    nodes = {}
    for n1, n2 in pairs:
        try:
            nodes[n1].add(n2)
        except KeyError:
            nodes[n1] = set([n2])
        try:
            nodes[n2].add(n1)
        except KeyError:
            nodes[n2] = set([n1])
    return nodes

# Thanks to Carleton Coffrin for proposing the embeddability metric
# used in the following function.
def maybe_embeddable(edges, adj):
    """Return (among other data) False if a QUBO cannot be embedded, True
    if there's a chance it can be embedded."""
    # Determine what we have and what we need in terms of graph connectivity.
    graph_needed = edges_to_neighbor_list(edges)
    num_edges_needed = len(edges)
    num_nodes_needed = len(graph_needed)
    graph_avail = edges_to_neighbor_list(adj)
    num_edges_avail = len(adj)
    num_nodes_avail = len(graph_avail)
    max_degree_avail = max([len(peers) for peers in graph_avail.values()])

    # Compute the minimum number of extra nodes/edges needed to map
    # the graph we need to the graph we have available.
    extras = 0
    for peers in graph_needed.values():
        deg = len(peers)
        if deg > max_degree_avail:
            extras += math.ceil(float(deg)/float(max_degree_avail)) - 1

    # If we need more nodes/edges that are available, the graph can't
    # be embedded.
    embed = True
    if num_nodes_needed + extras > num_nodes_avail:
        embed = False
    if num_edges_needed + extras > num_edges_avail:
        embed = False

    # Compute a histogram of node degrees.
    hist = {}
    for peers in graph_needed.values():
        deg = len(peers)
        try:
            hist[deg] += 1
        except KeyError:
            hist[deg] = 1

    # Return a set of useful information.
    return embed, extras, (num_nodes_needed, num_nodes_avail), (num_edges_needed, num_edges_avail), sorted(hist.items())
