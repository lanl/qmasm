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

def apply_prefix(sym, prefix=None, next_prefix=None):
    "Apply a prefix to a symbol name, replacing !next. with the next prefix."
    if prefix != None:
        sym = prefix + sym
    if "!next." in sym:
        if prefix == None or next_prefix == None:
            raise RemainingNextException
        sym = sym.replace(prefix + "!next.", next_prefix)
    return sym

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
        return qmasm.sym_map.to_number(sym)
    except KeyError:
        return qmasm.sym_map.new_symbol(sym)

def abend(str):
    "Abort the program on an error."
    sys.stderr.write("%s: %s\n" % (qmasm.progname, str))
    sys.exit(1)

def warn(str):
    "Issue a warning message but continue execution."
    sys.stderr.write("%s: Warning: %s\n" % (qmasm.progname, str))

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

class NonChimera(Exception):
    "Exception thrown when finding the topology of a non-Chimera graph."
    pass

def chimera_topology(solver):
    """Return the topology of the Chimera graph associated with a given solver.
    Throw NonChimera if the topology is not known to be a Chimera graph."""
    try:
        nominal_qubits = solver.properties["num_qubits"]
    except KeyError:
        # The Ising heuristic solver is an example of a solver that lacks a
        # fixed hardware representation.
        raise NonChimera
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

def canonicalize_strengths(strs):
    "Combine edges (A, B) and (B, A) into (A, B) with A < B."
    new_strs = defaultdict(lambda: 0.0)
    for (q1, q2), wt in strs.items():
        if q1 == q2:
            continue          # Discard vertex weights.
        if wt == 0.0:
            continue          # Discard zero weights.
        if q1 > q2:
            q1, q2 = q2, q1   # Canonicalize vertex order.
        new_strs[(q1, q2)] += wt
    return new_strs

class SymbolMapping:
    "Map between symbols and numbers."

    def __init__(self):
        self.sym2num = {}
        self.num2syms = {}
        self.next_sym_num = 0

    def new_symbol(self, sym):
        "Assign the next available number to a symbol."
        if sym in self.sym2num:
            raise Exception("Internal error: Symbol %s is already defined" % sym)
        self.sym2num[sym] = self.next_sym_num
        self.num2syms[self.next_sym_num] = set([sym])
        self.next_sym_num += 1
        return self.sym2num[sym]

    def to_number(self, sym):
        "Map a symbol to a single number."
        return self.sym2num[sym]

    def to_symbols(self, num):
        "Map a number to a set of one or more symbols."
        return self.num2syms[num]

    def max_number(self):
        "Return the maximum symbol number used so far."
        return self.next_sym_num - 1

    def all_numbers(self):
        "Return an unordered list of all numbers used."
        return self.num2syms.keys()

    def all_symbols(self):
        "Return an unordered list of all symbols used."
        return self.sym2num.keys()

    def symbol_number_items(self):
        "Return a list of {symbol, number} pairs."
        return self.sym2num.items()

    def alias(self, sym1, sym2):
        "Make two symbols point to the same number."
        # Ensure that both symbols are defined.
        n_new_defs = 0   # Number of new definitions made
        try:
            num1 = self.sym2num[sym1]
        except KeyError:
            num1 = self.new_symbol(sym1)
            n_new_defs += 1
        try:
            num2 = self.sym2num[sym2]
        except KeyError:
            num2 = self.new_symbol(sym2)
            n_new_defs += 1

        # Do nothing if the symbols are already aliased.
        if num1 == num2:
            return num1

        # Abort if both symbols were previously defined (and do not alias each
        # other).  In this case, we'd need to merge the associated weights and
        # strengths, and we don't currently do that.
        if n_new_defs == 0:
            abend("Unable to alias pre-existing variables %s and %s" % (sym1, sym2))

        # Replace all occurrences of the larger number with the smaller in
        # num2syms.
        new_num = min(num1, num2)
        old_num = max(num1, num2)
        self.num2syms[new_num].update(self.num2syms[old_num])
        del self.num2syms[old_num]

        # Regenerate sym2num from num2syms.
        self.sym2num = {}
        for n, ss in self.num2syms.items():
            for s in ss:
                self.sym2num[s] = n
        return new_num

    def overwrite_with(self, sym2num):
        "Overwrite the map's contents with a given symbol-to-number map."
        self.sym2num = sym2num
        self.num2syms = {}
        for s, n in self.sym2num.items():
            try:
                self.num2syms[n].add(s)
            except KeyError:
                self.num2syms[n] = set([s])
        if len(self.num2syms.keys()) == 0:
            self.next_sym_num = 0
        else:
            self.next_sym_num = max(self.num2syms.keys()) + 1
