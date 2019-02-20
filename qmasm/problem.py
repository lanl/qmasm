###################################
# Define an Ising or QUBO problem #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from collections import defaultdict
try:
    from dwave_sapi2.util import ising_to_qubo, qubo_to_ising
except ImportError:
    from .fake_dwave import *
import copy
import qmasm
import random
import string
import sys
try:
    import pulp
    have_pulp = True
except ImportError:
    have_pulp = False

def new_internal_sym():
    "Create a new internal symbol."
    while True:
        sym = "$"
        for i in range(5):
            sym += random.choice(string.ascii_lowercase)
        try:
            n = qmasm.sym_map.to_number(sym)
        except:
            # Symbol does not yet exist.
            return sym

class DisjointSet(object):
    "Set in a disjoint-set forest"

    def __init__(self, val=None):
        self.parent = self
        self.rank = 0
        self.contents = val

    def find(self):
        "Return an arbitrary element in the same set as us."
        if self.parent == self:
            return self
        else:
            self.parent = self.parent.find()
            return self.parent

    def union(self, other):
        "Destructively join two sets."
        self_root = self.find()
        other_root = other.find()
        if self_root.rank > other_root.rank:
            other_root.parent = self_root
        elif self_root.rank < other_root.rank:
            self_root.parent = other_root
        elif self_root != other_root:
            other_root.parent = self_root
            self_root.rank = self_root.rank + 1

class Problem(object):
    "Represent either an Ising or QUBO problem."

    def __init__(self, qubo):
        self.qubo = qubo     # True=QUBO; False=Ising
        self.weights = defaultdict(lambda: 0.0)    # Map from a spin to a point weight
        self.strengths = defaultdict(lambda: 0.0)  # Map from a pair of spins to a coupler strength
        self.chains = set()     # Subset of strengths keys that represents user-defined chains (always logical)
        self.antichains = set() # Subset of strengths keys that represents user-defined anti-chains (always logical)
        self.pin_chains = set() # Subset of strengths keys that represents {helper, pinned variable} pairs (always logical)
        self.pinned = []     # Pairs of {unique number, Boolean} to pin
        self.offset = 0.0    # Value to add to QUBO energy to convert to Ising energy or vice versa
        self.known_values = {}    # Map from symbol name to spin for values known a priori
        self.simple_offset = 0.0  # Value to add to Ising energy to compensate for problem simplification
        self.assertions = []      # List of assertions (as ASTs) to enforce

    def assign_chain_strength(self, ch_str):
        """Define a strength for each user-specified and automatically
        generated chain, and assign strengths to those chains (and negative
        strength to all anti-chains).  Return the computed chain strength."""
        chain_strength = ch_str
        if chain_strength == None:
            # Chain strength defaults to twice the maximum strength in the data.
            try:
                chain_strength = -2*max([abs(w) for w in self.strengths.values()])
            except ValueError:
                # No strengths -- use weights instead.
                try:
                    chain_strength = -2*max([abs(w) for w in self.weights.values()])
                except ValueError:
                    # No weights or strengths -- arbitrarily choose -1.
                    chain_strength = -1.0
        elif self.qubo:
            # With QUBO input we need to divide the chain strength by 4 for
            # consistency with the other coupler strengths.
            chain_strength /= 4.0
        for c in self.chains:
            self.strengths[c] += chain_strength
        for c in self.antichains:
            self.strengths[c] -= chain_strength
        return chain_strength

    def assign_pin_weight(self, pin_str, chain_str):
        """Define a strength for each explicitly pinned qubit.  Return the
        computed pin strength."""
        pin_weight = pin_str
        if pin_weight == None:
            # Pin strength defaults to the chain strength.
            pin_weight = chain_str
        elif self.qubo:
            # With QUBO input we need to divide the chain strength by 4 for
            # consistency with the other coupler strengths.
            pin_weight /= 4.0
        return pin_weight

    def pin_qubits(self, pin_str, chain_str):
        "Use a helper qubit to help pin values to true or false."
        for q_user, b in self.pinned:
            q_helper = qmasm.symbol_to_number(new_internal_sym())
            q1, q2 = q_helper, q_user
            if q1 > q2:
                q1, q2 = q2, q1
            if q_helper not in self.weights:
                self.weights[q_helper] = 0
            if b:
                self.weights[q_helper] -= pin_str
            else:
                self.weights[q_helper] += pin_str
            self.strengths[(q1, q2)] += -chain_str
            self.antichains.add((q1, q2))
            self.pin_chains.add((q_helper, q_user))

    def convert_to_ising(self):
        """Transform a QUBO problem into an Ising problem.  Return the new
        Ising problem."""
        if not self.qubo:
            raise TypeError("Can convert only QUBO problems to Ising problems")
        new_obj = copy.deepcopy(self)
        qmatrix = {(q, q): w for q, w in new_obj.weights.items()}
        qmatrix.update(new_obj.strengths)
        hvals, new_obj.strengths, qoffset = qubo_to_ising(qmatrix)
        new_obj.strengths = qmasm.canonicalize_strengths(new_obj.strengths)
        new_obj.weights = {i: hvals[i] for i in range(len(hvals))}
        new_obj.offset = qoffset
        new_obj.qubo = False
        return new_obj

    def convert_to_qubo(self):
        """Transform an Ising problem into a QUBO problem.  Return the new
        QUBO problem."""
        if self.qubo:
            raise TypeError("Can convert only Ising problems to QUBO problems")
        new_obj = copy.deepcopy(self)
        qmatrix, qoffset = ising_to_qubo(qmasm.dict_to_list(self.weights), self.strengths)
        new_obj.offset = qoffset
        new_obj.weights = defaultdict(lambda: 0.0,
                                      {q1: wt
                                       for (q1, q2), wt in qmatrix.items()
                                       if q1 == q2})
        new_obj.strengths = qmasm.canonicalize_strengths({(q1, q2): wt
                                                          for (q1, q2), wt in qmatrix.items()
                                                          if q1 != q2})
        new_obj.qubo = True
        return new_obj

    def convert_chains_to_aliases(self):
        "Replace user-specified chains with aliases."

        # Group qubits that can be aliased.
        num2alias = {}  # Map from a qubit number to a disjoint set (which maps to a qubit number)
        for q1, q2 in self.chains:
            if q1 not in num2alias:
                num2alias[q1] = DisjointSet(q1)
            if q2 not in num2alias:
                num2alias[q2] = DisjointSet(q2)
            num2alias[q1].union(num2alias[q2])

        # Regenerate our chains, discarding any that have been merged into a
        # single qubit.
        new_chains = set()
        for q1, q2 in self.chains:
            try:
                new_q1 = num2alias[q1].find().contents
            except KeyError:
                new_q1 = q1
            try:
                new_q2 = num2alias[q2].find().contents
            except KeyError:
                new_q2 = q2
            if new_q1 == new_q2:
                continue
            if new_q1 > new_q2:
                new_q1, new_q2 = new_q2, new_q1
            new_chains.add((new_q1, new_q2))
        self.chains = new_chains

        # Regenerate our anti-chains, discarding any that have been merged into
        # a single qubit.
        new_antichains = set()
        for q1, q2 in self.antichains:
            try:
                new_q1 = num2alias[q1].find().contents
            except KeyError:
                new_q1 = q1
            try:
                new_q2 = num2alias[q2].find().contents
            except KeyError:
                new_q2 = q2
            if new_q1 == new_q2:
                continue
            if new_q1 > new_q2:
                new_q1, new_q2 = new_q2, new_q1
            new_antichains.add((new_q1, new_q2))
        self.antichains = new_antichains

        # Regenerate our weights.
        new_weights = defaultdict(lambda: 0.0)
        for q, wt in self.weights.items():
            try:
                new_q = num2alias[q].find().contents
            except KeyError:
                new_q = q
            new_weights[new_q] += wt
        self.weights = new_weights

        # Regenerate our strengths.
        new_strengths = defaultdict(lambda: 0.0)
        for (q1, q2), wt in self.strengths.items():
            try:
                new_q1 = num2alias[q1].find().contents
            except KeyError:
                new_q1 = q1
            try:
                new_q2 = num2alias[q2].find().contents
            except KeyError:
                new_q2 = q2
            if new_q1 == new_q2:
                continue
            if new_q1 > new_q2:
                new_q1, new_q2 = new_q2, new_q1
            new_strengths[(new_q1, new_q2)] += wt
        self.strengths = new_strengths

        # Regenerate our pinned values.
        new_pinned = {}
        for q, b in self.pinned:
            try:
                new_q = num2alias[q].find().contents
            except KeyError:
                new_q = q
            new_pinned[new_q] = b
        self.pinned = sorted(new_pinned.items())

        # Regenerate the global symbol table.
        new_sym2num = {}
        for s, q in qmasm.sym_map.symbol_number_items():
            try:
                new_q = num2alias[q].find().contents
            except KeyError:
                new_q = q
            new_sym2num[s] = new_q
        qmasm.sym_map.overwrite_with(new_sym2num)

        # Renumber all of the above to compact the qubit numbers.
        qubits_used = set([q.find().contents for q in num2alias.values()])
        qubits_used.update(self.weights.keys())
        for q1, q2 in self.strengths.keys():
            qubits_used.add(q1)
            qubits_used.add(q2)
        qmap = dict(zip(sorted(qubits_used), range(len(qubits_used))))
        self.chains = set([(qmap[q1], qmap[q2]) for q1, q2 in self.chains])
        self.antichains = set([(qmap[q1], qmap[q2]) for q1, q2 in self.antichains])
        self.weights = defaultdict(lambda: 0.0,
                                   {qmap[q]: wt for q, wt in self.weights.items()})
        self.strengths = qmasm.canonicalize_strengths({(qmap[q1], qmap[q2]): wt for (q1, q2), wt in self.strengths.items()})
        self.pinned = [(qmap[q], b) for q, b in self.pinned]
        qmasm.sym_map.overwrite_with({s: qmap[q] for s, q in qmasm.sym_map.symbol_number_items()})

    def find_disconnected_variables(self):
        """Return a list of variables that are named but not coupled to any
        other variable."""
        # Construct a set of valid qubit numbers.
        valid_nums = set()
        for (a, b), str in self.strengths.items():
            if str != 0.0:
                valid_nums.add(a)
                valid_nums.add(b)

        # Complain about any variable whose number is not in the valid set.
        invalid_syms = set()
        for num in qmasm.sym_map.all_numbers():
            if num not in valid_nums:
                invalid_syms.update(qmasm.sym_map.to_symbols(num))
        return invalid_syms

    def couplers_to_strengths(self, couplers):
        """Convert a set of logical couplers to a set of physical strengths.
        Ignore couplers involving elided variables."""
        strengths = set()
        for lqs in couplers:
            try:
                pq1_list = self.embedding[lqs[0]]
                pq2_list = self.embedding[lqs[1]]
            except IndexError:
                continue   # One of the logical qubits was elided.
            for pq1 in pq1_list:
                for pq2 in pq2_list:
                    if pq1 > pq2:
                        pq1, pq2 = pq2, pq1
                    if (pq1, pq2) in self.strengths:
                        strengths.add((pq1, pq2))
        return strengths

    def estimate_energy(self):
        "Estimate minimum energy by constructing and solving a majorization problem."
        # Do nothing if the pulp module isn't available.
        if not have_pulp:
            return None

        # Convert from Ising to QUBO.
        if self.qubo:
            qubo = self
        else:
            qubo = self.convert_to_qubo()

        # Allocate all of the variables we'll need, one x per linear term and
        # one y per quadratic term.
        all_vars = set(qubo.weights.keys())
        for (q0, q1) in qubo.strengths.keys():
            all_vars.add(q0)
            all_vars.add(q1)
        x = {q: pulp.LpVariable("x_%d" % q, cat="Continuous", lowBound=0.0, upBound=1.0) for q in all_vars}
        y = {(q0, q1): pulp.LpVariable("y_%d_%d" % (q0, q1), cat="Continuous", lowBound=0.0, upBound=1.0)
             for q0 in all_vars
             for q1 in all_vars
             if q0 < q1}

        # Specify the objective function.
        prob = pulp.LpProblem("major", pulp.LpMinimize)
        lterms = pulp.LpAffineExpression([(xi, qubo.weights[i])
                                          for i, xi in x.items()
                                          if qubo.weights[i] != 0.0])
        qterms = pulp.LpAffineExpression([(yij, qubo.strengths[(i, j)])
                                          for (i, j), yij in y.items()
                                          if qubo.strengths[(i, j)] != 0.0])
        prob += qubo.offset + lterms + qterms

        # Specify constraints on the (linearized) quadratic terms.
        for (i, j), cij in qubo.strengths.items():
            yij = y[(i, j)]
            if cij > 0.0:
                prob += yij >= x[i] + x[j] - 1.0
                prob += yij >= 0.0
            elif cij < 0.0:
                prob += yij <= x[i]
                prob += yij <= x[j]

        # Solve the linear program.
        if prob.solve() != pulp.LpStatusOptimal:
            return None
        return pulp.value(prob.objective)

    def scale_weights_strengths(self, verbosity):
        "Manually scale the weights and strengths so Qubist doesn't complain."
        h_range = self.h_range
        j_range = self.j_range
        weight_list = qmasm.dict_to_list(self.weights)
        old_cap = max([abs(w) for w in weight_list + list(self.strengths.values())])
        new_cap = min(-h_range[0], h_range[1], -j_range[0], j_range[1])
        if old_cap == 0.0:
            # Handle the obscure case of a zero old_cap.
            old_cap = new_cap
        self.range_scale = new_cap/old_cap
        self.weights = qmasm.list_to_dict([w*self.range_scale for w in weight_list])
        self.strengths = {js: w*self.range_scale for js, w in self.strengths.items()}
        self.offset *= self.range_scale
        if verbosity >= 1 and old_cap != new_cap:
            sys.stderr.write("Scaling weights and strengths from [%.10g, %.10g] to [%.10g, %.10g].\n\n" % (-old_cap, old_cap, -new_cap, new_cap))

    def update_strengths_from_chains(self):
        "Update strengths using the chains introduced by embedding."
        self.strengths.update({c: qmasm.chain_strength for c in self.embedder_chains})

    def find_stray_variables(self):
        "Return a list of variables that are not coupled to any other variable."
        # Tally the number of times each variable is coupled.
        coupled_vars = defaultdict(lambda: 0.0)
        for q0, q1 in self.strengths:
            coupled_vars[q0] += 1
            coupled_vars[q1] += 1

        # Any variable with no couplings is a stray.
        strays = set()
        for q in self.weights:
            if q not in coupled_vars:
                strays.add(q)

        # Any pinned variable with exactly one coupling (to its helper
        # variable) is a stray.
        for q_helper, q_user in self.pin_chains:
            if coupled_vars[q_user] <= 1:
                strays.add(q_user)
        return strays

    def append_assertions_from_statements(self):
        "Convert user-specified chains, anti-chains, and pins to assertions."
        # Convert certain statement types to assertions.
        # TODO: Quote variables containing special characters.
        ap = qmasm.AssertParser()
        for stmt in qmasm.program:
            if stmt.__class__ == qmasm.AntiChain:
                ast = ap.parse("%s /= %s" % (stmt.sym1, stmt.sym2))
                self.assertions.append(ast)
            elif stmt.__class__ == qmasm.Chain:
                ast = ap.parse("%s = %s" % (stmt.sym1, stmt.sym2))
                self.assertions.append(ast)
            elif stmt.__class__ == qmasm.Pin:
                ast = ap.parse("%s = %d" % (stmt.sym, int(stmt.goal)))
                self.assertions.append(ast)
