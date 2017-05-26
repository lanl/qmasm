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

def new_internal_sym():
    "Create a new internal symbol."
    while True:
        sym = "$"
        for i in range(5):
            sym += random.choice(string.ascii_lowercase)
        if sym not in qmasm.sym2num:
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
        self.chains = {}     # Subset of strengths keys that represents chains
        self.pinned = []     # Pairs of {unique number, Boolean} to pin

    def assign_chain_strength(self, ch_str):
        """Define a strength for each user-specified and automatically generated
        chain, and assign strengths to those chains.  Return the computed
        chain strength."""
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
        for c in self.chains.keys():
            self.strengths[c] += chain_strength
        return chain_strength

    def assign_pin_strength(self, pin_str, chain_str):
        """Define a strength for each explicitly pinned qubit.  Return the
        computed pin strength."""
        pin_strength = pin_str
        if pin_strength == None:
            # Pin strength defaults to the chain strength.
            pin_strength = chain_str
        elif self.qubo:
            # With QUBO input we need to divide the chain strength by 4 for
            # consistency with the other coupler strengths.
            pin_strength /= 4.0
        return pin_strength

    def pin_qubits(self, pin_str, chain_str):
        "Use a helper qubit to help pin values to true or false."
        for q_user, b in self.pinned:
            q_helper = qmasm.symbol_to_number(new_internal_sym())
            q1, q2 = q_helper, q_user
            if q1 > q2:
                q1, q2 = q2, q1
            if b:
                self.weights[q_helper] -= pin_str
            else:
                self.weights[q_helper] += pin_str
            self.strengths[(q1, q2)] += -chain_str

    def convert_to_ising(self):
        """Transform a QUBO problem into an Ising problem.  Return the new
        Ising problem."""
        if not self.qubo:
            raise TypeError("Can convert only QUBO problems to Ising problems")
        new_obj = copy.deepcopy(self)
        qmatrix = {(q, q): w for q, w in new_obj.weights.items()}
        qmatrix.update(new_obj.strengths)
        hvals, new_obj.strengths, _ = qubo_to_ising(qmatrix)
        new_obj.strengths = qmasm.canonicalize_strengths(new_obj.strengths)
        new_obj.weights.update({i: hvals[i] for i in range(len(hvals))})
        new_obj.qubo = False
        return new_obj

    def convert_to_qubo(self):
        """Transform an Ising problem into a QUBO problem.  Return the new
        QUBO problem."""
        if self.qubo:
            raise TypeError("Can convert only Ising problems to QUBO problems")
        new_obj = copy.deepcopy(self)
        qmatrix, _ = ising_to_qubo(qmasm.dict_to_list(self.weights), self.strengths)
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
        """Convert chains to aliases where possible.  A chain is
        convertible if the qubits on either end have the same point weight
        applied to them."""

        # Group qubits that can be aliased.
        num2alias = {}  # Map from a qubit number to a disjoint set (which maps to a qubit number)
        for q1, q2 in self.chains:
            if self.weights[q1] == self.weights[q2]:
                if q1 not in num2alias:
                    num2alias[q1] = DisjointSet(q1)
                if q2 not in num2alias:
                    num2alias[q2] = DisjointSet(q2)
                num2alias[q1].union(num2alias[q2])

        # Regenerate our chains, discarding any that have been merged into a
        # single qubit.
        new_chains = {}
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
            new_chains[(new_q1, new_q2)] = None
        self.chains = new_chains

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
        for s, q in qmasm.sym2num.items():
            try:
                new_q = num2alias[q].find().contents
            except KeyError:
                new_q = q
            new_sym2num[s] = new_q
        qmasm.sym2num = new_sym2num

        # Renumber all of the above to compact the qubit numbers.
        qmap = dict(zip(self.weights.keys(), range(len(self.weights))))
        self.chains = {(qmap[q1], qmap[q2]): None for q1, q2 in self.chains.keys()}
        self.weights = defaultdict(lambda: 0.0,
                                   {qmap[q]: wt for q, wt in self.weights.items()})
        self.strengths = qmasm.canonicalize_strengths({(qmap[q1], qmap[q2]): wt for (q1, q2), wt in self.strengths.items()})
        self.pinned = [(qmap[q], b) for q, b in self.pinned]
        qmasm.sym2num = {s: qmap[q] for s, q in qmasm.sym2num.items()}
        qmasm.next_sym_num = max(qmasm.sym2num.values())

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
        for sym, num in qmasm.sym2num.items():
            if num not in valid_nums:
                invalid_syms.add(sym)
        return invalid_syms
