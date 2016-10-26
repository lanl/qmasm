###################################
# Define an Ising or QUBO problem #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from collections import defaultdict
from dwave_sapi2.util import ising_to_qubo, qubo_to_ising
import copy
import qmasm
import random
import string

def new_internal_sym():
    "Create a new internal symbol."
    while True:
        sym = "$"
        for i in range(5):
            sym += random.choice(string.lowercase)
        if not qmasm.sym2num.has_key(sym):
            return sym

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
            # Pin strength defaults to half the chain strength.
            pin_strength = chain_str/2.0
        elif self.qubo:
            # With QUBO input we need to divide the chain strength by 4 for
            # consistency with the other coupler strengths.
            pin_strength /= 4.0
        return pin_strength

    def pin_qubits(self, pin_str):
        "Use a helper qubit to help pin values to true or false."
        for q_user, b in self.pinned:
            q_helper = qmasm.symbol_to_number(new_internal_sym())
            q1, q2 = q_helper, q_user
            if q1 > q2:
                q1, q2 = q2, q1
            if b:
                self.weights[q_helper] += pin_str/2.0
                self.weights[q_user] += pin_str
                self.strengths[(q1, q2)] += -pin_str/2.0
            else:
                self.weights[q_helper] += pin_str/2.0
                self.weights[q_user] += -pin_str
                self.strengths[(q1, q2)] += pin_str/2.0

    def convert_to_ising(self):
        """Transform a QUBO problem into an Ising problem.  Return the new
        Ising problem."""
        if not self.qubo:
            raise TypeError("Can convert only QUBO problems to Ising problems")
        new_obj = copy.deepcopy(self)
        qmatrix = {(q, q): w for q, w in new_obj.weights.items()}
        qmatrix.update(new_obj.strengths)
        hvals, new_obj.strengths, _ = qubo_to_ising(qmatrix)
        new_obj.strengths = defaultdict(lambda: 0.0, new_obj.strengths)
        new_obj.weights.update({i: hvals[i] for i in range(len(hvals))})
        new_obj.qubo = False
        return new_obj

    def convert_to_qubo(self):
        """Transform an Ising problem into a QUBO problem.  Return the new
        QUBO problem."""
        if self.qubo:
            raise TypeError("Can convert only Ising problems to QUBO problems")
        new_obj = copy.deepcopy(self)
        qmatrix, _ = ising_to_qubo(self.weights, self.strengths)
        new_obj.weights = [0] * len(self.weights)
        for (q1, q2), wt in qmatrix.items():
            if q1 == q2:
                new_obj.weights[q1] = wt
        new_obj.strengths = {(q1, q2): wt for (q1, q2), wt in qmatrix.items() if q1 != q2}
        return new_obj

    def convert_chains_to_aliases(self):
        "Convert chains to aliases where possible."
        # Identify all chains that can be converted to aliases.  A chain is
        # convertible if the qubits on either end have the same point weight
        # applied to them.
        num2allsyms = [[] for _ in range(len(qmasm.sym2num))]
        for s, n in qmasm.sym2num.items():
            num2allsyms[n].append(s)
        make_aliases = []
        for q1, q2 in self.chains:
            if self.weights[q1] == self.weights[q2]:
                make_aliases.append((q1, q2))
        make_aliases.sort(reverse=True, key=lambda qs: (qs[1], qs[0]))

        # Replace each chain in make_aliases with an alias.  Work in reverse
        # order of qubit number and shift all greater qubit numbers downward.
        for q1, q2 in make_aliases:
            # Map q2's symbolic names to q1's.  Shift everything above q2
            # downwards.
            alias_sym2num = {}
            for s, sq in qmasm.sym2num.items():
                if sq == q2:
                    sq = q1
                elif sq > q2:
                    sq -= 1
                alias_sym2num[s] = sq
            qmasm.sym2num = alias_sym2num

            # Elide q2 from the list of weights.
            alias_weights = defaultdict(lambda: 0.0)
            for wq, wt in self.weights.items():
                if wq == q2:
                    continue
                if wq > q2:
                    wq -= 1
                if wt != 0.0:
                    alias_weights[wq] = wt
            alias_weights[q1] += self.weights[q2]   # Conserve overall energy.
            self.weights = alias_weights

            # Replace q2 with q1 in all strengths.  Shift everything above q2
            # downwards.
            alias_strengths = defaultdict(lambda: 0.0)
            for (sq1, sq2), wt in self.strengths.items():
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
            self.strengths = alias_strengths

            # Replace q2 with q1 in all strengths.  Shift everything above q2
            # downwards.
            alias_chains = {}
            for cq1, cq2 in self.chains.keys():
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
            self.chains = alias_chains

            # Replace q2 with q1 in all pinned qubits.  Shift everything above
            # q2 downwards.
            alias_pinned = []
            for pq, b in self.pinned:
                if pq == q2:
                    pq = q1
                if pq > q2:
                    pq -= 1
                alias_pinned.append((pq, b))
            self.pinned = alias_pinned

            # We now have one fewer symbol.
            qmasm.next_sym_num -= 1

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
