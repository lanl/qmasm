###################################
# Define an Ising or QUBO problem #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import dimod
from collections import defaultdict

# Problem is currently just a thin veneer over dimod.BinaryQuadraticModel.  If
# it turns out we don't even need this veneer, we may replace it with direct
# use of dimod.BinaryQuadraticModel.
class Problem(object):
    "Represent either an Ising or QUBO problem."

    def __init__(self, qubo):
        self.qubo = qubo     # True=QUBO; False=Ising
        self.weights = defaultdict(lambda: 0.0)    # Map from a spin to a point weight
        self.strengths = defaultdict(lambda: 0.0)  # Map from a pair of spins to a coupler strength
        self.chains = set()       # Subset of strengths keys that represents user-defined chains (always logical)
        self.antichains = set()   # Subset of strengths keys that represents user-defined anti-chains (always logical)
        self.pin_chains = set()   # Subset of strengths keys that represents {helper, pinned variable} pairs (always logical)
        self.assertions = []      # List of assertions (as ASTs) to enforce
        self.pending_asserts = [] # List of {string, op, string} tuples pending conversion to assertions
        self.pinned = []          # Pairs of {unique number, Boolean} to pin

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

    def as_bqm(self):
        "Return a BinaryQuadraticModel version of the Problem."
        # Create a BQM.
        btype = dimod.SPIN
        if self.qubo:
            btype = dimod.BINARY
        bqm = dimod.BinaryQuadraticModel(self.weights, self.strengths, 0, btype, problem=self)
        if self.qubo:
            bqm.change_vartype(dimod.SPIN, inplace=True)

        # Pin all variables the user asked to pin.
        bool2spin = {False: -1, True: +1}
        pins = {q: bool2spin[b] for q, b in self.pinned}
        bqm.fix_variables(pins)

        # Return the BQM.
        return bqm
