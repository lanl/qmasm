#####################################
# Update QASM weights and strengths #
# By Scott Pakin <pakin@lanl.gov>   #
#####################################

from qasm.globals import *
from qasm.utils import *
from dwave_sapi2.util import qubo_to_ising

def convert_to_ising():
    "Canonicalize a QUBO problem into an Ising problem."
    global weights, strengths
    qmatrix = {(q, q): w for q, w in weights.items()}
    qmatrix.update(strengths)
    hvals, strengths, _ = qubo_to_ising(qmatrix)
    strengths = defaultdict(lambda: 0.0, strengths)
    weights.update({i: hvals[i] for i in range(len(hvals))})

def assign_chain_strength(ch_str, qubo):
    """Define a strength for each user-specified chain, and assign strengths
    to those chains."""
    global chain_strength, strengths
    chain_strength = ch_str
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
    elif qubo:
        # With QUBO input we need to divide the chain strength by 4 for
        # consistency with the other coupler strengths.
        chain_strength /= 4.0
    for c in chains.keys():
        strengths[c] += chain_strength
