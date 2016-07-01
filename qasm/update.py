#####################################
# Update QASM weights and strengths #
# By Scott Pakin <pakin@lanl.gov>   #
#####################################

import qasm
from dwave_sapi2.util import qubo_to_ising

def convert_to_ising():
    "Canonicalize a QUBO problem into an Ising problem."
    qmatrix = {(q, q): w for q, w in weights.items()}
    qmatrix.update(strengths)
    hvals, qasm.strengths, _ = qubo_to_ising(qmatrix)
    qasm.strengths = defaultdict(lambda: 0.0, qasm.strengths)
    qasm.weights.update({i: hvals[i] for i in range(len(hvals))})

def assign_chain_strength(ch_str, qubo):
    """Define a strength for each user-specified and automatically generated
    chain, and assign strengths to those chains."""
    qasm.chain_strength = ch_str
    if qasm.chain_strength == None:
        # Chain strength defaults to the maximum strength in the data.
        try:
            qasm.chain_strength = -max([abs(w) for w in qasm.strengths.values()])
        except ValueError:
            # No strengths -- use weights instead.
            try:
                qasm.chain_strength = -max([abs(w) for w in qasm.weights.values()])
            except ValueError:
                # No weights or strengths -- arbitrarily choose -1.
                qasm.chain_strength = -1.0
    elif qubo:
        # With QUBO input we need to divide the chain strength by 4 for
        # consistency with the other coupler strengths.
        qasm.chain_strength /= 4.0
    for c in qasm.chains.keys():
        qasm.strengths[c] += qasm.chain_strength

def assign_pin_strength(pin_str, qubo):
    "Define a strength for each explicitly pinned qubit."
    qasm.pin_strength = pin_str
    if qasm.pin_strength == None:
        # Pin strength defaults to half the chain strength.
        qasm.pin_strength = qasm.chain_strength/2.0
    elif qubo:
        # With QUBO input we need to divide the chain strength by 4 for
        # consistency with the other coupler strengths.
        qasm.pin_strength /= 4.0
