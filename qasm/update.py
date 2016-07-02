#####################################
# Update QASM weights and strengths #
# By Scott Pakin <pakin@lanl.gov>   #
#####################################

import qasm
import random
import string
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

def new_internal_sym():
    "Create a new internal symbol."
    while True:
        sym = "$"
        for i in range(5):
            sym += random.choice(string.lowercase)
        if not qasm.sym2num.has_key(sym):
            return sym

def pin_qubits():
    "Use a helper qubit to help pin values to true or false."
    for q_user, b in qasm.pinned:
        q_helper = qasm.symbol_to_number(new_internal_sym())
        q1, q2 = q_helper, q_user
        if q1 > q2:
            q1, q2 = q2, q1
        if b:
            qasm.weights[q_helper] += qasm.pin_strength/2.0
            qasm.weights[q_user] += qasm.pin_strength
            qasm.strengths[(q1, q2)] += -qasm.pin_strength/2.0
        else:
            qasm.weights[q_helper] += qasm.pin_strength/2.0
            qasm.weights[q_user] += -qasm.pin_strength
            qasm.strengths[(q1, q2)] += qasm.pin_strength/2.0

def convert_chains_to_aliases():
    "Convert chains to aliases where possible."
    # Identify all chains that can be converted to aliases.
    num2allsyms = [[] for _ in range(len(qasm.sym2num))]
    for s, n in qasm.sym2num.items():
        num2allsyms[n].append(s)
    make_aliases = []
    for q1, q2 in qasm.chains:
        if qasm.weights[q1] == qasm.weights[q2]:
            make_aliases.append((q1, q2))
    make_aliases.sort(reverse=True, key=lambda qs: (qs[1], qs[0]))

    # Replace each chain in make_aliases with an alias.  Work in reverse order
    # of qubit number and shift all greater qubit numbers downward.
    for q1, q2 in make_aliases:
        # Map q2's symbolic names to q1's.  Shift everything above q2 downwards.
        alias_sym2num = {}
        for s, sq in qasm.sym2num.items():
            if sq == q2:
                sq = q1
            elif sq > q2:
                sq -= 1
            alias_sym2num[s] = sq
        qasm.sym2num = alias_sym2num

        # Elide q2 from the list of weights.
        alias_weights = defaultdict(lambda: 0.0)
        for wq, wt in qasm.weights.items():
            if wq == q2:
                continue
            if wq > q2:
                wq -= 1
            if wt != 0.0:
                alias_weights[wq] = wt
        alias_weights[q1] += qasm.weights[q2]   # Conserve overall energy.
        qasm.weights = alias_weights

        # Replace q2 with q1 in all strengths.  Shift everything above q2
        # downwards.
        alias_strengths = defaultdict(lambda: 0.0)
        for (sq1, sq2), wt in qasm.strengths.items():
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
        qasm.strengths = alias_strengths

        # Replace q2 with q1 in all strengths.  Shift everything above q2
        # downwards.
        alias_chains = {}
        for cq1, cq2 in qasm.chains.keys():
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
        qasm.chains = alias_chains

        # Replace q2 with q1 in all pinned qubits.  Shift everything above q2
        # downwards.
        alias_pinned = []
        for pq, b in qasm.pinned:
            if pq == q2:
                pq = q1
            if pq > q2:
                pq -= 1
            alias_pinned.append((pq, b))
        qasm.pinned = alias_pinned

        # We now have one fewer symbol.
        qasm.next_sym_num -= 1
