###############################################
# Fake a minimal set of D-Wave SAPI functions #
# By Scott Pakin <pakin@lanl.gov>             #
###############################################

import qmasm

class FakeSolver(object):
    properties = {}

class FakeConnection(object):
    def solver_names(self):
        return ["phony"]

    def get_solver(self, sname):
        if sname == "phony":
            return FakeSolver()
        else:
            raise KeyError

local_connection = FakeConnection()

def RemoteConnection(url, token, proxy=None):
    qmasm.abend("In the absence of D-Wave's libraries, remote connections are unavailable")

def ising_to_qubo(hs, js):
    "Convert a list of hs and a dictionary of Js to a dictionary of Qs."
    qs = {}

    # Compute the new point weights along the diagonal.
    for (i, j), s in js.items():
        qs[(i, i)] = 0.0
        qs[(j, j)] = 0.0
        qs[(i, j)] = 0.0
    for i in range(len(hs)):
        qs[(i, i)] = hs[i]*2

    # Compute the new strengths off the diagonal.
    for (i, j), s in js.items():
        qs[(i, i)] -= 2.0*s
        qs[(j, j)] -= 2.0*s
        qs[(i, j)] += 4.0*s

    # Discard zeroes.
    qs = {k: v for k, v in qs.items() if v != 0.0}

    # Compute an energy offset.
    qoffset = sum(hs) + sum(js.values())

    # Return the QUBO matrix and energy offset.
    return qs, qoffset

def qubo_to_ising(qs):
    "Convert a dictionary of Qs to a list of hs and a dictionary of Js."
    # Initialize the hs and js dictionaries.
    hs = {}  # We'll convert to a list later.
    js = {}
    for (i, j), s in qs.items():
        hs[i] = 0.0
        if i != j:
            hs[j] = 0.0
            js[(i, j)] = 0.0

    # Peform an initial conversion.
    for (i, j), s in qs.items():
        if i == j:
            # Point weight
            hs[i] += s/2.0
        else:
            # Coupler strength
            js[(i, j)] += s/4.0
            hs[i] += s/4.0
            hs[j] += s/4.0

    # Convert hs to a list and elide zeroes from js.
    try:
        mh = max(hs.keys())
    except ValueError:
        mh = -1
    hlist = [0] * (mh + 1)
    for i, s in hs.items():
        hlist[i] = s
    js = {k: v for k, v in js.items() if v != 0.0}

    # QMASM doesn't use an offset so ignore it.
    return hlist, js, None

def get_hardware_adjacency(solver):
    qmasm.abend("Without D-Wave's libraries, QMASM can do little more than output qbsolv, MiniZinc, and flattened QMASM files")

def fix_variables(Q, method="optimized"):
    "Return the original QUBO unmodified."
    return {"new_Q": Q,
            "offset": 0.0,
            "fixed_variables": {}}

def unembed_answer(solutions, embeddings, broken_chains=None, h=None, j=None):
    "Return the original solutions unmodified."
    return solutions
