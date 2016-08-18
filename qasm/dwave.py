########################################
# Manage D-Wave communication for QASM #
# By Scott Pakin <pakin@lanl.gov>      #
########################################

from dwave_sapi2.core import solve_ising
from dwave_sapi2.embedding import find_embedding, embed_problem, unembed_answer
from dwave_sapi2.local import local_connection
from dwave_sapi2.remote import RemoteConnection
from dwave_sapi2.util import get_hardware_adjacency
import copy
import math
import os
import qasm
import re
import sys

def connect_to_dwave():
    """
    Establish a connection to the D-Wave, and use this to talk to a solver.
    We rely on the qOp infrastructure to set the environment variables properly.
    """
    try:
        url = os.environ["DW_INTERNAL__HTTPLINK"]
        token = os.environ["DW_INTERNAL__TOKEN"]
        conn = RemoteConnection(url, token)
    except KeyError:
        url = "<local>"
        token = "<N/A>"
        conn = local_connection
    except IOError as e:
        qasm.abend("Failed to establish a remote connection (%s)" % e)
    try:
        qasm.solver_name = os.environ["DW_INTERNAL__SOLVER"]
    except:
        # Solver was not specified: Use the first available solver.
        qasm.solver_name = conn.solver_names()[0]
    try:
        qasm.solver = conn.get_solver(qasm.solver_name)
    except KeyError:
        qasm.abend("Failed to find solver %s on connection %s" % (qasm.solver_name, url))

def find_dwave_embedding(logical, verbosity):
    """Find an embedding of a logical problem in the D-Wave's physical topology.
    Store the embedding within the Problem object."""
    edges = logical.strengths.keys()
    edges.sort()
    try:
        hw_adj = get_hardware_adjacency(qasm.solver)
    except KeyError:
        # The Ising heuristic solver is an example of a solver that lacks a fixed
        # hardware representation.  We therefore assert that the hardware exactly
        # matches the problem'input graph.
        hw_adj = edges
    if verbosity >= 2:
        sys.stderr.write("Embedding the logical adjacency within the physical topology.\n\n")
        stdout_fd = os.dup(sys.stdout.fileno())
        os.dup2(sys.stderr.fileno(), sys.stdout.fileno())
        embedding = find_embedding(edges, hw_adj, verbose=1)
        os.dup2(stdout_fd, sys.stdout.fileno())
        sys.stderr.write("\n")
    else:
        embedding = find_embedding(edges, hw_adj, verbose=0)
    if embedding == []:
        # We received an empty embedding.  I've seen this happen with the
        # ising-heuristic solver.  A workaround seems to be to fabricate a
        # trivial embedding in which logical qubit X maps to physical qubit X.
        embedding = [[q] for q in range(qasm.next_sym_num + 1)]
    logical.embedding = embedding
    logical.hw_adj = hw_adj
    logical.edges = edges

def embed_problem_on_dwave(logical, verbosity):
    """Embed a logical problem in the D-Wave's physical topology.  Return a
    physical Problem object."""
    # Embed the problem.  Abort on failure.
    find_dwave_embedding(logical, verbosity)
    try:
        h_range = qasm.solver.properties["h_range"]
        j_range = qasm.solver.properties["j_range"]
    except KeyError:
        h_range = [-1.0, 1.0]
        j_range = [-1.0, 1.0]
    weight_list = [logical.weights[q] for q in range(qasm.next_sym_num + 1)]
    smearable = any([s != 0.0 for s in logical.strengths.values()])
    try:
        [new_weights, new_strengths, new_chains, new_embedding] = embed_problem(
            weight_list, logical.strengths, logical.embedding, logical.hw_adj,
            True, smearable, h_range, j_range)
    except ValueError as e:
        qasm.abend("Failed to embed the problem in the solver (%s)" % e)

    # Construct a physical Problem object.
    physical = copy.deepcopy(logical)
    physical.chains = new_chains
    physical.embedding = new_embedding
    physical.h_range = h_range
    physical.j_range = j_range
    physical.strengths = new_strengths
    physical.weight_list = weight_list
    physical.weights = new_weights
    physical.pinned = []
    for l, v in logical.pinned:
        physical.pinned.extend([(p, v) for p in physical.embedding[l]])
    return physical

def optimize_dwave_layout(logical, physical, verbosity):
    "Iteratively search for a better embedding.  Return a new physical Problem."
    try:
        # Say what we're about to do
        if verbosity >= 2:
            sys.stderr.write("Optimizing the embedding.\n\n")

        # Determine the edges of a rectangle of cells we want to use.
        L, M, N = qasm.chimera_topology(qasm.solver)
        L2 = 2*L
        ncells = (qasm.next_sym_num + 1) // L2
        edgey = max(int(math.sqrt(ncells)), 1)
        edgex = max((ncells + edgey - 1) // edgey, 1)

        # Repeatedly expand edgex and edgey until the embedding works.
        while edgex <= M and edgey <= N:
            # Retain adjacencies only within the rectangle.
            alt_hw_adj = []
            for q1, q2 in physical.hw_adj:
                c1 = q1//L2
                if c1 % M >= edgex:
                    continue
                if c1 // M >= edgey:
                    continue
                c2 = q2//L2
                if c2 % M >= edgex:
                    continue
                if c2 // M >= edgey:
                    continue
                alt_hw_adj.append((q1, q2))
            alt_hw_adj = set(alt_hw_adj)

            # Embed again.
            try:
                if verbosity >= 2:
                    sys.stderr.write("  Trying a %dx%d unit-cell embedding ... " % (edgex, edgey))
                alt_embedding = find_embedding(physical.edges, alt_hw_adj, verbose=0)
                [new_weights, new_strengths, new_chains, new_embedding] = embed_problem(
                    physical.weight_list, logical.strengths, alt_embedding,
                    alt_hw_adj, True, True, physical.h_range, physical.j_range)
                if verbosity >= 2:
                    sys.stderr.write("succeeded\n")
                break  # Success!
            except ValueError as e:
                # Failure -- increase edgex or edgey and try again.
                if verbosity >= 2:
                    sys.stderr.write("failed\n")
                if edgex < edgey:
                    edgex += 1
                else:
                    edgey += 1
        if not(edgex <= M and edgey <= N):
            qasm.abend("Failed to optimize the embedding")
    except KeyError:
        if verbosity >= 2:
            sys.stderr.write("  - Failed to query the machine topology for embedding parameters\n")
    if verbosity >= 2:
        sys.stderr.write("\n")
    new_physical = copy.deepcopy(physical)
    new_physical.weights = new_weights
    new_physical.strengths = new_strengths
    new_physical.chains = new_chains
    new_physical.embedding = new_embedding
    return new_physical

def update_strengths_from_chains(physical):
    """Update strengths using the chains introduced by embedding.  Return a new
    physical Problem object."""
    new_physical = copy.deepcopy(physical)
    new_physical.chains = {c: qasm.chain_strength for c in physical.chains.keys()}
    new_physical.strengths = physical.strengths.copy()
    new_physical.strengths.update(new_physical.chains)
    return new_physical

def scale_weights_strengths(physical, verbosity):
    "Manually scale the weights and strengths so Qubist doesn't complain."
    h_range = physical.h_range
    j_range = physical.j_range
    old_cap = max([abs(w) for w in physical.weights + physical.strengths.values()])
    new_cap = min(-h_range[0], h_range[1], -j_range[0], j_range[1])
    if old_cap == 0.0:
        # Handle the obscure case of a zero old_cap.
        old_cap = new_cap
    new_weights = [w*new_cap/old_cap for w in physical.weights]
    new_strengths = {js: w*new_cap/old_cap for js, w in physical.strengths.items()}
    if verbosity >= 1 and old_cap != new_cap:
        sys.stderr.write("Scaling weights and strengths from [%.10g, %.10g] to [%.10g, %.10g].\n\n" % (-old_cap, old_cap, -new_cap, new_cap))
    new_physical = copy.deepcopy(physical)
    new_physical.weights = new_weights
    new_physical.strengths = new_strengths
    return new_physical

# Define a function that says whether a solution contains no broken pins and no
# broken (user-specified) chains.
def solution_is_intact(physical, soln):
    # Reject broken pins.
    bool2spin = [-1, +1]
    for pnum, pin in physical.pinned:
        if soln[pnum] != bool2spin[pin]:
            return False

    # Reject broken chains.
    for q1, q2 in physical.chains.keys():
        if soln[q1] != soln[q2]:
            return False

    # The solution looks good!
    return True

def submit_dwave_problem(physical, samples, anneal_time):
    "Submit a QMI to the D-Wave."
    # Submit a QMI to the D-Wave and get back a list of solution vectors.
    solver_params = dict(chains=physical.embedding, num_reads=samples, annealing_time=anneal_time)
    while True:
        # Repeatedly remove parameters the particular solver doesn't like until
        # it actually works -- or fails for a different reason.
        try:
            answer = solve_ising(qasm.solver, physical.weights, physical.strengths, **solver_params)
            break
        except ValueError as e:
            # Is there a better way to extract the failing symbol than a regular
            # expression match?
            bad_name = re.match(r'"(.*?)"', str(e))
            if bad_name == None:
                raise e
            del solver_params[bad_name.group(1)]
        except RuntimeError as e:
            qasm.abend(e)

    # Discard solutions with broken pins or broken chains.
    solutions = answer["solutions"]
    valid_solns = [s for s in solutions if solution_is_intact(physical, s)]
    final_answer = unembed_answer(valid_solns, physical.embedding, broken_chains="discard")
    return answer, final_answer
