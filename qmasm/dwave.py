#########################################
# Manage D-Wave communication for QMASM #
# By Scott Pakin <pakin@lanl.gov>       #
#########################################

from collections import defaultdict
try:
    from dwave_sapi2.core import solve_ising
    from dwave_sapi2.embedding import find_embedding, embed_problem, unembed_answer
    from dwave_sapi2.local import local_connection
    from dwave_sapi2.remote import RemoteConnection
    from dwave_sapi2.util import get_hardware_adjacency
except ImportError:
    from fake_dwave import *
import copy
import hashlib
import marshal
import math
import os
import qmasm
import re
import sys
import tempfile

def connect_to_dwave():
    """
    Establish a connection to the D-Wave, and use this to talk to a solver.
    We rely on the qOp infrastructure to set the environment variables properly.
    """
    try:
        url = os.environ["DW_INTERNAL__HTTPLINK"]
        token = os.environ["DW_INTERNAL__TOKEN"]
        try:
            proxy = os.environ["DW_INTERNAL__HTTPPROXY"]
        except KeyError:
            proxy = ""
        conn = RemoteConnection(url, token, proxy)
    except KeyError:
        url = "<local>"
        token = "<N/A>"
        conn = local_connection
    except IOError as e:
        qmasm.abend("Failed to establish a remote connection (%s)" % e)
    try:
        qmasm.solver_name = os.environ["DW_INTERNAL__SOLVER"]
    except:
        # Solver was not specified: Use the first available solver.
        qmasm.solver_name = conn.solver_names()[0]
    try:
        qmasm.solver = conn.get_solver(qmasm.solver_name)
    except KeyError:
        qmasm.abend("Failed to find solver %s on connection %s" % (qmasm.solver_name, url))

class EmbeddingCache(object):
    "Read and write an embedding cache file."

    def __init__(self, edges, adj):
        # Ensure we have a valid cache directory.
        self.hash = None
        try:
            self.cachedir = os.environ["QMASMCACHE"]
        except KeyError:
            self.cachedir = None
            return None
        if not os.path.isdir(self.cachedir):
            qmasm.abend("QMASMCACHE is set to %s, which is not an extant directory" % self.cachedir)

        # Compute an MD5 sum of our inputs.
        sha = hashlib.sha1()
        sha.update(str(sorted(edges)))
        sha.update(str(sorted(adj)))
        self.hash = sha.hexdigest()

    def read(self):
        "Read an embedding from an embedding cache or None on a cache miss."
        if self.hash == None:
            return None
        try:
            h = open(os.path.join(self.cachedir, self.hash))
        except IOError:
            return None
        embedding = marshal.load(h)
        h.close()
        return embedding

    def write(self, embedding):
        "Write an embedding to an embedding cache."
        if self.hash == None:
            return
        try:
            h = open(os.path.join(self.cachedir, self.hash), "w")
        except IOError:
            return None
        marshal.dump(embedding, h)
        h.close()

def find_dwave_embedding(logical, optimize, verbosity):
    """Find an embedding of a logical problem in the D-Wave's physical topology.
    Store the embedding within the Problem object."""
    # SAPI tends to choke when embed_problem is told to embed a problem
    # containing a zero-weight node whose adjacent couplers all have zero
    # strength.  (Tested with SAPI 2.4.)  To help out SAPI, we simply remove
    # all zero-strength couplers.
    edges = [e for e in logical.strengths.keys() if logical.strengths[e] != 0.0]
    edges.sort()
    logical.edges = edges
    try:
        hw_adj = get_hardware_adjacency(qmasm.solver)
    except KeyError:
        # The Ising heuristic solver is an example of a solver that lacks a
        # fixed hardware representation.  We therefore assert that the hardware
        # is an all-to-all network that connects every node to every other node.
        endpoints = set([a for a, b in edges] + [b for a, b in edges])
        hw_adj = [(a, b) for a in endpoints for b in endpoints if a != b]

    # Determine the edges of a rectangle of cells we want to use.
    L, M, N = qmasm.chimera_topology(qmasm.solver)
    L2 = 2*L
    ncells = (qmasm.next_sym_num + L2) // L2   # Round up the number of cells.
    if optimize:
        edgey = max(int(math.sqrt(ncells)), 1)
        edgex = max((ncells + edgey - 1) // edgey, 1)
    else:
        edgey = N
        edgex = M

    # Announce what we're about to do.
    if verbosity >= 2:
        sys.stderr.write("Embedding the logical adjacency within the physical topology.\n\n")

    # Repeatedly expand edgex and edgey until the embedding works.
    while edgex <= M and edgey <= N:
        # Retain adjacencies only within the rectangle.
        alt_hw_adj = []
        for q1, q2 in hw_adj:
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
        logical.hw_adj = alt_hw_adj

        # See if we already have an embedding in the embedding cache.
        ec = EmbeddingCache(edges, alt_hw_adj)
        if verbosity >= 2:
            if ec.cachedir == None:
                sys.stderr.write("  No embedding cache directory was specified ($QMASMCACHE).\n")
            else:
                sys.stderr.write("  Using %s as the embedding cache directory ...\n" % ec.cachedir)
        embedding = ec.read()
        if embedding == []:
            # Cache hit, but embedding had failed
            if verbosity >= 2:
                sys.stderr.write("  Found failed embedding %s in the embedding cache.\n\n" % ec.hash)
        elif embedding != None:
            # Successful cache hit!
            if verbosity >= 2:
                sys.stderr.write("  Found successful embedding %s in the embedding cache.\n\n" % ec.hash)
            logical.embedding = embedding
            return
        if verbosity >= 2 and ec.cachedir != None:
            sys.stderr.write("  No existing embedding found in the embedding cache.\n")

        # Try to find an embedding, unless we previously determined that it had
        # failed.
        if embedding != []:
            if verbosity >= 2:
                sys.stderr.write("  Trying a %dx%d unit-cell embedding ... " % (edgex, edgey))
                status_file = tempfile.NamedTemporaryFile(mode="w", delete=False)
                stdout_fd = os.dup(sys.stdout.fileno())
                os.dup2(status_file.fileno(), sys.stdout.fileno())
                embedding = find_embedding(edges, alt_hw_adj, verbose=1)
                sys.stdout.flush()
                os.dup2(stdout_fd, sys.stdout.fileno())
                status_file.close()
                if len(embedding) > 0:
                    sys.stderr.write("succeeded\n\n")
                else:
                    sys.stderr.write("failed\n\n")
                with open(status_file.name, "r") as status:
                    for line in status:
                        sys.stderr.write("    %s" % line)
                sys.stderr.write("\n")
                os.remove(status_file.name)
            else:
                embedding = find_embedding(edges, alt_hw_adj, verbose=0)
            ec.write(embedding)
            if len(embedding) > 0:
                # Success!
                break

        # Continued failure -- increase edgex or edgey and try again.
        if edgex < edgey:
            edgex += 1
        else:
            edgey += 1
    if not(edgex <= M and edgey <= N):
        qmasm.abend("Failed to embed the problem")
    logical.embedding = embedding

def embed_problem_on_dwave(logical, optimize, verbosity):
    """Embed a logical problem in the D-Wave's physical topology.  Return a
    physical Problem object."""
    # Embed the problem.  Abort on failure.
    find_dwave_embedding(logical, optimize, verbosity)
    try:
        h_range = qmasm.solver.properties["h_range"]
        j_range = qmasm.solver.properties["j_range"]
    except KeyError:
        h_range = [-1.0, 1.0]
        j_range = [-1.0, 1.0]
    weight_list = qmasm.dict_to_list(logical.weights)
    smearable = any([s != 0.0 for s in logical.strengths.values()])
    try:
        [new_weights, new_strengths, new_chains, new_embedding] = embed_problem(
            weight_list, logical.strengths, logical.embedding, logical.hw_adj,
            True, smearable, h_range, j_range)
    except ValueError as e:
        qmasm.abend("Failed to embed the problem in the solver (%s)" % e)

    # Construct a physical Problem object.
    physical = copy.deepcopy(logical)
    physical.chains = new_chains
    physical.embedding = new_embedding
    physical.h_range = h_range
    physical.j_range = j_range
    physical.strengths = new_strengths
    physical.weights = defaultdict(lambda: 0.0,
                                   {q: new_weights[q]
                                    for q in range(len(new_weights))
                                    if new_weights[q] != 0.0})
    physical.pinned = []
    for l, v in logical.pinned:
        physical.pinned.extend([(p, v) for p in physical.embedding[l]])
    return physical

def update_strengths_from_chains(physical):
    """Update strengths using the chains introduced by embedding.  Return a new
    physical Problem object."""
    new_physical = copy.deepcopy(physical)
    new_physical.chains = {c: qmasm.chain_strength for c in physical.chains.keys()}
    new_physical.strengths = physical.strengths.copy()
    new_physical.strengths.update(new_physical.chains)
    return new_physical

def scale_weights_strengths(physical, verbosity):
    "Manually scale the weights and strengths so Qubist doesn't complain."
    h_range = physical.h_range
    j_range = physical.j_range
    weight_list = qmasm.dict_to_list(physical.weights)
    old_cap = max([abs(w) for w in weight_list + physical.strengths.values()])
    new_cap = min(-h_range[0], h_range[1], -j_range[0], j_range[1])
    if old_cap == 0.0:
        # Handle the obscure case of a zero old_cap.
        old_cap = new_cap
    new_weights = qmasm.list_to_dict([w*new_cap/old_cap for w in weight_list])
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

def submit_dwave_problem(verbosity, physical, samples, anneal_time, spin_revs, postproc):
    "Submit a QMI to the D-Wave."
    # Map abbreviated to full names for postprocessing types.
    postproc = {"": "", "opt": "optimization", "sample": "sampling"}[postproc]

    # Submit a QMI to the D-Wave and get back a list of solution vectors.
    solver_params = dict(chains=physical.embedding,
                         num_reads=samples,
                         annealing_time=anneal_time,
                         num_spin_reversal_transforms=spin_revs,
                         postprocess=postproc)
    unused_params = dict()
    while True:
        # Repeatedly remove parameters the particular solver doesn't like until
        # it actually works -- or fails for a different reason.
        try:
            weight_list = qmasm.dict_to_list(physical.weights)
            answer = solve_ising(qmasm.solver, weight_list, physical.strengths, **solver_params)
            break
        except ValueError as e:
            # Is there a better way to extract the failing symbol than a regular
            # expression match?
            bad_name_match = re.match(r'"(.*?)"', str(e))
            if bad_name_match == None:
                raise e
            bad_name = bad_name_match.group(1)
            unused_params[bad_name] = solver_params[bad_name]
            del solver_params[bad_name]
        except RuntimeError as e:
            qmasm.abend(e)
    if verbosity >= 2:
        # Output parameters we kept and those we discarded
        sys.stderr.write("Parameters accepted by the %s solver:\n" % qmasm.solver_name)
        if len(solver_params) > 0:
            for k, v in solver_params.items():
                sys.stderr.write("    %s = %s\n" % (k, v))
        else:
            sys.stderr.write("    [none]\n")
        sys.stderr.write("\n")
        sys.stderr.write("Parameters rejected by the %s solver:\n" % qmasm.solver_name)
        if len(unused_params) > 0:
            for k, v in unused_params.items():
                sys.stderr.write("    %s = %s\n" % (k, v))
        else:
            sys.stderr.write("    [none]\n")
        sys.stderr.write("\n")

    # Tally the occurrences of each solution
    solutions = answer["solutions"]
    semifinal_answer = unembed_answer(solutions, physical.embedding, broken_chains="vote")
    try:
        num_occurrences = {tuple(k): v
                           for k, v in zip(semifinal_answer, answer["num_occurrences"])}
    except KeyError:
        num_occurrences = {tuple(a): 1 for a in semifinal_answer}

    # Discard solutions with broken pins or broken chains.
    valid_solns = [s for s in solutions if solution_is_intact(physical, s)]
    final_answer = unembed_answer(valid_solns, physical.embedding, broken_chains="discard")
    return answer, final_answer, num_occurrences
