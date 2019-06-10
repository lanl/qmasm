#########################################
# Manage D-Wave communication for QMASM #
# By Scott Pakin <pakin@lanl.gov>       #
#########################################

from collections import defaultdict
try:
    from dwave_sapi2.core import async_solve_ising, await_completion
    from dwave_sapi2.embedding import find_embedding, embed_problem, unembed_answer
    from dwave_sapi2.fix_variables import fix_variables
    from dwave_sapi2.local import local_connection
    from dwave_sapi2.remote import RemoteConnection
    from dwave_sapi2.util import get_hardware_adjacency, ising_to_qubo, qubo_to_ising
except ImportError:
    from .fake_dwave import *
import copy
import hashlib
import json
import marshal
import math
import operator
import os
import qmasm
import re
import sys
import tempfile
import time

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

        # Compute a SHA-1 sum of our inputs.
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

def report_embeddability(edges, adj):
    """Output some metrics on how likely a set of edges can be embedded in
    a given adjacency graph."""
    # Output a histogram of node degrees.
    embed, extras, (nvars, nqubits), (nstrs, ncouplers), deg_hist = qmasm.maybe_embeddable(edges, adj)
    sys.stderr.write("Embeddability metrics:\n\n")
    sys.stderr.write("    Degree  Tally\n")
    sys.stderr.write("    ------  -----\n")
    for d_h in deg_hist:
        sys.stderr.write("    %6d  %5d\n" % d_h)
    sys.stderr.write("\n")

    # Output the largest possible clique.
    max_d_h = None
    for d, h in deg_hist:
        if h > d:
            max_d_h = (d, h)
    if max_d_h == None:
        sys.stderr.write("    Largest possible clique: none\n\n")
    else:
        sys.stderr.write("    Largest possible clique: K_%d\n\n" % max_d_h[0])

    # Report whether embedding is impossible.
    sys.stderr.write("    Minimum qubit usage:     (%5d + %5d) / %5d = %8.2f%%\n" % \
                     (nvars, extras, nqubits, (nvars + extras)*100.0/nqubits))
    sys.stderr.write("    Minimum coupler usage:   (%5d + %5d) / %5d = %8.2f%%\n" % \
                     (nstrs, extras, ncouplers, (nstrs + extras)*100.0/ncouplers))
    if embed:
        sys.stderr.write("    Embedding is impossible: NO\n")
    else:
        sys.stderr.write("    Embedding is impossible: YES\n")
    sys.stderr.write("\n")

def read_hardware_adjacency(fname, verbosity):
    """Read a hardware adjacency list from a file.  Each line must contain
    a space-separated pair of vertex numbers."""
    adj = set()
    lineno = 0
    if verbosity >= 2:
        sys.stderr.write("Reading hardware adjacency from %s ... " % fname)
    with open(fname) as f:
        for orig_line in f:
            # Discard comments then parse the line into exactly two vertices.
            lineno += 1
            line = orig_line.partition("#")[0]
            try:
                verts = [int(v) for v in line.split()]
            except ValueError:
                qmasm.abend('Failed to parse line %d of file %s ("%s")' % (lineno, fname, orig_line.strip()))
            if len(verts) == 0:
                continue
            if len(verts) != 2 or verts[0] == verts[1]:
                qmasm.abend('Failed to parse line %d of file %s ("%s")' % (lineno, fname, orig_line.strip()))

            # Canonicalize the vertex numbers and add the result to the set.
            if verts[1] > verts[0]:
                verts[0], verts[1] = verts[1], verts[0]
            adj.add((verts[0], verts[1]))
    if verbosity >= 2:
        sys.stderr.write("%d unique edges found\n\n" % len(adj))
    return sorted(adj)

def qubo_vars(Q):
    "Return a set of all variables named by a QUBO."
    vars = set()
    for q1, q2 in Q:
        vars.add(q1)
        vars.add(q2)
    return vars

def simplify_problem(logical, verbosity):
    """Try to find spins that can be removed from the problem because their
    value is known a priori."""
    # SAPI's fix_variables function works only on QUBOs so we have to convert.
    # We directly use SAPI's ising_to_qubo function instead of our own
    # convert_to_qubo because the QUBO has to be in matrix form.
    hs = qmasm.dict_to_list(logical.weights)
    Js = logical.strengths
    Q, qubo_offset = ising_to_qubo(hs, Js)

    # Simplify the problem if possible.
    simple = fix_variables(Q, method="standard")
    new_Q = simple["new_Q"]
    fixed_vars = simple["fixed_variables"]
    if verbosity >= 2:
        # Also determine if we could get rid of more qubits if we care about
        # only *a* solution rather than *all* solutions.
        alt_simple = fix_variables(Q, method="optimized")
        all_gone = len(alt_simple["new_Q"]) == 0

    # Work around the rare case in which fix_variables drops a variable
    # entirely, leaving it neither in new_Q nor in fixed_variables.  If this
    # happenes, we explicitly re-add the variable from Q to new_Q and
    # transitively everything it touches (removing from fixed_vars if a
    # variable appears there).
    old_vars = qubo_vars(Q)
    new_vars = qubo_vars(new_Q)
    new_vars.update(fixed_vars)
    missing_vars = sorted(old_vars.difference(new_vars))
    while len(missing_vars) > 0:
        q = missing_vars.pop()
        for (q1, q2), val in Q.items():
            if q1 == q or q2 == q:
                new_Q[(q1, q2)] = val
                fixed_vars.pop(q1, None)
                fixed_vars.pop(q2, None)
                if q1 == q and q2 > q:
                    missing_vars.append(q2)
                elif q2 == q and q1 > q:
                    missing_vars.append(q1)

    # At high verbosity levels, list all of the known symbols and their value.
    if verbosity >= 2:
        # Map each logical qubit to one or more symbols.
        num2syms = [[] for _ in range(qmasm.sym_map.max_number() + 1)]
        max_sym_name_len = 7
        for q, n in qmasm.sym_map.symbol_number_items():
            num2syms[n].append(q)
            max_sym_name_len = max(max_sym_name_len, len(repr(num2syms[n])) - 1)

        # Output a table of know values
        sys.stderr.write("Elided qubits whose low-energy value can be determined a priori:\n\n")
        if len(fixed_vars) > 0:
            sys.stderr.write("    Logical  %-*s  Value\n" % (max_sym_name_len, "Name(s)"))
            sys.stderr.write("    -------  %s  -----\n" % ("-" * max_sym_name_len))
            truval = {0: "False", +1: "True"}
            for q, b in sorted(fixed_vars.items()):
                try:
                    syms = qmasm.sym_map.to_symbols(q)
                except KeyError:
                    continue
                name_list = " ".join(sorted(syms))
                sys.stderr.write("    %7d  %-*s  %-s\n" % (q, max_sym_name_len, name_list, truval[b]))
            sys.stderr.write("\n")

    # Return the original problem if no qubits could be elided.
    if verbosity >= 2:
        sys.stderr.write("  %6d logical qubits before elision\n" % (qmasm.sym_map.max_number() + 1))
    if len(fixed_vars) == 0:
        if verbosity >= 2:
            sys.stderr.write("  %6d logical qubits after elision\n\n" % (qmasm.sym_map.max_number() + 1))
            if all_gone:
                sys.stderr.write("    Note: A complete solution can be found classically using roof duality and strongly connected components.\n\n")
        return logical

    # Construct a simplified problem, renumbering so as to compact qubit
    # numbers.
    new_obj = copy.deepcopy(logical)
    new_obj.known_values = {s: 2*fixed_vars[n] - 1
                            for s, n in qmasm.sym_map.symbol_number_items()
                            if n in fixed_vars}
    new_obj.simple_offset = simple["offset"]
    hs, Js, ising_offset = qubo_to_ising(new_Q)
    qubits_used = set([i for i in range(len(hs)) if hs[i] != 0.0])
    for q1, q2 in Js.keys():
        qubits_used.add(q1)
        qubits_used.add(q2)
    qmap = dict(zip(sorted(qubits_used), range(len(qubits_used))))
    new_obj.chains = set([(qmap[q1], qmap[q2])
                          for q1, q2 in new_obj.chains
                          if q1 in qmap and q2 in qmap])
    new_obj.antichains = set([(qmap[q1], qmap[q2])
                              for q1, q2 in new_obj.antichains
                              if q1 in qmap and q2 in qmap])
    new_obj.weights = defaultdict(lambda: 0.0,
                                  {qmap[i]: hs[i]
                                   for i in range(len(hs))
                                   if hs[i] != 0.0})
    new_obj.strengths = qmasm.canonicalize_strengths({(qmap[q1], qmap[q2]): wt
                                                      for (q1, q2), wt in Js.items()})
    new_obj.pinned = [(qmap[q], b)
                      for q, b in new_obj.pinned
                      if q in qmap]
    qmasm.sym_map.overwrite_with({s: qmap[q]
                                  for s, q in qmasm.sym_map.symbol_number_items()
                                  if q in qmap})
    if verbosity >= 2:
        # Report the number of logical qubits that remain, but compute the
        # number that could be removed if only a single solution were required.
        sys.stderr.write("  %6d logical qubits after elision\n\n" % (qmasm.sym_map.max_number() + 1))
        if qmasm.sym_map.max_number() > -1 and all_gone:
            sys.stderr.write("    Note: A complete solution can be found classically using roof duality and strongly connected components.\n\n")
    return new_obj

def find_dwave_embedding(logical, optimization, verbosity, hw_adj_file):
    """Find an embedding of a logical problem in the D-Wave's physical topology.
    Store the embedding within the Problem object."""
    # SAPI tends to choke when embed_problem is told to embed a problem
    # containing a zero-weight node whose adjacent couplers all have zero
    # strength.  (Tested with SAPI 2.4.)  To help out SAPI, we simply remove
    # all zero-strength couplers.
    edges = [e for e in logical.strengths.keys() if logical.strengths[e] != 0.0]
    edges.sort()
    logical.edges = edges
    if hw_adj_file == None:
        try:
            hw_adj = get_hardware_adjacency(qmasm.solver)
        except KeyError:
            # The Ising heuristic solver is an example of a solver that lacks a
            # fixed hardware representation.  We therefore assert that the
            # hardware is an all-to-all network that connects every node to
            # every other node.
            endpoints = set([a for a, b in edges] + [b for a, b in edges])
            hw_adj = [(a, b) for a in endpoints for b in endpoints if a != b]
    else:
        hw_adj = read_hardware_adjacency(hw_adj_file, verbosity)

    # Tell the user if we have any hope at all of embedding the problem.
    if verbosity >= 2:
        report_embeddability(edges, hw_adj)

    # Determine the edges of a rectangle of cells we want to use.  If we read
    # the topology from a file or otherwise can't prove that we have a Chimera
    # graph, we call this rectangle 0x0 and force the main embedding loop to
    # exit after a single iteration because we don't know the topology is even
    # rectangular.
    edgex = 0
    edgey = 0
    M = 0
    N = 0
    num_vars = len(qmasm.sym_map.all_numbers())
    try:
        if hw_adj_file == None:
            L, M, N = qmasm.chimera_topology(qmasm.solver)
            L2 = 2*L
            ncells = (num_vars + L2) // L2   # Round up the number of cells.
            if optimization >= 2:
                edgey = max(int(math.sqrt(ncells)), 1)
                edgex = max((ncells + edgey - 1) // edgey, 1)
            else:
                edgey = N
                edgex = M
    except qmasm.NonChimera:
        pass

    # Announce what we're about to do.
    if verbosity >= 2:
        sys.stderr.write("Embedding the logical adjacency within the physical topology.\n\n")

    # Repeatedly expand edgex and edgey until the embedding works.
    while edgex <= M and edgey <= N:
        if edgex == M and edgey == N:
            alt_hw_adj = hw_adj
        else:
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
                # SAPI's find_embedding is hard-wired to write to stdout.
                # Trick it into writing into a pipe instead.
                if edgex == 0 and edgey == 0:
                    sys.stderr.write("  Trying to embed ... ")
                else:
                    sys.stderr.write("  Trying a %dx%d unit-cell embedding ...\n\n" % (edgex, edgey))
                sepLine = "=== EMBEDDING ===\n"
                r, w = os.pipe()
                pid = os.fork()
                if pid == 0:
                    # Child -- perform the embedding.
                    os.close(r)
                    os.dup2(w, sys.stdout.fileno())
                    embedding = find_embedding(edges, alt_hw_adj, verbose=1)
                    sys.stdout.flush()
                    os.write(w, sepLine)
                    os.write(w, json.dumps(embedding) + "\n")
                    os.close(w)
                    os._exit(0)
                else:
                    # Parent -- report the embedding's progress.
                    os.close(w)
                    pipe = os.fdopen(r, "r", 10000)
                    while True:
                        try:
                            rstr = pipe.readline()
                            if rstr == sepLine:
                                break
                            if rstr == "":
                                qmasm.abend("Embedder failed to terminate properly")
                            sys.stderr.write("      %s" % rstr)
                        except:
                            pass

                    # Receive the embedding from the child.
                    embedding = json.loads(pipe.readline())
                    sys.stderr.write("\n")
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

def embed_problem_on_dwave(logical, optimization, verbosity, hw_adj_file):
    """Embed a logical problem in the D-Wave's physical topology.  Return a
    physical Problem object."""
    # Embed the problem.  Abort on failure.
    find_dwave_embedding(logical, optimization, verbosity, hw_adj_file)
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
    physical.embedding = new_embedding
    physical.embedder_chains = set(new_chains)
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

# Determine a suitable annealing time to use if none was specified.
def get_default_annealing_time():
    try:
        # Use the default value.
        anneal_time = qmasm.solver.properties["default_annealing_time"]
    except KeyError:
        try:
            # If the default value is undefined, use the minimum allowed
            # value.
            anneal_time = qmasm.solver.properties["annealing_time_range"][0]
        except KeyError:
            # If all else fails, use 20 as a reasonable default.
            anneal_time = 20
    return anneal_time

def compute_sample_counts(samples, anneal_time):
    "Return a list of sample counts to request of the hardware."
    # The number of samples multiplied by the time per sample can't exceed the
    # maximum run duration.
    try:
        # Newest SAPI
        max_run_duration = qmasm.solver.properties["problem_run_duration_range"][1]
    except KeyError:
        try:
            # Slightly older SAPI versions
            max_run_duration = qmasm.solver.properties["max_run_duration"]
        except KeyError:
            # Very old SAPI versions
            max_run_duration = 3000000
    try:
        therm_time = qmasm.solver.properties["default_programming_thermalization"]
    except KeyError:
        therm_time = 0   # As good a guess as any
    max_samples = (max_run_duration - therm_time)//anneal_time

    # The number of samples can't exceed the maximum the hardware allows.
    try:
        max_samples = min(max_samples, qmasm.solver.properties["num_reads_range"][1])
    except KeyError:
        pass

    # Split the number of samples into pieces of size max_samples.
    if samples <= max_samples:
        samples_list = [samples]
    else:
        samples_list = []
        s = samples
        while s > 0:
            if s >= max_samples:
                samples_list.append(max_samples)
                s -= max_samples
            else:
                samples_list.append(s)
                s = 0
    return samples_list

def compute_spin_rev_counts(spin_revs, samples_list):
    "Divide a total number of spin reversals across a set of samples."
    samples = sum(samples_list)
    spin_rev_frac = float(spin_revs)/float(samples)  # We'll evenly divide our spin reversals among samples.
    spin_rev_list = [int(samps*spin_rev_frac) for samps in samples_list]
    missing_srs = spin_revs - sum(spin_rev_list)
    while missing_srs > 0:
        # Account for rounding errors by adding back in some missing spin
        # reversals.
        for i in range(samples):
            if samples_list[i] > spin_rev_list[i]:
                spin_rev_list[i] += 1
                missing_srs -= 1
                if missing_srs == 0:
                    break
    return spin_rev_list

def report_parameters_used(solver_params, unused_params):
    "Output parameters we kept and those we discarded."
    sys.stderr.write("Parameters accepted by the %s solver:\n\n" % qmasm.solver_name)
    if len(solver_params) > 0:
        for k in solver_params.keys():
            sys.stderr.write("    %s\n" % k)
    else:
        sys.stderr.write("    [none]\n")
    sys.stderr.write("\n")
    sys.stderr.write("Parameters rejected by the %s solver:\n\n" % qmasm.solver_name)
    if len(unused_params) > 0:
        for k in unused_params.keys():
            sys.stderr.write("    %s\n" % k)
    else:
        sys.stderr.write("    [none]\n")
    sys.stderr.write("\n")

def report_subproblems_submitted(nqmis, problems, samples_list, spin_rev_list):
    "Output the problem ID for each subproblem submitted."
    sys.stderr.write("Subproblems submitted:\n\n")
    tot_samps = sum([samples_list[i] for i in range(nqmis)])
    tot_sp_revs = sum([spin_rev_list[i] for i in range(nqmis)])
    samp_digs = max(len("Samples"), len(str(tot_samps)))
    sr_digs = max(len("Spin reversals"), len(str(tot_sp_revs)))
    sys.stderr.write("    Problem ID                            %-*s  %-*s\n" %
                     (samp_digs, "Samples", sr_digs, "Spin reversals"))
    sys.stderr.write("    ------------------------------------  %-*s  %-*s\n" %
                     (samp_digs, "-"*samp_digs, sr_digs, "-"*sr_digs))
    for i in range(nqmis):
        status = problems[i].status()
        prob_id = status["problem_id"]
        samps = samples_list[i]
        sp_revs = spin_rev_list[i]
        sys.stderr.write("    %-36s  %*d  %*d\n" %
                         (prob_id, samp_digs, samps, sr_digs, sp_revs))
    if nqmis > 1:
        all_str = "All %d problem IDs" % nqmis
        sys.stderr.write("    %-36s  %*d  %*d\n" %
                         (all_str, samp_digs, tot_samps, sr_digs, tot_sp_revs))
    sys.stderr.write("\n")

def merge_answers(answers):
    "Merge a list of answers into a single, combined piece of information."
    # Handle the trivial case of a single answer.
    n_ans = len(answers)
    if n_ans == 1:
        return answers[0]

    # Merge identical solutions.  Update energies and num_occurrences
    # accordingly.
    solutions = []
    energies = []
    num_occurrences = []
    nremaining = [len(ans["energies"]) for ans in answers]
    idx_list = [0] * n_ans
    while sum(nremaining) > 0:
        # Find a minimum energy.
        min_energy = 2**30
        min_soln = []
        for i in range(n_ans):
            ans = answers[i]
            idx = idx_list[i]
            try:
                if ans["energies"][idx] < min_energy:
                    min_energy = ans["energies"][idx]
                    min_soln = ans["solutions"][idx]
            except IndexError:
                pass

        # Merge all equal solutions.
        n_occ = 0
        energies.append(min_energy)
        solutions.append(min_soln)
        for i in range(n_ans):
            ans = answers[i]
            idx = idx_list[i]
            try:
                if ans["energies"][idx] == min_energy and ans["solutions"][idx] == min_soln:
                    n_occ += ans["num_occurrences"][idx]
                    idx_list[i] += 1
                    nremaining[i] -= 1
            except IndexError:
                pass
        num_occurrences.append(n_occ)

    # Timing measurements that represent <something> per <something> are
    # averaged.  All other timing measurements are summed.
    timing = {}
    for ans in answers:
        for k, v in ans["timing"].items():
            try:
                timing[k] += v
            except KeyError:
                timing[k] = v
    for k, v in timing.items():
        if "_per_" in k:
            timing[k] = int(v/float(n_ans) + 0.5)

    # Construct a unified answer dictionary and return it.
    merged_answers = {
        "solutions": solutions,
        "energies": energies,
        "num_occurrences": num_occurrences,
        "timing": timing
    }
    return merged_answers

def submit_dwave_problem(verbosity, physical, samples, anneal_time, spin_revs, postproc):
    "Submit a QMI to the D-Wave."
    # Map abbreviated to full names for postprocessing types.
    postproc = {"none": "", "opt": "optimization", "sample": "sampling"}[postproc]

    # Determine the annealing time to use.
    if anneal_time == None:
        anneal_time = get_default_annealing_time()

    # Compute a list of the number of samples to take each iteration
    # and the number of spin reversals to perform.
    samples_list = compute_sample_counts(samples, anneal_time)
    spin_rev_list = compute_spin_rev_counts(spin_revs, samples_list)
    nqmis = len(samples_list)   # Number of (non-unique) QMIs to submit

    # Submit one or more QMIs to the D-Wave.
    problems = []
    for i in range(nqmis):
        solver_params = dict(chains=physical.embedding,
                             num_reads=samples_list[i],
                             annealing_time=anneal_time,
                             num_spin_reversal_transforms=spin_rev_list[i],
                             postprocess=postproc)
        unused_params = dict()
        while True:
            # Repeatedly remove parameters the particular solver doesn't like
            # until it actually works -- or fails for a different reason.
            try:
                weight_list = qmasm.dict_to_list(physical.weights)
                p = async_solve_ising(qmasm.solver, weight_list, physical.strengths, **solver_params)
                problems.append(p)
                break
            except ValueError as e:
                # Is there a better way to extract the failing symbol than a
                # regular expression match?
                bad_name_match = re.match(r'"(.*?)"', str(e))
                if bad_name_match == None:
                    raise e
                bad_name = bad_name_match.group(1)
                unused_params[bad_name] = solver_params[bad_name]
                del solver_params[bad_name]
            except RuntimeError as e:
                qmasm.abend(e)
    if verbosity >= 2:
        report_parameters_used(solver_params, unused_params)

    # Output problem IDs as soon as they become available.
    if verbosity >= 1:
        try:
            while any([problems[i].status()["problem_id"] == "" for i in range(nqmis)]):
                await_completion(problems, nqmis, 1)
            report_subproblems_submitted(nqmis, problems, samples_list, spin_rev_list)
        except KeyError:
            pass   # Not all solvers support "problem_id".

    # Wait for the solver to complete.
    if verbosity >= 2:
        sys.stderr.write("Number of subproblems completed:\n\n")
        cdigits = len(str(nqmis))     # Digits in the number of completed QMIs
        tdigits = len(str(nqmis*5))   # Estimate 5 seconds per QMI submission
        start_time = time.time()
    done = False
    while not done:
        done = await_completion(problems, nqmis, 10)
        if verbosity >= 2:
            ncomplete = sum([problems[i].status()["state"] == "DONE" for i in range(nqmis)])
            sys.stderr.write("    %*d of %d (%3.0f%%) after %*.0f seconds\n" %
                             (cdigits, ncomplete, nqmis,
                              100.0*float(ncomplete)/float(nqmis),
                              tdigits, time.time() - start_time))
    if verbosity >= 2:
        sys.stderr.write("\n")
    answers = [p.result() for p in problems]

    # Merge the result of seperate runs into a composite answer.
    answer = merge_answers(answers)

    # Return a Solutions object for further processing.
    return qmasm.Solutions(answer, physical, verbosity >= 2)
