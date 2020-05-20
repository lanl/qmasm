###################################
# Solve a binary quadratic model  #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import copy
import hashlib
import json
import marshal
import minorminer
import os
import sys
from dimod import ExactSolver
from dwave.cloud import Client
from dwave.embedding import embed_bqm
from dwave.system import DWaveSampler
from neal import SimulatedAnnealingSampler
from tabu import TabuSampler

class EmbeddingCache(object):
    "Read and write an embedding cache file."

    def __init__(self, qmasm, edges, adj):
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
        sha.update(str(sorted(edges)).encode("utf-8"))
        sha.update(str(sorted(adj)).encode("utf-8"))
        self.hash = sha.hexdigest()

    def read(self):
        "Read an embedding from an embedding cache or None on a cache miss."
        if self.hash == None:
            return None
        try:
            h = open(os.path.join(self.cachedir, self.hash), "rb")
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
            h = open(os.path.join(self.cachedir, self.hash), "wb")
        except IOError:
            return None
        marshal.dump(embedding, h)
        h.close()

class Sampler(object):
    "Interface to Ocean samplers."

    def __init__(self, qmasm, profile=None, solver=None):
        "Acquire either a software sampler or a sampler representing a hardware solver."
        # Store our non-None parameters for later use.
        self.qmasm = qmasm
        self.sampler, self.client_info = self._get_sampler(profile, solver)

    def _get_sampler(self, profile, solver):
        "Return a dimod.Sampler object."
        # Handle built-in software samplers as special cases.
        info = {}
        if solver != None:
            info["solver_name"] = solver
        if solver == "exact":
            return ExactSolver(), info
        if solver == "neal":
            return SimulatedAnnealingSampler(), info
        if solver == "tabu":
            return TabuSampler(), info

        # In the common case, read the configuration file, either the
        # default or the one named by the DWAVE_CONFIG_FILE environment
        # variable.
        if profile != None:
            info["profile"] = profile
        try:
            with Client.from_config(profile=profile) as client:
                if solver == None:
                    solver = client.default_solver
                else:
                    solver = {"name": solver}
                sampler = DWaveSampler(profile=profile, solver=solver)
                info = {"solver_name": sampler.solver.name,
                        "endpoint": client.endpoint}
                return sampler, info
        except Exception as err:
            self.qmasm.abend("Failed to construct a sampler (%s)" % str(err))

    def show_properties(self, verbose):
        "Output either short or all solver properties."
        if verbose == 0:
            return

        # Combine properties from various sources into a single dictionary.
        props = self.client_info.copy()
        props.update(self.sampler.properties)

        # Determine the width of the widest key.
        max_key_len = len("Parameter")
        prop_keys = sorted(props)
        for k in prop_keys:
            max_key_len = max(max_key_len, len(k))

        # Output either "short" values (if verbose = 1) or all values (if
        # verbose > 1).
        short_value_len = 70 - max_key_len
        sys.stderr.write("Encountered the following solver properties:\n\n")
        sys.stderr.write("    %-*s  Value\n" % (max_key_len, "Parameter"))
        sys.stderr.write("    %s  %s\n" % ("-" * max_key_len, "-" * max(5, short_value_len)))
        for k in prop_keys:
            if isinstance(props[k], str):
                val_str = props[k]
            else:
                val_str = repr(props[k])
            if verbose >= 2 or len(val_str) <= short_value_len:
                sys.stderr.write("    %-*s  %s\n" % (max_key_len, k, val_str))
            elif len(val_str) > short_value_len:
                # Value is too long for a single -v.
                if isinstance(props[k], list):
                    val_str = "(%d items; use -v -v to view)" % len(props[k])
                else:
                    val_str = "(%d characters; use -v -v to view)" % len(val_str)
                sys.stderr.write("    %-*s  %s\n" % (max_key_len, k, val_str))
        sys.stderr.write("\n")

    def _find_adjacency(self, sampler):
        "Perform a depth-first search for an adjacency attribute."
        # Successful base case: The given sampler reports its adjacency
        # structure.
        try:
            return sampler.adjacency
        except AttributeError:
            pass

        # Failed base case: The given sampler has no children.
        try:
            children = sampler.children
        except AttributeError:
            return None

        # Recursive case: Search each child for its adjacency.
        for c in children:
            adj = self._find_adjacency(c)
            if adj != None:
                return adj
        return None

    def get_hardware_adjacency(self):
        "Return the hardware adjacency structure, if any."
        return self._find_adjacency(self.sampler)

    def read_hardware_adjacency(self, fname, verbosity):
        """Read a hardware adjacency list from a file.  Each line must contain
        a space-separated pair of vertex numbers."""
        # Read from a file to a set of vertex pairs.
        adj = set()
        lineno = 0
        if verbosity >= 2:
            sys.stderr.write("Reading hardware adjacency from %s.\n" % fname)
        with open(fname) as f:
            for orig_line in f:
                # Discard comments then parse the line into exactly two vertices.
                lineno += 1
                line = orig_line.partition("#")[0]
                try:
                    verts = [int(v) for v in line.split()]
                except ValueError:
                    self.qmasm.abend('Failed to parse line %d of file %s ("%s")' % (lineno, fname, orig_line.strip()))
                if len(verts) == 0:
                    continue
                if len(verts) != 2 or verts[0] == verts[1]:
                    self.qmasm.abend('Failed to parse line %d of file %s ("%s")' % (lineno, fname, orig_line.strip()))

                # Canonicalize the vertex numbers and add the result to the set.
                if verts[1] > verts[0]:
                    verts[0], verts[1] = verts[1], verts[0]
                adj.add((verts[0], verts[1]))
        if verbosity >= 2:
            sys.stderr.write("%d unique edges found\n\n" % len(adj))

        # Convert from vertex pairs to a dictionary from a vertex to its
        # neighbors.  Note that we treat all edges as bidirectional.
        adj_dict = {}
        for u, v in adj:
            try:
                adj_dict[u].append(v)
            except KeyError:
                adj_dict[u] = [v]
            try:
                adj_dict[v].append(u)
            except KeyError:
                adj_dict[v] = [u]
        return adj_dict

    def _find_embedding(self, edges, adj, **kwargs):
        "Wrap minorminer.find_embedding with a version that intercepts its output."
        # Minorminer accepts the hardware adjacency as a list of
        # pairs, not a map from each node to its neighbors.
        mm_adj = [(u, v) for u in adj for v in adj[u]]

        # Given verbose=0, invoke minorminer directly.
        if "verbose" not in kwargs or kwargs["verbose"] == 0:
            # Convert all keys to integers for consistency.
            embedding = minorminer.find_embedding(edges, mm_adj, **kwargs)
            return {int(k): v for k, v in embedding.items()}

        # minorminer's find_embedding is hard-wired to write to stdout.
        # Trick it into writing into a pipe instead.
        sepLine = "=== EMBEDDING ===\n"
        r, w = os.pipe()
        pid = os.fork()
        if pid == 0:
            # Child -- perform the embedding.
            os.close(r)
            os.dup2(w, sys.stdout.fileno())
            embedding = minorminer.find_embedding(edges, mm_adj, **kwargs)
            sys.stdout.flush()
            os.write(w, sepLine.encode())
            os.write(w, (json.dumps(embedding) + "\n").encode())
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
                        self.qmasm.abend("Embedder failed to terminate properly")
                    sys.stderr.write("      %s" % rstr)
                except:
                    pass

            # Receive the embedding from the child.  Convert all keys to
            # integers for consistency.
            embedding = json.loads(pipe.readline())
            return {int(k): v for k, v in embedding.items()}

    def find_problem_embedding(self, logical, topology_file, verbosity):
        """Find an embedding of a problem on a physical topology, if
        necessary.  Return a physical Sampler object."""
        # Create a physical Problem.
        physical = copy.deepcopy(logical)
        physical.embedding = {}

        # Acquire the hardware topology unless we were given a specific
        # topology to use.
        if topology_file == None:
            hw_adj = self.get_hardware_adjacency()
        else:
            hw_adj = self.read_hardware_adjacency(topology_file, verbosity)
        physical.hw_adj = hw_adj

        # See if we already have an embedding in the embedding cache.
        edges = logical.bqm.quadratic
        if hw_adj == None or len(edges) == 0:
            # Either the sampler does not require embedding, or we have no work
            # to do.
            return physical
        if verbosity >= 2:
            sys.stderr.write("Minor-embedding the logical problem onto the physical topology:\n\n")
        ec = EmbeddingCache(self.qmasm, edges, hw_adj)
        if verbosity >= 2:
            if ec.cachedir == None:
                sys.stderr.write("  No embedding cache directory was specified ($QMASMCACHE).\n")
            else:
                sys.stderr.write("  Using %s as the embedding cache directory.\n" % ec.cachedir)
        embedding = ec.read()
        if embedding == {}:
            # Cache hit, but embedding had failed
            if verbosity >= 2:
                sys.stderr.write("  Found failed embedding %s in the embedding cache.\n\n" % ec.hash)
        elif embedding != None:
            # Successful cache hit!
            if verbosity >= 2:
                sys.stderr.write("  Found successful embedding %s in the embedding cache.\n\n" % ec.hash)
            physical.embedding = embedding
            return physical
        if verbosity >= 2 and ec.cachedir != None:
            sys.stderr.write("  No existing embedding found in the embedding cache.\n")

        # Minor-embed the logical problem onto the hardware topology.
        if verbosity < 2:
            physical.embedding = self._find_embedding(edges, hw_adj)
        else:
            sys.stderr.write("  Running the embedder.\n\n")
            physical.embedding = self._find_embedding(edges, hw_adj, verbose=1)
            sys.stderr.write("\n")

        # Cache the embedding for next time.
        ec.write(physical.embedding)
        if verbosity >= 2 and ec.cachedir != None:
            sys.stderr.write("  Caching the embedding as %s.\n\n" % ec.hash)
        return physical

    def embed_problem(self, logical, topology_file, verbosity):
        "Embed a problem on a physical topology, if necessary."
        # Embed the problem.
        physical = self.find_problem_embedding(logical, topology_file, verbosity)
        physical.bqm = embed_bqm(physical.bqm, physical.embedding,
                                 physical.hw_adj, -physical.qmasm.chain_strength)

        # Update weights and strengths.  Maintain a reference to the
        # logical problem.
        physical.logical = logical
        physical.weights = physical.bqm.linear
        physical.strengths = physical.bqm.quadratic

        # Some problem parameters are not relevant in the physical problem.
        physical.chains = None
        physical.antichains = None
        physical.pinned = None
        physical.known_values = None
        return physical
