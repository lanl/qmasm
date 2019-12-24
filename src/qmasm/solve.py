###################################
# Solve a binary quadratic model  #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import copy
import json
import minorminer
import os
import sys
from dimod import ExactSolver
from neal import SimulatedAnnealingSampler
from tabu import TabuSampler
from dwave.cloud import Client
from dwave.system import DWaveSampler

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

    def read_hardware_adjacency(self, fname, verbosity):
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

    def get_hardware_adjacency(self):
        "Return the hardware adjacency structure, if any."
        hw_adj = self._find_adjacency(self.sampler)
        if hw_adj != None:
            hw_adj = [(u, v) for u in hw_adj for v in hw_adj[u]]
        return hw_adj

    def embed_problem(self, logical, topology_file, verbosity):
        "Embed a problem on a physical topology, if necessary."
        # Create a physical Problem.
        physical = copy.deepcopy(logical)
        physical.embedding = []

        # Minor-embed the logical problem onto the hardware topology.
        if topology_file == None:
            hw_adj = self.get_hardware_adjacency()
        else:
            hw_adj = self.read_hardware_adjacency(topology_file, verbosity)
        if hw_adj == None or len(logical.bqm.quadratic) == 0:
            return physical
        if verbosity < 2:
            physical.embedding = minorminer.find_embedding(logical.bqm.quadratic, hw_adj, verbose=0)
        else:
            # minorminer's find_embedding is hard-wired to write to stdout.
            # Trick it into writing into a pipe instead.
            sys.stderr.write("Minor-embedding the logical problem onto the physical topology:\n\n")
            sepLine = "=== EMBEDDING ===\n"
            r, w = os.pipe()
            pid = os.fork()
            if pid == 0:
                # Child -- perform the embedding.
                os.close(r)
                os.dup2(w, sys.stdout.fileno())
                embedding = minorminer.find_embedding(logical.bqm.quadratic, hw_adj, verbose=1)
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

                # Receive the embedding from the child.
                physical.embedding = json.loads(pipe.readline())
                sys.stderr.write("\n")
        return physical
