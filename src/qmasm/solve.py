###################################
# Solve a binary quadratic model  #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import copy
import dimod
import hashlib
import json
import marshal
import minorminer
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dimod import ExactSolver, SampleSet
from dwave.cloud import Client
from dwave.embedding import embed_bqm
from dwave.system import DWaveSampler
from neal import SimulatedAnnealingSampler
from qmasm.solutions import Solutions
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
        elif solver == "neal":
            return SimulatedAnnealingSampler(), info
        elif solver == "tabu":
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
        # Embed the problem.  We first filter out isolated variables that don't
        # appear in the embedding graph to prevent embed_bqm from complaining.
        physical = self.find_problem_embedding(logical, topology_file, verbosity)
        sub_bqm = physical.bqm.copy()
        sub_bqm.remove_variables_from([q
                                       for q in physical.bqm.linear
                                       if q not in physical.embedding and physical.bqm.linear[q] == 0.0])
        for q, wt in sub_bqm.linear.items():
            if q not in physical.embedding and wt != 0.0:
                abend("Logical qubit %d has a nonzero weight (%.5g) but was not embedded" % (q, wt))
        physical.bqm = embed_bqm(sub_bqm, physical.embedding,
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

    def _get_default_annealing_time(self):
        "Determine a suitable annealing time to use if none was specified."
        try:
            # Use the default value.
            anneal_time = self.sampler.properties["default_annealing_time"]
        except KeyError:
            try:
                # If the default value is undefined, use the minimum allowed
                # value.
                anneal_time = self.sampler.properties["annealing_time_range"][0]
            except KeyError:
                # If all else fails, use 20 as a reasonable default.
                anneal_time = 20
        return anneal_time

    def _compute_sample_counts(self, samples, anneal_time):
        "Return a list of sample counts to request of the hardware."
        # The formula in the D-Wave System Documentation (under
        # "problem_run_duration_range") is Duration = (annealing_time +
        # readout_thermalization)*num_reads + programming_thermalization.  We
        # split the number of samples so that each run time is less than the
        # maximum duration.
        props = self.sampler.properties
        try:
            max_run_duration = props["problem_run_duration_range"][1]
            prog_therm = props["default_programming_thermalization"]
            read_therm = props["default_readout_thermalization"]
        except KeyError:
            # Assume we're on a software solver.
            return [samples]
        max_samples = (max_run_duration - prog_therm)//(anneal_time + read_therm)
        # Independent of the maximum sample count just computed, the number of
        # samples can't exceed the maximum the hardware allows.
        try:
            max_samples = min(max_samples, props["num_reads_range"][1])
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

    def _compute_spin_rev_counts(self, spin_revs, samples_list):
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

    def _merge_results(self, results):
        "Merge results into a single SampleSet."
        sum_keys = ["total_real_time",
                    "qpu_access_overhead_time",
                    "post_processing_overhead_time",
                    "qpu_sampling_time",
                    "total_post_processing_time",
                    "qpu_programming_time",
                    "run_time_chip",
                    "qpu_access_time"]
        avg_keys = ["anneal_time_per_run",
                    "readout_time_per_run",
                    "qpu_delay_time_per_sample",
                    "qpu_anneal_time_per_sample",
                    "qpu_readout_time_per_sample"]
        timing = {}
        merged = dimod.concatenate(results)
        nsamples = len(merged)
        for sk in sum_keys:
            timing[sk] = sum([r.info["timing"][sk] for r in results])
        for ak in avg_keys:
            timing[ak] = sum([r.info["timing"][ak]*len(r) for r in results])//nsamples
        merged.info["timing"] = timing
        return merged

    def _submit_and_block(self, bqm, **params):
        "Submit a job and wait for it to complete"
        # This is a workaround for a bug in dwave-system.  See
        # https://github.com/dwavesystems/dwave-system/issues/297#issuecomment-632384524
        result = self.sampler.sample(bqm, **params)
        result.resolve()
        return result

    def acquire_samples(self, verbosity, physical, samples, anneal_time, spin_revs, postproc):
        "Acquire a number of samples from either a hardware or software sampler."
        # Map abbreviated to full names for postprocessing types.
        postproc = {"none": "", "opt": "optimization", "sample": "sampling"}[postproc]

        # Determine the annealing time to use.
        if anneal_time == None:
            anneal_time = self._get_default_annealing_time()

        # Compute a list of the number of samples to take each iteration
        # and the number of spin reversals to perform.
        samples_list = self._compute_sample_counts(samples, anneal_time)
        spin_rev_list = self._compute_spin_rev_counts(spin_revs, samples_list)
        nqmis = len(samples_list)   # Number of (non-unique) QMIs to submit

        # Submit all of our QMIs asynchronously.
        results = [None for i in range(nqmis)]
        executor = ThreadPoolExecutor()
        for i in range(nqmis):
            solver_params = dict(chains=list(physical.embedding.values()),
                                 num_reads=samples_list[i],
                                 annealing_time=anneal_time,
                                 num_spin_reversal_transforms=spin_rev_list[i],
                                 postprocess=postproc)
            future = executor.submit(self._submit_and_block, physical.bqm, **solver_params)
            results[i] = SampleSet.from_future(future)

        # Wait for the QMIs to finish.
        executor.shutdown(wait=False)
        if verbosity == 1:
            if nqmis == 1:
                sys.stderr.write("Waiting for the problem to complete.\n\n")
            else:
                sys.stderr.write("Waiting for all %d subproblems to complete.\n\n" % nqmis)
        elif verbosity >= 2:
            sys.stderr.write("Number of subproblems completed:\n\n")
            cdigits = len(str(nqmis))     # Digits in the number of completed QMIs
            tdigits = len(str(nqmis*5))   # Estimate 5 seconds per QMI submission
            start_time = time.time()
        ncomplete = 0
        prev_ncomplete = 0
        while ncomplete < nqmis:
            ncomplete = sum([int(r.done()) for r in results])
            if verbosity >= 2 and ncomplete > prev_ncomplete:
                end_time = time.time()
                sys.stderr.write("    %*d of %d (%3.0f%%) after %*.0f seconds\n" %
                                 (cdigits, ncomplete, nqmis,
                                  100.0*float(ncomplete)/float(nqmis),
                                  tdigits, end_time - start_time))
                prev_ncomplete = ncomplete
            if ncomplete < nqmis:
                time.sleep(1)
        if verbosity >= 2:
            sys.stderr.write("\n")
            sys.stderr.write("    Average time per subproblem: %.2g seconds\n\n" % ((end_time - start_time)/nqmis))
            sys.stderr.write("IDs of completed subproblems:\n\n")
            for i in range(nqmis):
                sys.stderr.write("    %s\n" % results[i].info["problem_id"])
            sys.stderr.write("\n")

        # Merge the result of seperate runs into a composite answer.
        answer = self._merge_results(results)

        # Return a Solutions object for further processing.
        return Solutions(answer, physical, verbosity >= 2)
