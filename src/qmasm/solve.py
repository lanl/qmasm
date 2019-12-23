###################################
# Solve a binary quadratic model  #
# By Scott Pakin <pakin@lanl.gov> #
###################################

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
        sys.stderr.write("\n")
