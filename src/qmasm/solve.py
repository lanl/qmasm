###################################
# Solve a binary quadratic model  #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from dimod import ExactSolver
from neal import SimulatedAnnealingSampler
from tabu import TabuSampler
from dwave.cloud import Client
from dwave.system import DWaveSampler

class Sampler(object):
    "Interface to Ocean samplers."

    def __init__(self, qmasm, profile=None, solver=None):
        "Acquire either a software sampler or a sampler representing a hardware solver."
        self.qmasm = qmasm
        self.sampler = self._get_sampler(profile, solver)

        # Handle built-in software samplers as special cases.
        self.sampler = None
        if solver == "exact":
            self.sampler = ExactSampler()
        elif solver == "sim-anneal":
            self.sampler = SimulatedAnnealingSampler()
        elif solver == "tabu":
            self.sampler = TabuSampler()
        else:
            # In the common case, read the configuration file, either the
            # default or the one named by the DWAVE_CONFIG_FILE environment
            # variable.
            client = Client.from_config(profile=profile)
            if solver == None:
                solver = client.default_solver
            else:
                solver = {"name": solver}
            self.sampler = DWaveSampler(solver=solver)
            client.close()

    def _get_sampler(self, profile, solver):
        "Return a dimod.Sampler object."
        # Handle built-in software samplers as special cases.
        if solver == "exact":
            return ExactSampler()
        if solver == "sim-anneal":
            return SimulatedAnnealingSampler()
        if solver == "tabu":
            return TabuSampler()

        # In the common case, read the configuration file, either the
        # default or the one named by the DWAVE_CONFIG_FILE environment
        # variable.
        try:
            with Client.from_config(profile=profile) as client:
                if solver == None:
                    solver = client.default_solver
                else:
                    solver = {"name": solver}
                    return DWaveSampler(solver=solver)
        except Exception as err:
            self.qmasm.abend("Failed to construct a sampler (%s)" % str(err))
