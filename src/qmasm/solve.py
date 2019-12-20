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

    def __init__(self, profile=None, solver=None):
        "Acquire either a software sampler or a sampler representing a hardware solver."
        # Handle built-in software samplers as special cases.
        self.sampler = None
        if solver == "exact":
            self.sampler = ExactSampler()
        elif solver == "sim-anneal":
            self.sampler = SimulatedAnnealingSampler()
        elif solver == "tabu":
            self.sampler = TabuSampler()
        else:
            # Read the configuration file, either the default or the one named
            # by the DWAVE_CONFIG_FILE environment variable.
            client = Client.from_config(profile=profile)
            if solver == None:
                solver = client.default_solver                
            else:
                solver = {"name": solver}
            self.sampler = DWaveSampler(solver=solver)
            client.close()
