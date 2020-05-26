########################################
# Store and filter solutions for QMASM #
# By Scott Pakin <pakin@lanl.gov>      #
########################################

import sys
from dwave.embedding import unembed_sampleset, chain_breaks

# Define a class to represent a single solution.
class Solution:
    "Represent a near-minimal state of a spin system."

    def __init__(self, problem, num2syms, all_num2syms, phys2log, all_vars,
                 phys_soln, log_soln, tally, energy):
        # Map named variables to spins.
        self.problem = problem
        self.soln_spins = log_soln
        self.raw_spins = phys_soln
        self.tally = tally
        self.energy = energy
        self.phys2log = phys2log
        self.names = []       # List of names for each spin
        self.spins = []       # Spin for each named variable
        self.all_names = []   # List of names for each spin, including "$" names
        self.all_spins = []   # Spin for each named variable, including "$" name
        self.id = 0           # Map from spins to an int
        self._checked_asserts = None  # Memoized result of check_assertions

        # Populate various lists describing aspects of the problem/solution.
        for q in range(len(log_soln)):
            # Add all names to all_names and all_spins.
            if all_num2syms[q] == []:
                continue
            self.all_names.append(" ".join(all_num2syms[q]))
            self.all_spins.append(log_soln[q])

            # Add only non-"$" names to names and spins.
            if num2syms[q] == []:
                continue
            self.names.append(" ".join(num2syms[q]))
            self.spins.append(log_soln[q])
            self.id = self.id*2 + (log_soln[q] + 1)/2

        # Additionally map the spins computed during simplification.
        for nm, s in problem.logical.known_values.items():
            self.all_names.append(nm)
            self.all_spins.append(s)
            if "$" in nm and not all_vars:
                continue
            self.names.append(nm)
            self.spins.append(s)
            self.id = self.id*2 + (s + 1)/2

class Solutions(object):
    "Represent all near-minimal states of a spin system."

    def __init__(self, answer, problem, all_vars):
        # Store our arguments.
        self.answer = answer
        self.problem = problem
        self.all_vars = all_vars

        # Unembed the solutions.  Fix rather than discard invalid solutions.
        fixed_answer = unembed_sampleset(self.answer, self.problem.embedding,
                                         self.problem.logical.bqm,
                                         chain_break_method=chain_breaks.majority_vote)

        # Define a mapping of physical to logical qubits.
        phys2log = {}
        for lq in range(len(self.problem.embedding)):
            for pq in self.problem.embedding[lq]:
                phys2log[pq] = lq

        # Establish mappings from qubit numbers to symbols.
        self.qmasm = self.problem.qmasm
        max_num = self.qmasm.sym_map.max_number()
        num2syms = [[] for _ in range(max_num + 1)]
        all_num2syms = [[] for _ in range(max_num + 1)]
        for s, n in self.qmasm.sym_map.symbol_number_items():
            all_num2syms[n].append(s)
            if all_vars or "$" not in s:
                num2syms[n].append(s)

        # Construct one solution object per solution.
        self.solutions = []
        energies = fixed_answer.record.energy
        tallies = fixed_answer.record.num_occurrences
        fixed_solns = fixed_answer.record.sample
        raw_solns = self.answer.record.sample
        for i in range(len(fixed_solns)):
            self.solutions.append(Solution(self.problem, num2syms, all_num2syms,
                                           phys2log, all_vars, raw_solns[i],
                                           fixed_solns[i], tallies[i], energies[i]))

    def report_timing_information(self, verbosity):
        "Output solver timing information."
        if verbosity == 0:
            return
        timing_info = sorted(list(self.answer.info["timing"].items()))
        sys.stderr.write("Timing information:\n\n")
        sys.stderr.write("    %-30s %-10s\n" % ("Measurement", "Value (us)"))
        sys.stderr.write("    %s %s\n" % ("-" * 30, "-" * 10))
        for timing_value in timing_info:
            sys.stderr.write("    %-30s %10d\n" % timing_value)
        sys.stderr.write("\n")
