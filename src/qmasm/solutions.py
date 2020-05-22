########################################
# Store and filter solutions for QMASM #
# By Scott Pakin <pakin@lanl.gov>      #
########################################

import sys
from dwave.embedding import unembed_sampleset, chain_breaks

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
