########################################
# Store and filter solutions for QMASM #
# By Scott Pakin <pakin@lanl.gov>      #
########################################

import copy
import sys
from dwave.embedding import unembed_sampleset, chain_breaks, chain_break_frequency

# Define a class to represent a single solution.
class Solution:
    "Represent a near-minimal state of a spin system."

    def __init__(self, problem, num2syms, all_num2syms, phys2log, all_vars,
                 log2idx, raw_ss, fixed_soln, tally, energy):

        # Map named variables to spins.
        self.problem = problem
        self.phys2log = phys2log
        self.log2idx = log2idx
        self.raw_ss = raw_ss
        self.fixed_soln = fixed_soln
        self.tally = tally
        self.energy = energy
        self.names = []       # List of names for each spin
        self.spins = []       # Spin for each named variable
        self.all_names = []   # List of names for each spin, including "$" names
        self.all_spins = []   # Spin for each named variable, including "$" name
        self.id = 0           # Map from spins to an int
        self._checked_asserts = None  # Memoized result of check_assertions

        # Populate various lists describing aspects of the problem/solution.
        for q in range(len(fixed_soln)):
            # Add all names to all_names and all_spins.
            if all_num2syms[q] == []:
                continue
            self.all_names.append(" ".join(all_num2syms[q]))
            self.all_spins.append(fixed_soln[q])

            # Add only non-"$" names to names and spins.
            if num2syms[q] == []:
                continue
            self.names.append(" ".join(num2syms[q]))
            self.spins.append(fixed_soln[q])
            self.id = self.id*2 + (fixed_soln[q] + 1)/2

        # Additionally map the spins computed during simplification.
        for nm, s in problem.logical.known_values.items():
            self.all_names.append(nm)
            self.all_spins.append(s)
            if "$" in nm and not all_vars:
                continue
            self.names.append(nm)
            self.spins.append(s)
            self.id = self.id*2 + (s + 1)/2

    def broken_chains(self):
        "Return True if the solution contains broken chains."
        unbroken = unembed_sampleset(self.raw_ss, self.problem.embedding,
                                     self.problem.logical.bqm,
                                     chain_break_method=chain_breaks.discard)
        return len(unbroken) == 0

    def broken_user_chains(self):
        "Return True if the solution contains broken user-specified chains or anti-chains."
        # Compare logical qubits for equal values.
        for lq1, lq2 in self.problem.logical.chains:
            try:
                idx1, idx2 = self.log2idx[lq1], self.log2idx[lq1]
            except KeyError:
                pass  # Elided logical qubit
            if self.fixed_soln[idx1] != self.fixed_soln[idx2]:
                return True
        for lq1, lq2 in self.problem.logical.antichains:
            try:
                idx1, idx2 = self.log2idx[lq1], self.log2idx[lq1]
            except KeyError:
                pass  # Elided logical qubit
            if self.fixed_soln[idx1] == self.fixed_soln[idx2]:
                return True
        return False

    def broken_pins(self):
        "Return True if the solution contains broken pins."
        bool2spin = [-1, +1]
        for pq, pin in self.problem.logical.pinned:
            try:
                idx = self.log2idx[pq]
                if self.fixed_soln[idx] != bool2spin[pin]:
                    return True
            except KeyError:
                pass  # Elided logical qubit
        return False

    def _check_assertions(self, stop_on_fail=False):
        "Return the result of applying each assertion."
        # Return the previous result, if any.
        if self._checked_asserts != None:
            return self._checked_asserts

        # Construct a mapping from names to bits.
        name2bit = {}
        for i in range(len(self.all_spins)):
            names = self.all_names[i].split()
            spin = self.all_spins[i]
            if spin in [-1, 1]:
                bit = (spin + 1)//2
            else:
                bit = None
            for nm in names:
                name2bit[nm] = bit
        for name, spin in self.problem.logical.known_values.items():
            name2bit[name] = (spin + 1)//2

        # Test each assertion in turn.
        results = []
        for a in self.problem.assertions:
            results.append((str(a), a.evaluate(name2bit)))
            if stop_on_fail and not results[-1][1]:
                return results
        self._checked_asserts = results
        return self._checked_asserts

    def failed_assertions(self, stop_on_fail):
        "Return True if the solution contains failed assertions."
        return not all([a[1] for a in self._check_assertions(stop_on_fail)])

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
            try:
                for pq in self.problem.embedding[lq]:
                    phys2log[pq] = lq
            except KeyError:
                pass  # Pinned or elided logical qubit

        # Establish a mapping from qubit numbers to symbols.
        self.qmasm = self.problem.qmasm
        max_num = self.qmasm.sym_map.max_number()
        num2syms = [[] for _ in range(max_num + 1)]
        all_num2syms = [[] for _ in range(max_num + 1)]
        for s, n in self.qmasm.sym_map.symbol_number_items():
            all_num2syms[n].append(s)
            if all_vars or "$" not in s:
                num2syms[n].append(s)

        # Establish a mapping from logical qubit number to index into a sample.
        log2idx = {}
        ss_vars = fixed_answer.variables
        for i in range(len(ss_vars)):
            log2idx[ss_vars[i]] = i

        # Construct one Solution object per solution.
        self.solutions = []
        energies = fixed_answer.record.energy
        tallies = fixed_answer.record.num_occurrences
        fixed_solns = fixed_answer.record.sample
        raw_solns = self.answer.record.sample
        for i in range(len(fixed_solns)):
            sset = self.answer.slice(i, i + 1)
            self.solutions.append(Solution(self.problem, num2syms, all_num2syms,
                                           phys2log, all_vars, log2idx, sset,
                                           fixed_solns[i], tallies[i], energies[i]))

        # Store the frequency of chain breaks across the entire SampleSet.
        cbf = chain_break_frequency(self.answer, self.problem.embedding)
        self.chain_breaks = [(num2syms[n], f) for n, f in cbf.items()]

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

    def report_chain_break_information(self, verbosity):
        "Output information about common chain breaks."
        # Ensure we have something to report.
        if verbosity < 2:
            return
        sys.stderr.write("Chain-break frequencies:\n\n")
        total_breakage = sum([cb[1] for cb in self.chain_breaks])
        if total_breakage == 0.0:
            sys.stderr.write("    [No broken chains encountered]\n\n")
            return

        # Report only chains that have ever been broken.
        chain_breaks = []
        max_name_len = 11
        for vs, f in self.chain_breaks:
            vstr = " ".join(sorted(vs))
            chain_breaks.append((vstr, f))
            max_name_len = max(max_name_len, len(vstr))
        chain_breaks.sort()
        sys.stderr.write("    %-*s  Broken\n" % (max_name_len, "Variable(s)"))
        sys.stderr.write("    %s  -------\n" % ("-" * max_name_len))
        for vs, f in chain_breaks:
            if f > 0.0:
                sys.stderr.write("    %-*s  %6.2f%%\n" % (max_name_len, vs, f*100.0))
        sys.stderr.write("\n")

    def discard_broken_chains(self):
        "Discard solutions with broken chains.  Return the remaining solutions."
        return [s for s in self.solutions if not s.broken_chains()]

    def discard_broken_user_chains(self):
        "Discard solutions with broken user-specified chains.  Return the remaining solutions."
        return [s for s in self.solutions if not s.broken_user_chains()]

    def discard_broken_pins(self):
        "Discard solutions with broken pins.  Return the remaining solutions."
        return [s for s in self.solutions if not s.broken_pins()]

    def discard_failed_assertions(self, stop_on_fail):
        "Discard solutions with failed assertions.  Return the remaining solutions."
        return [s for s in self.solutions if not s.failed_assertions(stop_on_fail)]

    def discard_non_minimal(self):
        "Discard solutions with non-minimal energy.  Return the remaining solutions."
        min_energy = self.solutions[0].energy
        return [s for s in self.solutions if s.energy == min_energy]

    def filter(self, show, verbose, nsamples):
        '''Return solutions as filtered according to the "show" parameter.
        Output information about the filtering based on the "verbose"
        parameter.'''
        # Prepare various views of the solutions.
        all_solns = copy.copy(self)
        valid_solns = copy.copy(self)
        best_solns = copy.copy(self)

        # Output the total number and number of unique solutions.
        if verbose >= 1:
            ndigits = len(str(nsamples))
            sys.stderr.write("Number of solutions found (filtered cumulatively):\n\n")
            sys.stderr.write("    %*d total solutions\n" % (ndigits, nsamples))
            sys.stderr.write("    %*d unique solutions\n" % (ndigits, len(all_solns.solutions)))

        # Filter out broken chains.
        if verbose >= 1 or show == "valid":
            valid_solns.solutions = valid_solns.discard_broken_chains()
            if verbose >= 1:
                sys.stderr.write("    %*d with no broken chains\n" % (ndigits, len(valid_solns.solutions)))
        if show == "best":
            filtered_best_solns = best_solns.discard_broken_chains()
            if len(filtered_best_solns) > 1:
                best_solns.solutions = filtered_best_solns

        # Filter out broken user-defined chains.
        if verbose >= 1 or show == "valid":
            valid_solns.solutions = valid_solns.discard_broken_user_chains()
            if verbose >= 1:
                sys.stderr.write("    %*d with no broken user-specified chains or anti-chains\n" % (ndigits, len(valid_solns.solutions)))
        if show == "best":
            filtered_best_solns = best_solns.discard_broken_user_chains()
            if len(filtered_best_solns) > 1:
                best_solns.solutions = filtered_best_solns

        # Filter out broken pins.
        if verbose >= 1 or show == "valid":
            valid_solns.solutions = valid_solns.discard_broken_pins()
            if verbose >= 1:
                sys.stderr.write("    %*d with no broken pins\n" % (ndigits, len(valid_solns.solutions)))
        if show == "best":
            filtered_best_solns = best_solns.discard_broken_pins()
            if len(filtered_best_solns) > 1:
                best_solns.solutions = filtered_best_solns

        # Filter out failed assertions.
        stop_on_fail = show == "valid" and verbose < 2
        if verbose >= 1 or show == "valid":
            valid_solns.solutions = valid_solns.discard_failed_assertions(stop_on_fail)
            if verbose >= 1:
                sys.stderr.write("    %*d with no failed assertions\n" % (ndigits, len(valid_solns.solutions)))
        if show == "best":
            filtered_best_solns = best_solns.discard_failed_assertions(stop_on_fail)
            if len(filtered_best_solns) > 1:
                best_solns.solutions = filtered_best_solns

        # Filter out solutions that are not at minimal energy.
        if verbose >= 1 or show == "valid":
            valid_solns.solutions = valid_solns.discard_non_minimal()
            if verbose >= 1:
                sys.stderr.write("    %*d at minimal energy\n" % (ndigits, len(valid_solns.solutions)))
        if show == "best":
            filtered_best_solns = best_solns.discard_non_minimal()
            if len(filtered_best_solns) > 1:
                best_solns.solutions = filtered_best_solns
