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

    def __init__(self, problem, sym2col, num2col, raw_ss, fixed_soln, tally, energy, all_vars):
        # Map named variables to spins.
        self.problem = problem
        self.sym2col = sym2col
        self.num2col = num2col
        self.raw_ss = raw_ss
        self.fixed_soln = fixed_soln
        self.tally = tally
        self.energy = energy
        self.id = 0                   # Map from spins to an int
        self._checked_asserts = None  # Memoized result of check_assertions

        # Compute an ID for the solution.
        for s, c in sorted(self.sym2col.items()):
            if "$" in s and not all_vars:
                continue
            self.id = self.id*2 + int((self.fixed_soln[c] + 1)//2)

        # Ensure every symbol has an associated value.
        self.sym2bool = self._all_symbol_values()

    def _all_symbol_values(self):
        """Return a mapping from every symbol to a value, including symbols
        with no physical representation."""
        # Assign values to all symbols that were eventually embedded.
        sym2bool = {}
        sym_map = self.problem.qmasm.sym_map
        all_sym2num = sym_map.symbol_number_items()
        for s, n in all_sym2num:
            try:
                c = self.num2col[n]
                sym2bool[s] = self.fixed_soln[c] == 1
            except KeyError:
                pass  # We handle non-embedded symbols below.

        # Include values inferred from roof duality.
        for s, spin in self.problem.logical.known_values.items():
            sym2bool[s] = spin == 1

        # Include variables pinned explicitly by the program or user.
        for num, bval in self.problem.logical.pinned:
            for s in sym_map.to_symbols(num):
                sym2bool[s] = bval

        # Include contracted variables (those chained to other another
        # variable's value).
        for num1, num2 in self.problem.logical.contractions.items():
            for s1 in sym_map.to_symbols(num1):
                s2 = list(sym_map.to_symbols(num2))[0]
                sym2bool[s1] = sym2bool[s2]
        return sym2bool

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
                idx1, idx2 = self.num2col[lq1], self.num2col[lq1]
            except KeyError:
                pass  # Elided logical qubit
            if self.fixed_soln[idx1] != self.fixed_soln[idx2]:
                return True
        for lq1, lq2 in self.problem.logical.antichains:
            try:
                idx1, idx2 = self.num2col[lq1], self.num2col[lq1]
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
                idx = self.num2col[pq]
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

        # Test each assertion in turn.
        results = []
        sym2bit = {s: int(bval) for s, bval in self.sym2bool.items()}
        for a in self.problem.assertions:
            results.append((str(a), a.evaluate(sym2bit)))
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

        # Construct one Solution object per solution.
        self.solutions = []
        energies = fixed_answer.record.energy
        tallies = fixed_answer.record.num_occurrences
        fixed_solns = fixed_answer.record.sample
        raw_solns = self.answer.record.sample
        for i in range(len(fixed_solns)):
            sset = self.answer.slice(i, i + 1)
            sym2col, num2col = self._map_to_column(fixed_answer,
                                                   self.problem.embedding,
                                                   self.problem.qmasm.sym_map)
            self.solutions.append(Solution(self.problem, sym2col, num2col, sset,
                                           fixed_solns[i], tallies[i], energies[i],
                                           all_vars))

        # Store the frequency of chain breaks across the entire SampleSet.
        self.chain_breaks = chain_break_frequency(self.answer, self.problem.embedding)

    def _map_to_column(self, sset, embedding, sym_map):
        """Return a mapping from a symbol name and from a logical qubit number
        to a column number in a SampleSet."""
        # Compute a mapping from logical qubit to column number.
        num2col = {}
        for i in range(len(sset.variables)):
            num2col[sset.variables[i]] = i

        # Compute a mapping from symbol to a column number.  We need to go
        # symbol --> logical qubit --> column number.
        sym2col = {}
        for sym, num in sym_map.symbol_number_items():
            try:
                sym2col[sym] = num2col[num]
            except KeyError:
                pass
        return sym2col, num2col

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
        total_breakage = sum(self.chain_breaks.values())
        if total_breakage == 0.0:
            sys.stderr.write("    [No broken chains encountered]\n\n")
            return

        # Report only chains that have ever been broken.
        chain_breaks = []
        max_name_len = 11
        sym_map = self.problem.qmasm.sym_map
        for n, f in self.chain_breaks.items():
            ss = sym_map.to_symbols(n)
            sstr = " ".join(sorted(ss))
            chain_breaks.append((sstr, f))
            max_name_len = max(max_name_len, len(sstr))
        chain_breaks.sort()
        sys.stderr.write("    %-*s  Broken\n" % (max_name_len, "Variable(s)"))
        sys.stderr.write("    %s  -------\n" % ("-" * max_name_len))
        for sstr, f in chain_breaks:
            if f > 0.0:
                sys.stderr.write("    %-*s  %6.2f%%\n" % (max_name_len, sstr, f*100.0))
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
        if len(self.solutions) == 0:
            return []
        min_energy = self.solutions[0].energy
        return [s for s in self.solutions if s.energy == min_energy]

    def merge_duplicates(self):
        """Merge duplicate solutions (same assignments to all non-ignored
        variables).  Replace the merged solutions."""
        id2soln = {}
        for s in copy.deepcopy(self.solutions):  # Deep copy because of destructive update of tally
            try:
                id2soln[s.id].tally += s.tally
            except KeyError:
                id2soln[s.id] = s
        solutions = list(id2soln.values())
        solutions.sort(key=lambda s: (s.energy, s.id))
        return solutions

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

        # Merge duplicate solutions.
        if verbose >= 1 or show == "valid":
            valid_solns.solutions = valid_solns.merge_duplicates()
            if verbose >= 1:
                sys.stderr.write("    %*d excluding duplicate variable assignments\n" % (ndigits, len(valid_solns.solutions)))
                sys.stderr.write("\n")
        if show == "best":
            filtered_best_solns = best_solns.merge_duplicates()
            if len(filtered_best_solns) > 1:
                best_solns.solutions = filtered_best_solns
