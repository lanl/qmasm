########################################
# Store and filter solutions for QMSAM #
# By Scott Pakin <pakin@lanl.gov>      #
########################################

try:
    from dwave_sapi2.embedding import unembed_answer
except ImportError:
    from .fake_dwave import *
import qmasm

# Specify the minimum distinguishable difference between energy readings.
min_energy_delta = 0.005

# Define a class to represent a single solution.
class Solution:
    "Represent a near-minimal state of a spin system."

    def __init__(self, problem, num2syms, all_num2syms, phys2log, all_vars, soln, tally, energy):
        # Map named variables to spins.
        self.problem = problem
        self.soln_spins = soln
        self.tally = tally
        self.energy = energy
        self.phys2log = phys2log
        self.names = []       # List of names for each named row
        self.spins = []       # Spin for each named row
        self.all_names = []   # List of names for each named row, including "$" names
        self.all_spins = []   # Spin for each named row, including "$" names
        self.id = 0           # Map from spins to an int
        self._checked_asserts = None  # Memoized result of check_assertions
        for q in range(len(soln)):
            # Add all names to all_names and all_spins.
            if all_num2syms[q] == []:
                continue
            self.all_names.append(" ".join(all_num2syms[q]))
            self.all_spins.append(soln[q])

            # Add only non-"$" names to names and spins.
            if num2syms[q] == []:
                continue
            self.names.append(" ".join(num2syms[q]))
            self.spins.append(soln[q])
            self.id = self.id*2 + (soln[q] + 1)/2

        # Additionally map the spins computed during simplification.
        for nm, s in problem.known_values.items():
            self.all_names.append(nm)
            self.all_spins.append(s)
            if "$" in nm and not all_vars:
                continue
            self.names.append(nm)
            self.spins.append(s)
            self.id = self.id*2 + (s + 1)/2

    def broken_pins(self):
        "Return True if the solution contains broken pins."
        bool2spin = [-1, +1]
        for pq, pin in self.problem.pinned:
            lq = self.phys2log[pq]
            if self.soln_spins[lq] != bool2spin[pin]:
                return True
        return False

    def broken_user_chains(self):
        "Return True if the solution contains broken user-specified chains."
        # Compare logical qubits for equal values.
        for pq1, pq2 in self.problem.chains.keys():
            lq1 = self.phys2log[pq1]
            lq2 = self.phys2log[pq2]
            if self.soln_spins[lq1] != self.soln_spins[lq2]:
                return True
        return False

    def failed_assertions(self):
        "Return True if the solution contains failed assertions."
        return not all(self.check_assertions())

    def check_assertions(self):
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
                spin = (spin + 1)//2
            else:
                spin = None
            for nm in names:
                name2bit[nm] = spin

        # Test each assertion in turn.
        results = []
        for a in self.problem.assertions:
            results.append((str(a), a.evaluate(name2bit)))
        self._checked_asserts = results
        return self._checked_asserts

class Solutions:
    "Represent all near-minimal states of a spin system."

    def __init__(self, answer, problem, all_vars):
        # Store our arguments.
        self.answer = answer
        self.problem = problem
        self.all_vars = all_vars

        # Unembed the solutions.  Fix rather than discard invalid solutions.
        raw_solns = answer["solutions"]
        fixed_solns = unembed_answer(raw_solns, problem.embedding,
                                     broken_chains="minimize_energy",
                                     h=problem.weights, j=problem.strengths)

        # Define a mapping of physical to logical qubits.
        phys2log = {}
        for lq in range(len(self.problem.embedding)):
            for pq in self.problem.embedding[lq]:
                phys2log[pq] = lq

        # Establish mappings from qubit numbers to symbols.
        max_num = qmasm.sym_map.max_number()
        num2syms = [[] for _ in range(max_num + 1)]
        all_num2syms = [[] for _ in range(max_num + 1)]
        for s, n in qmasm.sym_map.symbol_number_items():
            all_num2syms[n].append(s)
            if all_vars or "$" not in s:
                num2syms[n].append(s)

        # Construct one solution object per solution.
        self.solutions = []
        tallies = answer["num_occurrences"]
        energies = [e + problem.simple_offset for e in answer["energies"]]
        for i in range(len(fixed_solns)):
            self.solutions.append(Solution(problem, num2syms, all_num2syms,
                                           phys2log, all_vars, fixed_solns[i],
                                           tallies[i], energies[i]))

    def discard_broken_chains(self):
        "Discard solutions with broken chains."
        # Tally the number of occurrences of each solution, because we lose the
        # association of solution to tally when we discard solutions.
        tallies = {tuple(s.soln_spins): s.tally for s in self.solutions}

        # Re-embed the raw solution, discarding broken chains.
        raw_solns = self.answer["solutions"]
        good_solns = unembed_answer(raw_solns, self.problem.embedding,
                                    broken_chains="discard",
                                    h=self.problem.weights,
                                    j=self.problem.strengths)

        # Convert from a list to a set for rapid lookup.
        good_soln_set = set()
        for s in good_solns:
            good_soln_set.add(tuple(s))

        # Filter the solution objects, retaining only those appearing in
        # good_solns.
        solutions = []
        for s in self.solutions:
            if tuple(s.soln_spins) in good_soln_set:
                solutions.append(s)
        self.solutions = solutions

    def discard_broken_pins(self):
        "Discard solutions with broken pins."
        self.solutions = [s for s in self.solutions if not s.broken_pins()]

    def discard_broken_user_chains(self):
        "Discard solutions with broken user-specified chains."
        self.solutions = [s for s in self.solutions if not s.broken_user_chains()]
    def discard_failed_assertions(self):
        "Discard solutions with failed assertions."
        self.solutions = [s for s in self.solutions if not s.failed_assertions()]

    def discard_non_minimal(self):
        "Discard all solutions with non-minimal energy."
        if len(self.solutions) == 0:
            return
        min_e = self.solutions[0].energy  # Assume solutions are sorted by energy.
        self.solutions = [s for s in self.solutions if abs(s.energy - min_e) < min_energy_delta]

    def discard_duplicates(self):
        "Discard all duplicate solutions (same assignments to all non-ignored variables)."
        id2soln = {}
        for s in self.solutions:
            if s.id in id2soln:
                id2soln[s.id].tally += s.tally  # Accumulate the tallies of merged solutions.
            else:
                id2soln[s.id] = s
            self.solutions = id2soln.values()
        self.solutions.sort(key=lambda s: s.energy)
