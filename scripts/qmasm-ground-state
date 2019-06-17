#! /usr/bin/env python

############################################
# Output the ground state of a QMASM macro #
# By Scott Pakin <pakin@lanl.gov>          #
############################################

import argparse
import multiprocessing
import qmasm
import threading

# Parse the command line.
cl_parser = argparse.ArgumentParser(description="Compute the ground state of a QMASM macro")
cl_parser.add_argument("input", nargs="*", default=[],
                           help="file from which to read a symbolic Hamiltonian")
cl_parser.add_argument("-m", "--macro", metavar="MACRO", default="",
                       help="name of a macro whose ground state should be reported")
cl_parser.add_argument("-a", "--all", action="store_true",
                       help="output all states, not just the ground state")
cl_parser.add_argument("-p", "--precision", type=float, default=0.005,
                       help="minimum difference between two floating-point values to be considered different")
cl_args = cl_parser.parse_args()
if cl_args.macro == "":
    qmasm.abend("A macro must be specified with --macro")
prec = cl_args.precision

# Parse the original input file(s) into an internal representation.
fparse = qmasm.FileParser()
fparse.process_files(cl_args.input)
if cl_args.macro not in fparse.macros:
    qmasm.abend('Macro "%s" not found' % cl_args.macro)

def macro_to_coeffs(macro):
    "Convert a macro to h and J coefficient maps."
    h = {}
    J = {}
    syms = set()
    for m in macro:
        if m.__class__ == qmasm.Weight:
            try:
                h[m.sym] += m.weight
            except KeyError:
                h[m.sym] = m.weight
            syms.add(m.sym)
        elif m.__class__ == qmasm.Strength:
            try:
                J[(m.sym1, m.sym2)] += m.strength
            except KeyError:
                J[(m.sym1, m.sym2)] = m.strength
            try:
                J[(m.sym2, m.sym1)] += m.strength
            except KeyError:
                J[(m.sym2, m.sym1)] = m.strength
            syms.add(m.sym1)
            syms.add(m.sym2)
        elif m.__class__ == qmasm.Assert:
            pass
        else:
            qmasm.abend("Only weights and strengths are currently supported")
    return sorted(syms, key=lambda s: ("$" in s, s)), h, J

def similar(a, b):
    "Return True if two floating-point numbers are nearly equal."
    return abs(a - b) < prec

class FindGroundState(threading.Thread):
    "FindGroundState uses a team of threads to find a truth table's ground states."
    # Define state that is shared among threads.
    class SharedState(object):
        def __init__(self):
            self.min_energy = 2**30    # Minimum energy observed
            self.all_energies = set()  # All energies, rounded to the given precision
            self.table = []            # List of {spins, energy} tuples
            self.lock = threading.Lock()

        def update(self, min_energy, all_energies, table):
            "Update all of our state at once."
            self.lock.acquire()
            self.all_energies.update(all_energies)
            if min_energy < self.min_energy:
                if not cl_args.all and not similar(self.min_energy, min_energy):
                    # New minimum -- clear the table.
                    self.table = []
                self.min_energy = min_energy
            if cl_args.all or similar(self.min_energy, min_energy):
                self.table.extend(table)
            self.lock.release()

    state = SharedState()

    def __init__(self, ncols, begin, end, step):
        threading.Thread.__init__(self)
        self.ncols = ncols
        self.begin = begin
        self.end = end
        self.step = step

    def run(self):
        "Find the ground state by exhaustively examining a range of rows."
        # In our thread, we locally examine 1/Nth of the rows.
        min_energy = 2**30
        all_energies = set()
        table = []
        ncols = self.ncols
        for row in range(self.begin, self.end, self.step):
            # Precompute a vector of spins.
            spins = [(row>>(ncols - i - 1)&1)*2 - 1 for i in range(ncols)]

            # Compute the current row's energy
            energy = 0.0
            for i in range(ncols):
                try:
                    hi = h[syms[i]]
                    energy += hi*spins[i]
                except KeyError:
                    pass
            for i in range(ncols - 1):
                for j in range(i + 1, ncols):
                    try:
                        Jij = J[(syms[i], syms[j])]
                        energy += Jij*spins[i]*spins[j]
                    except KeyError:
                        pass

            # Keep track of minimum energy, within precision limits.
            all_energies.add(round(energy/prec)*prec)
            if energy < min_energy:
                min_energy = energy
                if not cl_args.all:
                    table = []    # We don't need excited states unless we're outputting them.
            if cl_args.all or similar(energy, min_energy):
                table.append((spins, energy))

        # Merge our findings into the global state.
        self.state.update(min_energy, all_energies, table)

def output_ground_state(syms, h, J):
    "Exhaustively evaluate the ground states of a truth table."
    # Compute the ground states using multiple threads, one per
    # hardware thread.
    ncols = len(syms)
    nrows = 2**ncols
    nthreads = min(multiprocessing.cpu_count(), nrows)
    tobjs = []
    for i in range(nthreads):
        tobjs.append(FindGroundState(ncols, i, nrows, nthreads))
        tobjs[i].start()
    for t in tobjs:
        t.join()
    min_energy = FindGroundState.state.min_energy
    all_energies = FindGroundState.state.all_energies
    table = FindGroundState.state.table
    table.sort()

    # Output the ground-state rows (or all rows if verbosity is enabled).
    width = max([len(s) for s in syms])
    width = max (width, 2)  # Need room for "+1" and "-1"
    print(" ".join(["%*s" % (width, s) for s in syms]) + "     Energy GS")
    print(" ".join(["-"*width for s in syms]) + " " + "-"*10 + " --")
    for t in table:
        spins, energy = t
        spin_str = " ".join(["%+*d" % (width, s) for s in spins])
        energy_str = "%10.2f" % energy
        gs_str = " N"
        if similar(energy, min_energy):
            gs_str = " Y"
        print("%s %s %s" % (spin_str, energy_str, gs_str))
    print("")

    # Also output the gap between the ground state and first excited state.
    all_energies = sorted(set([round(e/prec)*prec for e in all_energies]))
    if len(all_energies) == 1:
        print("=== GAP: N/A ===")
    else:
        print("=== GAP: %.5g ===" % (all_energies[1] - all_energies[0]))

# Process each macro in turn.
print("=== MACRO: %s ===" % cl_args.macro)
print("")
syms, h, J = macro_to_coeffs(fparse.macros[cl_args.macro])
output_ground_state(syms, h, J)
