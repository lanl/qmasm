#! /usr/bin/env python

###################################
# Quantum Macro Assembler         #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import re
import sys
from qmasm.cmdline import ParseCommandLine
from qmasm.parse import FileParser
from qmasm.problem import Problem, BQMMixins
from qmasm.utils import Utilities, SymbolMapping
from qmasm.solve import Sampler

class QMASM(ParseCommandLine, Utilities, BQMMixins):
    "QMASM represents everything the program can do."

    def __init__(self):
        # List of Statement objects
        self.program = []

        # Map between symbols and numbers.
        self.sym_map = SymbolMapping()

        # Multiple components of QMASM require a definition of an identifier.
        self.ident_re = re.compile(r'[^-+*/%&\|^~()<=>#,\s]+')

    def run(self):
        "Execute the entire QMASM processing sequence."

        # Parse the command line.
        cl_args = self.parse_command_line()
        self.report_command_line(cl_args)

        # Parse the original input file(s) into an internal representation.
        fparse = FileParser(self)
        fparse.process_files(cl_args.input)

        # Parse the variable pinnings specified on the command line.  Append
        # these to the program.
        if cl_args.pin != None:
            for pin in cl_args.pin:
                self.program.extend(fparse.process_pin("[command line]", 1, pin))

        # Walk the statements in the program, processing each in turn.
        logical = Problem(cl_args.qubo)
        for stmt in self.program:
            stmt.update_qmi("", "<ERROR>", logical)

        # Define a strength for each user-specified chain and anti-chain, and
        # assign strengths to those chains.
        self.chain_strength = logical.assign_chain_strength(cl_args.chain_strength)
        if cl_args.verbose >= 1:
            sys.stderr.write("Chain strength: %7.4f\n\n" % self.chain_strength)

        # Work from now on with an Ocean BinaryQuadraticModel.
        bqm = logical.as_bqm()

        # Convert chains to aliases where possible.
        if cl_args.O >= 1:
            # Say what we're about to do
            if cl_args.verbose >= 2:
                sys.stderr.write("Replaced user-defined chains with aliases:\n\n")
                sys.stderr.write("  %6d logical qubits before optimization\n" % len(self.set_of_all_variables(bqm)))

            # Replace chains with aliases wherever we can.
            self.convert_chains_to_aliases(bqm)

            # Summarize what we just did.
            if cl_args.verbose >= 2:
                sys.stderr.write("  %6d logical qubits after optimization\n\n" % len(self.set_of_all_variables(bqm)))

        # Establish a connection to a D-Wave or software sampler.
        sampler = Sampler(self, profile=cl_args.profile, solver=cl_args.solver)

def main():
    "Run QMASM."
    q = QMASM()
    q.run()

if __name__ == "__main__":
    main()
