#! /usr/bin/env python

###################################
# Quantum Macro Assembler         #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import re
import sys
from qmasm.cmdline import ParseCommandLine
from qmasm.parse import FileParser
from qmasm.problem import Problem
from qmasm.utils import Utilities, SymbolMapping
from qmasm.solve import Sampler

class QMASM(ParseCommandLine, Utilities):
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

        # We now have enough information to produce an Ocean
        # BinaryQuadraticModel.
        logical.generate_bqm()

        # Convert chains to aliases where possible.
        if cl_args.O >= 1:
            logical.convert_chains_to_aliases(cl_args.verbose)

        # Simplify the problem if possible.
        if cl_args.O >= 1:
            logical.simplify_problem(cl_args.verbose)

        # Establish a connection to a D-Wave or software sampler.
        sampler = Sampler(self, profile=cl_args.profile, solver=cl_args.solver)
        sampler.show_properties(cl_args.verbose)

def main():
    "Run QMASM."
    q = QMASM()
    q.run()

if __name__ == "__main__":
    main()
