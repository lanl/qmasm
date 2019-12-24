#! /usr/bin/env python

###################################
# Quantum Macro Assembler         #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import re
import sys
from qmasm.cmdline import ParseCommandLine
from qmasm.output import OutputMixin
from qmasm.parse import FileParser
from qmasm.problem import Problem
from qmasm.solve import Sampler
from qmasm.utils import Utilities, SymbolMapping

class QMASM(ParseCommandLine, Utilities, OutputMixin):
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
        logical = Problem(self, cl_args.qubo)
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

        # Convert user-specified chains, anti-chains, and pins to assertions.
        logical.append_assertions_from_statements()

        # Determine if we're expected to write an output file.  If --run was
        # specified, we write a file only if --output was also specified.
        write_output_file = not (cl_args.output == "<stdout>" and cl_args.run)

        # If the user requested QMASM output, always output it here.
        if write_output_file and cl_args.format == "qmasm":
            self.write_output(logical, cl_args.output, cl_args.format, cl_args.qubo)
            if not cl_args.run:
                sys.exit(0)

        # If the user requested bqpjson output, output it here unless
        # --always-embed was specified.
        if write_output_file and cl_args.format == "bqpjson" and not cl_args.always_embed:
            self.write_output(logical, cl_args.output, cl_args.format, cl_args.qubo)
            if not cl_args.run:
                sys.exit(0)

        # Embed the problem on the physical topology.
        physical = sampler.embed_problem(logical, cl_args.topology_file, cl_args.verbose)

def main():
    "Run QMASM."
    q = QMASM()
    q.run()

if __name__ == "__main__":
    main()
