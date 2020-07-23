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

        # Define a scale factor for converting floats to ints for MiniZinc's
        # sake.
        self.minizinc_scale_factor = 10000.0

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
        logical = Problem(self)
        for stmt in self.program:
            stmt.update_qmi("", "<ERROR>", logical)

        # Define a strength for each user-specified chain and anti-chain, and
        # assign strengths to those chains.
        self.chain_strength = logical.assign_chain_strength(cl_args.chain_strength)
        if cl_args.verbose >= 1:
            sys.stderr.write("Chain strength: %7.4f\n\n" % self.chain_strength)

        # We now have enough information to produce an Ocean BinaryQuadraticModel.
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

        # Determine when to write an output file: either pre- or post-embedding.
        write_time = None
        if write_output_file:
            # Start with always_embed --> post, otherwise pre.
            write_time = "pre"
            if cl_args.always_embed:
                write_time = "post"

            # QMASM output is always pre-embedding.
            if cl_args.format == "qmasm" and write_time == "post":
                self.warn("Ignoring --always embed; incompatible with --format=qmasm")
                write_time = "pre"

            # Qubist output is always post-embedding.
            if cl_args.format == "qubist":
                write_time = "post"

        # Produce pre-embedding output files.
        if write_output_file and write_time == "pre":
            self.write_output(logical, cl_args.output, cl_args.format, cl_args.qubo, sampler)
            if not cl_args.run:
                sys.exit(0)


        # Embed the problem on the physical topology.
        physical = sampler.embed_problem(logical, cl_args.topology_file,
                                         cl_args.pack_qubits, cl_args.physical,
                                         cl_args.verbose)

        # Map each logical qubit to one or more symbols.
        max_num = self.sym_map.max_number()
        num2syms = [[] for _ in range(max_num + 1)]
        all_num2syms = [[] for _ in range(max_num + 1)]
        max_sym_name_len = 7
        for s, n in self.sym_map.symbol_number_items():
            all_num2syms[n].append(s)
            if cl_args.verbose >= 2 or "$" not in s:
                num2syms[n].append(s)
                max_sym_name_len = max(max_sym_name_len, len(repr(num2syms[n])) - 1)

        # Output the embedding.
        physical.output_embedding(cl_args.verbose, max_sym_name_len, all_num2syms)

        # Abort if any variables failed to embed.
        if physical.embedding != {}:
            danglies = physical.dangling_variables(num2syms)
            if len(danglies) > 0:
                self.abend("Disconnected variables encountered: %s" % str(sorted(danglies)))

        # Output some problem statistics.
        if cl_args.verbose > 0:
            physical.output_embedding_statistics()

        # Produce post-embedding output files.
        if write_output_file and write_time == "post":
            self.write_output(physical, cl_args.output, cl_args.format, cl_args.qubo, sampler)

        # If we weren't told to run anything we can exit now.
        if not cl_args.run:
            return

        # Solve the specified problem.
        if cl_args.verbose >= 1:
            sys.stderr.write("Submitting the problem to the %s solver.\n\n" % sampler.client_info["solver_name"])
        composites = self.parse_composite_string(cl_args.composites)
        anneal_sched = self.parse_anneal_sched_string(cl_args.schedule)
        solutions = sampler.acquire_samples(cl_args.verbose,
                                            composites,
                                            physical,
                                            anneal_sched,
                                            cl_args.samples,
                                            cl_args.anneal_time,
                                            cl_args.spin_revs,
                                            cl_args.postproc)
        solutions.report_timing_information(cl_args.verbose)
        solutions.report_chain_break_information(cl_args.verbose)

        # Filter the solutions as directed by the user.
        filtered_solns = solutions.filter(cl_args.show, cl_args.verbose, cl_args.samples)

        # Output energy tallies.
        if cl_args.verbose >= 2:
            solutions.report_energy_tallies(filtered_solns, True)
        elif cl_args.verbose >= 1:
            solutions.report_energy_tallies(filtered_solns, False)

        # Output the solution to the standard output device.
        show_asserts = cl_args.verbose >= 2 or cl_args.show in ["best", "all"]
        filtered_solns.output_solutions(cl_args.values, cl_args.verbose, show_asserts)

def main():
    "Run QMASM."
    q = QMASM()
    q.run()

if __name__ == "__main__":
    main()
