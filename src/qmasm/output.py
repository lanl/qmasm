###################################
# Output QUBOs in various formats #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import sys

class OutputMixin(object):
    "Provide functions for outputting problems and solutions."

    def open_output_file(self, oname):
        "Open a file or standard output."
        if oname == "<stdout>":
            outfile = sys.stdout
        else:
            try:
                outfile = open(oname, "w")
            except IOError:
                self.abend('Failed to open %s for output' % oname)
        return outfile

    def output_qmasm(self, outfile):
        "Output weights and strengths as a flattened QMASM source file."
        for p in self.program:
            outfile.write("%s\n" % p.as_str())

    def write_output(self, problem, oname, oformat, as_qubo):
        "Write an output file in one of a variety of formats."

        # Open the output file.
        outfile = self.open_output_file(oname)

        # Output the weights and strengths in the specified format.
        if oformat == "qubist":
            self.output_qubist(outfile, as_qubo, problem)
        elif oformat == "dw":
            self.output_dw(outfile, problem)
        elif oformat == "qbsolv":
            self.output_qbsolv(outfile, problem)
        elif oformat == "qmasm":
            self.output_qmasm(outfile)
        elif oformat == "minizinc":
            self.output_minizinc(outfile, problem)
        elif oformat == "bqpjson":
            self.output_bqpjson(outfile, as_qubo, problem)

        # Close the output file.
        if oname != "<stdout>":
            outfile.close()
