#! /usr/bin/env python

###################################
# Quantum Macro Assembler         #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from qmasm.cmdline import ParseCommandLine

class QMASM(ParseCommandLine):
    "QMASM represents everything the program can do."
    
    def run(self):
        "Execute the entire QMASM processing sequence."

        # Parse the command line.
        cl_args = self.parse_command_line()
        self.report_command_line(cl_args)

def main():
    "Run QMASM."
    q = QMASM()
    q.run()

if __name__ == "__main__":
    main()
