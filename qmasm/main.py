#! /usr/bin/env python

###################################
# Quantum Macro Assembler         #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import cmdline

class QMASM(cmdline.ParseCommandLine):
    "QMASM represents everything the program can do."
    
    def run(self):
        "Execute the entire QMASM processing sequence."

        # Parse the command line.
        cl_args = self.parse_command_line()
        self.report_command_line(cl_args)

# Run QMASM as a script.
if __name__== "__main__":
    q = QMASM()
    q.run()
