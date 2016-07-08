###################################
# Parse a QASM source file        #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import os
import qasm
import shlex
import string
import sys

# Define a function that aborts the program, reporting an invalid
# input line as part of the error message.
filename = "<stdin>"
lineno = 0
def error_in_line(str):
    sys.stderr.write('%s:%d: error: %s\n' % (filename, lineno, str))
    sys.exit(1)

# Define synonyms for "true" and "false".
str2bool = {s: True for s in ["1", "+1", "T", "TRUE"]}
str2bool.update({s: False for s in ["0", "-1", "F", "FALSE"]})

# Define a function that searches a list of directories for a file.
def find_file_in_path(pathnames, filename):
    for pname in pathnames:
        fname = os.path.join(pname, filename)
        if os.path.exists(fname):
            return fname
        fname_qasm = fname + ".qasm"
        if os.path.exists(fname_qasm):
            return fname_qasm
    return None

# Define a function that says if a string can be treated as a float.
def is_float(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

class Statement(object):
    "One statement in a QASM source file."

    def __init__(self, lineno):
        self.lineno = lineno

    def error_in_line(self, msg):
        if self.lineno == None:
            qasm.abend(msg)
        else:
            sys.stderr.write('%s:%d: error: %s\n' % (cl_args.input, self.lineno, msg))
        sys.exit(1)

class Weight(Statement):
    "Represent a point weight on a qubit."
    def __init__(self, lineno, sym, weight):
        super(Weight, self).__init__(lineno)
        self.sym = sym
        self.weight = weight

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        num = qasm.symbol_to_number(prefix + self.sym)
        weights[num] += self.weight

class Chain(Statement):
    "Chain between qubits."
    def __init__(self, lineno, sym1, sym2):
        super(Chain, self).__init__(lineno)
        self.sym1 = sym1
        self.sym2 = sym2

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        num1 = qasm.symbol_to_number(prefix + self.sym1)
        num2 = qasm.symbol_to_number(prefix + self.sym2)
        if num1 == num2:
            self.error_in_line("A chain cannot connect a spin to itself")
        elif num1 > num2:
            num1, num2 = num2, num1
        chains[(num1, num2)] = None   # Value is a don't-care.

class Pin(Statement):
    "Pinning of a qubit to true or false."
    def __init__(self, lineno, sym, goal):
        super(Pin, self).__init__(lineno)
        self.sym = sym
        self.goal = goal

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        num = qasm.symbol_to_number(prefix + self.sym)
        pinned.append((num, self.goal))

class Alias(Statement):
    "Alias of one symbol to another."
    def __init__(self, lineno, sym1, sym2):
        super(Alias, self).__init__(lineno)
        self.sym1 = sym1
        self.sym2 = sym2

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        sym1 = prefix + self.sym1
        sym2 = prefix + self.sym2
        try:
            qasm.sym2num[sym1] = qasm.sym2num[sym2]
        except KeyError:
            self.error_in_line("Cannot make symbol %s an alias of undefined symbol %s" % (sym1, sym2))
        if sym1 == sym2:
            self.error_in_line("Fields cannot alias themselves")

class Strength(Statement):
    "Coupler strength between two qubits."
    def __init__(self, lineno, sym1, sym2, strength):
        super(Strength, self).__init__(lineno)
        self.sym1 = sym1
        self.sym2 = sym2
        self.strength = strength

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        num1 = qasm.symbol_to_number(prefix + self.sym1)
        num2 = qasm.symbol_to_number(prefix + self.sym2)
        if num1 == num2:
            self.error_in_line("A coupler cannot connect a spin to itself")
        elif num1 > num2:
            num1, num2 = num2, num1
        strengths[(num1, num2)] += self.strength

class MacroUse(Statement):
    "Instantiation of a macro definition."
    def __init__(self, lineno, name, body, prefix):
        super(MacroUse, self).__init__(lineno)
        self.name = name
        self.body = body
        self.prefix = prefix

    def update_qmi(self, prefix, weights, strengths, chains, pinned):
        for stmt in self.body:
            stmt.update_qmi(prefix + self.prefix, weights, strengths, chains, pinned)

# Define a function that parses an input file into an internal representation.
# This function can be called recursively (due to !include directives).
macros = {}        # Map from a macro name to a list of Statement objects
current_macro = (None, [])   # Macro currently being defined (name and statements)
aliases = {}       # Map from a symbol to its textual expansion
target = qasm.program   # Reference to either the program or the current macro
def parse_file(infilename, infile):
    global macros, current_macro, aliases, target, filename, lineno
    filename = infilename
    for line in infile:
        # Split the line into fields and apply text aliases.
        lineno += 1
        if line.strip() == "":
            continue
        fields = shlex.split(line, True)
        for i in range(len(fields)):
            try:
                fields[i] = aliases[fields[i]]
            except KeyError:
                pass

        # Process the line.
        if len(fields) == 0:
            # Ignore empty lines.
            continue
        elif len(fields) >= 2 and fields[0] == "!include":
            # "!include" "<filename>" -- process a named auxiliary file.
            incname = string.join(fields[1:], " ")
            if len(incname) >= 2 and incname[0] == "<" and incname[-1] == ">":
                # Search QASMPATH for the filename.
                incname = incname[1:-1]
                try:
                    qasmpath = string.split(os.environ["QASMPATH"], ":")
                    qasmpath.append(".")
                except KeyError:
                    qasmpath = ["."]
                found_incname = find_file_in_path(qasmpath, incname)
                if found_incname != None:
                    incname = found_incname
            elif len(incname) >= 2:
                # Search only the current directory for the filename.
                found_incname = find_file_in_path(["."], incname)
                if found_incname != None:
                    incname = found_incname
            try:
                incfile = open(incname)
            except IOError:
                qasm.abend('Failed to open %s for input' % incname)
            parse_file(incname, incfile)
            incfile.close()
        elif len(fields) == 2:
            if fields[0] == "!begin_macro":
                # "!begin_macro" <name> -- begin a macro definition.
                name = fields[1]
                if macros.has_key(name):
                    error_in_line("Macro %s is multiply defined" % name)
                if current_macro[0] != None:
                    error_in_line("Nested macros are not supported")
                current_macro = (name, [])
                target = current_macro[1]
            elif fields[0] == "!end_macro":
                # "!end_macro" <name> -- end a macro definition.
                name = fields[1]
                if current_macro[0] == None:
                    error_in_line("Ended macro %s with no corresponding begin" % name)
                if current_macro[0] != name:
                    error_in_line("Ended macro %s after beginning macro %s" % (name, current_macro[0]))
                macros[name] = current_macro[1]
                target = qasm.program
                current_macro = (None, [])
            else:
                # <symbol> <weight> -- increment a symbol's point weight.
                try:
                    val = float(fields[1])
                except ValueError:
                    error_in_line('Failed to parse "%s %s" as a symbol followed by a numerical weight' % (fields[0], fields[1]))
                target.append(Weight(lineno, fields[0], val))
        elif len(fields) == 3:
            if fields[1] == "=":
                # <symbol_1> = <symbol_2> -- create a chain between <symbol_1> and
                # <symbol_2>.
                target.append(Chain(lineno, fields[0], fields[2]))
            elif fields[1] == ":=":
                # <symbol> := <value> -- force symbol <symbol> to have value
                # <value>.
                try:
                    goal = str2bool[fields[2].upper()]
                except KeyError:
                    error_in_line('Right-hand side ("%s") must be a Boolean value' % fields[2])
                target.append(Pin(lineno, fields[0], goal))
            elif fields[1] == "<->":
                # <symbol_1> <-> <symbol_2> -- make <symbol_1> an alias of
                # <symbol_2>.
                target.append(Alias(lineno, fields[0], fields[2]))
            elif fields[0] == "!use_macro":
                # "!use_macro" <macro_name> <instance_name> -- instantiate
                # a macro using <instance_name> as each variable's prefix.
                name = fields[1]
                try:
                    target.append(MacroUse(lineno, name, macros[name], fields[2] + "."))
                except KeyError:
                    error_in_line("Unknown macro %s" % name)
            elif fields[0] == "!alias":
                # "!alias" <symbol> <text> -- replace a field of <symbol> with
                # <text>.
                aliases[fields[1]] = fields[2]
            elif is_float(fields[2]):
                # <symbol_1> <symbol_2> <strength> -- increment a coupler strength.
                try:
                    strength = float(fields[2])
                except ValueError:
                    error_in_line('Failed to parse "%s" as a number' % fields[2])
                target.append(Strength(lineno, fields[0], fields[1], strength))
            else:
                # Three fields but none of the above cases
                error_in_line('Cannot parse "%s"' % line)
        else:
            # Neither two nor three fields
            error_in_line('Cannot parse "%s"' % line)

def parse_files(file_list):
    "Parse a list of file(s) into an internal representation."
    if file_list == []:
        # No files were specified: Read from standard input.
        parse_file("<stdin>", sys.stdin)
        if current_macro[0] != None:
            error_in_line("Unterminated definition of macro %s" % current_macro[0])
    else:
        # Files were specified: Process each in turn.
        for infilename in file_list:
            try:
                infile = open(infilename)
            except IOError:
                qasm.abend('Failed to open %s for input' % infilename)
            parse_file(infilename, infile)
            if current_macro[0] != None:
                error_in_line("Unterminated definition of macro %s" % current_macro[0])
            infile.close()

def parse_pin(pin):
    "Parse a pin statement passed on the command line."
    for pstr in pin:
        lhs_rhs = pstr.split(":=")
        if len(lhs_rhs) != 2:
            qasm.abend('Failed to parse --pin="%s"' % pstr)
        lhs = lhs_rhs[0].split()
        rhs = []
        for r in lhs_rhs[1].upper().split():
            try:
                rhs.append(str2bool[r])
            except KeyError:
                for subr in r:
                    try:
                        rhs.append(str2bool[subr])
                    except KeyError:
                        qasm.abend('Failed to parse --pin="%s"' % pstr)
            if len(lhs) != len(rhs):
                qasm.abend('Different number of left- and right-hand-side values in --pin="%s" (%d vs. %d)' % (pstr, len(lhs), len(rhs)))
            for l, r in zip(lhs, rhs):
                qasm.program.append(Pin(None, l, r))
