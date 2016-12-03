###################################
# Parse a QMASM source file       #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import os
import qmasm
import re
import shlex
import string
import sys

# Define a function that aborts the program, reporting an invalid
# input line as part of the error message.
filename = "<stdin>"
lineno = 0
def error_in_line(str):
    global filename, lineno
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
        fname_qmasm = fname + ".qmasm"
        if os.path.exists(fname_qmasm):
            return fname_qmasm
    return None

# Define a function that says if a string can be treated as a float.
def is_float(str):
    try:
        float(str)
        return True
    except ValueError:
        return False

class Statement(object):
    "One statement in a QMASM source file."

    def __init__(self, filename, lineno):
        self.filename = filename
        self.lineno = lineno

    def error_in_line(self, msg):
        if self.lineno == None:
            qmasm.abend(msg)
        else:
            sys.stderr.write('%s:%d: error: %s\n' % (self.filename, self.lineno, msg))
        sys.exit(1)

class Weight(Statement):
    "Represent a point weight on a qubit."
    def __init__(self, filename, lineno, sym, weight):
        super(Weight, self).__init__(filename, lineno)
        self.sym = sym
        self.weight = weight

    def as_str(self, prefix=""):
        return "%s%s %s" % (prefix, self.sym, self.weight)

    def update_qmi(self, prefix, problem):
        num = qmasm.symbol_to_number(prefix + self.sym)
        problem.weights[num] += self.weight

class Chain(Statement):
    "Chain between qubits."
    def __init__(self, filename, lineno, sym1, sym2):
        super(Chain, self).__init__(filename, lineno)
        self.sym1 = sym1
        self.sym2 = sym2

    def as_str(self, prefix=""):
        return "%s%s = %s%s" % (prefix, self.sym1, prefix, self.sym2)

    def update_qmi(self, prefix, problem):
        num1 = qmasm.symbol_to_number(prefix + self.sym1)
        num2 = qmasm.symbol_to_number(prefix + self.sym2)
        if num1 == num2:
            self.error_in_line("A chain cannot connect a spin to itself")
        elif num1 > num2:
            num1, num2 = num2, num1
        problem.chains[(num1, num2)] = None   # Value is a don't-care.

class Pin(Statement):
    "Pinning of a qubit to true or false."
    def __init__(self, filename, lineno, sym, goal):
        super(Pin, self).__init__(filename, lineno)
        self.sym = sym
        self.goal = goal

    def as_str(self, prefix=""):
        return "%s%s := %s" % (prefix, self.sym, self.goal)

    def update_qmi(self, prefix, problem):
        num = qmasm.symbol_to_number(prefix + self.sym)
        problem.pinned.append((num, self.goal))

class Alias(Statement):
    "Alias of one symbol to another."
    def __init__(self, filename, lineno, sym1, sym2):
        super(Alias, self).__init__(filename, lineno)
        self.sym1 = sym1
        self.sym2 = sym2

    def as_str(self, prefix=""):
        return "%s%s <-> %s%s" % (prefix, self.sym1, prefix, self.sym2)

    def update_qmi(self, prefix, problem):
        sym1 = prefix + self.sym1
        sym2 = prefix + self.sym2
        try:
            qmasm.sym2num[sym1] = qmasm.sym2num[sym2]
        except KeyError:
            self.error_in_line("Cannot make symbol %s an alias of undefined symbol %s" % (sym1, sym2))
        if sym1 == sym2:
            self.error_in_line("Fields cannot alias themselves")

class Strength(Statement):
    "Coupler strength between two qubits."
    def __init__(self, filename, lineno, sym1, sym2, strength):
        super(Strength, self).__init__(filename, lineno)
        self.sym1 = sym1
        self.sym2 = sym2
        self.strength = strength

    def as_str(self, prefix=""):
        return "%s%s %s%s %s" % (prefix, self.sym1, prefix, self.sym2, self.strength)

    def update_qmi(self, prefix, problem):
        num1 = qmasm.symbol_to_number(prefix + self.sym1)
        num2 = qmasm.symbol_to_number(prefix + self.sym2)
        if num1 == num2:
            self.error_in_line("A coupler cannot connect a spin to itself")
        elif num1 > num2:
            num1, num2 = num2, num1
        problem.strengths[(num1, num2)] += self.strength

class MacroUse(Statement):
    "Instantiation of a macro definition."
    def __init__(self, filename, lineno, name, body, prefix):
        super(MacroUse, self).__init__(filename, lineno)
        self.name = name
        self.body = body
        self.prefix = prefix

    def as_str(self, prefix=""):
        stmt_strs = []
        for stmt in self.body:
            stmt_strs.append(stmt.as_str(self.prefix + prefix))
        return "\n".join(stmt_strs)

    def update_qmi(self, prefix, problem):
        for stmt in self.body:
            stmt.update_qmi(prefix + self.prefix, problem)

# Define a function that parses an input file into an internal representation.
# This function can be called recursively (due to !include directives).
macros = {}        # Map from a macro name to a list of Statement objects
current_macro = (None, [])   # Macro currently being defined (name and statements)
aliases = {}       # Map from a symbol to its textual expansion
target = qmasm.program   # Reference to either the program or the current macro
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
                # Search QMASMPATH for the filename.
                incname = incname[1:-1]
                try:
                    qmasmpath = string.split(os.environ["QMASMPATH"], ":")
                    qmasmpath.append(".")
                except KeyError:
                    qmasmpath = ["."]
                found_incname = find_file_in_path(qmasmpath, incname)
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
                qmasm.abend('Failed to open %s for input' % incname)
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
                target = qmasm.program
                current_macro = (None, [])
            else:
                # <symbol> <weight> -- increment a symbol's point weight.
                try:
                    val = float(fields[1])
                except ValueError:
                    error_in_line('Failed to parse "%s %s" as a symbol followed by a numerical weight' % (fields[0], fields[1]))
                target.append(Weight(filename, lineno, fields[0], val))
        elif len(fields) == 3:
            if fields[1] == "=":
                # <symbol_1> = <symbol_2> -- create a chain between <symbol_1>
                # and <symbol_2>.
                target.extend(process_chain(filename, lineno, " ".join(fields[:3])))
            elif fields[1] == ":=":
                # <symbol> := <value> -- force symbol <symbol> to have value
                # <value>.
                target.extend(process_pin(filename, lineno, " ".join(fields[:3])))
            elif fields[1] == "<->":
                # <symbol_1> <-> <symbol_2> -- make <symbol_1> an alias of
                # <symbol_2>.
                target.extend(process_alias(filename, lineno, " ".join(fields[:3])))
            elif fields[0] == "!use_macro":
                # "!use_macro" <macro_name> <instance_name> -- instantiate
                # a macro using <instance_name> as each variable's prefix.
                name = fields[1]
                try:
                    target.append(MacroUse(filename, lineno, name, macros[name], fields[2] + "."))
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
                target.append(Strength(filename, lineno, fields[0], fields[1], strength))
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
                qmasm.abend('Failed to open %s for input' % infilename)
            parse_file(infilename, infile)
            if current_macro[0] != None:
                error_in_line("Unterminated definition of macro %s" % current_macro[0])
            infile.close()

class PinParser(object):
    "Provide methods for parsing a pin statement."

    def __init__(self):
        self.bracket_re = re.compile(r'^\s*(\d+)(\s*(?:\.\.|:)\s*(\d+))?\s*$')
        self.bool_re = re.compile(r'TRUE|FALSE|T|F|0|[-+]?1', re.IGNORECASE)

    def expand_brackets(self, vars, expr):
        """Repeat one or more variables for each bracketed expression.  For
        example, expanding ("hello", "1 .. 3") should produce
        ("hello[1]", "hello[2]", "hello[3]")."""
        # Determine the starting and ending numbers and the step.
        bmatch = self.bracket_re.search(expr)
        if bmatch == None:
            return ["%s[%s]" % (v, expr) for v in vars]
        bmatches = bmatch.groups()
        num1 = int(bmatches[0])
        if bmatches[2] == None:
            num2 = num1
        else:
            num2 = int(bmatches[2])
        if num1 <= num2:
            step = 1
        else:
            step = -1

        # Append the same bracketed constant to each variable.
        new_vars = []
        for v in vars:
            for i in range(num1, num2 + step, step):
                new_vars.append("%s[%d]" % (v, i))
        return new_vars

    def parse_lhs(self, lhs):
        "Parse the left-hand side of a pin statement."
        variables = [""]
        group_len = 1    # Number of variables produced from the same bracketed expression
        bracket_expr = ""
        in_bracket = False
        for c in lhs:
            if c == "[":
                if in_bracket:
                    qmasm.abend("Nested brackets are not allowed")
                in_bracket = True
            elif c == "]":
                if not in_bracket:
                    qmasm.abend('Encountered "]" before seeing a "["')
                old_vars = variables[:-group_len]
                current_vars = variables[-group_len:]
                new_vars = self.expand_brackets(current_vars, bracket_expr)
                variables = old_vars + new_vars
                group_len = len(new_vars)
                in_bracket = False
                bracket_expr = ""
            elif in_bracket:
                bracket_expr += c
            elif c == " " or c == "\t":
                if in_bracket:
                    qmasm.abend("Unterminated bracketed expression")
                if variables[-1] != "":
                    variables.append("")
                group_len = 1
            else:
                for i in range(1, group_len + 1):
                    variables[-i] += c
        if in_bracket:
            qmasm.abend("Unterminated bracketed expression")
        if variables[-1] == "":
            variables.pop()
        return variables

    def parse_rhs(self, rhs):
        "Parse the right-hand side of a pin statement."
        for inter in [t.strip() for t in self.bool_re.split(rhs)]:
            if inter != "":
                qmasm.abend('Unexpected "%s" in pin right-hand side "%s"' % (inter, rhs))
        return [qmasm.str2bool[t.upper()] for t in self.bool_re.findall(rhs)]

def process_pin(filename, lineno, pin_str):
    "Parse a pin statement into one or more Pin objects and add these to the program."
    lhs_rhs = pin_str.split(":=")
    if len(lhs_rhs) != 2:
        qmasm.abend('Failed to parse pin statement "%s"' % pin_str)
    pin_parser = PinParser()
    lhs_list = pin_parser.parse_lhs(lhs_rhs[0])
    rhs_list = pin_parser.parse_rhs(lhs_rhs[1])
    if len(lhs_list) != len(rhs_list):
        qmasm.abend('Different number of left- and right-hand-side values in "%s" (%d vs. %d)' % (pin_str, len(lhs_list), len(rhs_list)))
    return [Pin(filename, lineno, l, r) for l, r in zip(lhs_list, rhs_list)]

def process_chain(filename, lineno, chain_str):
    "Parse a chain statement into one or more Chain objects and add these to the program."
    # We use the LHS parser from PinParser to parse both sides of the chain.
    lhs_rhs = chain_str.split("=")
    if len(lhs_rhs) != 2:
        qmasm.abend('Failed to parse chain statement "%s"' % chain_str)
    pin_parser = PinParser()
    lhs_list = pin_parser.parse_lhs(lhs_rhs[0])
    rhs_list = pin_parser.parse_lhs(lhs_rhs[1])  # Note use of parse_lhs to parse the RHS.
    if len(lhs_list) != len(rhs_list):
        qmasm.abend('Different number of left- and right-hand-side values in "%s" (%d vs. %d)' % (chain_str, len(lhs_list), len(rhs_list)))
    return [Chain(filename, lineno, l, r) for l, r in zip(lhs_list, rhs_list)]

def process_alias(filename, lineno, alias_str):
    "Parse an alias statement into one or more Alias objects and add these to the program."
    # We use the LHS parser from PinParser to parse both sides of the alias.
    lhs_rhs = alias_str.split("<->")
    if len(lhs_rhs) != 2:
        qmasm.abend('Failed to parse alias statement "%s"' % alias_str)
    pin_parser = PinParser()
    lhs_list = pin_parser.parse_lhs(lhs_rhs[0])
    rhs_list = pin_parser.parse_lhs(lhs_rhs[1])  # Note use of parse_lhs to parse the RHS.
    if len(lhs_list) != len(rhs_list):
        qmasm.abend('Different number of left- and right-hand-side values in "%s" (%d vs. %d)' % (alias_str, len(lhs_list), len(rhs_list)))
    return [Alias(filename, lineno, l, r) for l, r in zip(lhs_list, rhs_list)]
