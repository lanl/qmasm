###################################
# Parse a QMASM source file       #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import copy
import os
import qmasm
import re
import string
import sys

# Define a function that aborts the program, reporting an invalid
# input line as part of the error message.
def error_in_line(filename, lineno, str):
    sys.stderr.write('%s:%d: error: %s\n' % (filename, lineno, str))
    sys.exit(1)

class Environment(object):
    "Maintain a variable environment as a stack of scopes."

    toks_re = re.compile(r'([^-+*/%&\|^~!()<=>#\s\[:\]]+)')  # Regex to split a symbol into tokens (cf. AssertParser.ident_re with square brackets and colons added)

    def __init__(self):
        self.stack = [{}]
        self.self_copy = None   # Regenerate a copy of self if equal to None

    def __getitem__(self, key):
        "Search each scope in turn for the given key."
        for i in range(len(self.stack) - 1, -1, -1):
            try:
                return self.stack[i][key]
            except KeyError:
                pass
        raise KeyError(key)

    def __setitem__(self, key, val):
        self.stack[-1][key] = val
        self.self_copy = None

    def push(self):
        "Push a new scope on the environment stack."
        self.stack.append({})
        self.self_copy = None

    def pop(self):
        "Discard the top of the environment stack."
        self.stack.pop()
        self.self_copy = None

    def keys(self):
        "Return all keys in all scopes."
        d = set()
        for s in self.stack:
            d.update(s.keys())
        return list(d)

    def __copy__(self):
        "Return an independent copy of the environment stack."
        if self.self_copy == None:
            self.self_copy = copy.deepcopy(self)
        return self.self_copy

    def sub_syms(self, sym):
        """Substitute values for variables encountered in a given symbol name
        or list of symbol names."""
        if isinstance(sym, (list,)):
            # Recursively process each list element in turn.
            return [self.sub_syms(s) for s in sym]
        elif isinstance(sym, (str,)):
            # Process a single word.
            toks = self.toks_re.split(sym)
            for i in range(len(toks)):
                try:
                    toks[i] = str(self[toks[i]])
                except KeyError:
                    pass
            return "".join(toks)
        else:
            # Non-strings are returned unmodified.
            return sym

# I'm too lazy to write another parser so I'll simply define an
# alternative entry point to the assertion parser.
class ExprParser(qmasm.AssertParser):
    "Parse an arithmetic expression."

    def parse(self, filename, lineno, s):
        "Parse an arithmetic expression into an AST"
        self.tokens = self.lex(s)
        self.tokidx = -1
        self.advance()
        try:
            ast = self.expression()
            if self.sym[0] != "EOF":
                raise self.ParseError('Parse error at "%s"' % self.sym[1])
        except self.ParseError as e:
            sys.stderr.write('%s:%d: error: %s in "%s"\n' % (filename, lineno, e, s))
            sys.exit(1)
        return ast

# I'm too lazy to write another parser so I'll simply define an
# alternative entry point to the assertion parser.
class RelationParser(qmasm.AssertParser):
    "Parse a relational expression."

    def parse(self, filename, lineno, s):
        "Parse a relational expression into an AST"
        self.tokens = self.lex(s)
        self.tokidx = -1
        self.advance()
        try:
            ast = self.conjunction()
            if self.sym[0] != "EOF":
                raise self.ParseError('Parse error at "%s"' % self.sym[1])
        except self.ParseError as e:
            sys.stderr.write('%s:%d: error: %s in "%s"\n' % (filename, lineno, e, s))
            sys.exit(1)
        return ast

class LoopIterator(object):
    '''Iterate over arithmetic and geometric integer sequences.  It is
    assumed that the penultimate element of the input list is the string
    "...".'''

    def __init__(self, filename, lineno, rhs):
        # Prepare the iterations we intend to perform.
        self.first_val = rhs[0]
        last_val = rhs[-1]
        if len(rhs) == 3:
            # "<x_0> ... <x_n>" indicates an arithmetic progression with a
            # delta of +/- 1.
            if self.first_val <= last_val:
                self.increment = lambda x: x + 1
            else:
                self.increment = lambda x: x - 1
        else:
            # For "<x_0> <x_1> <x_2> ... <x_n>", compute the progression
            # type and increment.
            adeltas = set([rhs[i + 1] - rhs[i] for i in range(len(rhs) - 3)])
            if len(adeltas) == 1:
                # Arithmetic progression
                delta = list(adeltas)[0]
                self.increment = lambda x: x + delta
            else:
                # Geometric progression
                try:
                    mdeltas = set([rhs[i + 1] // rhs[i] for i in range(len(rhs) - 3)])
                except ZeroDivisionError:
                    mdeltas = set(["multiple", "values"])  # Force the next test to fail.
                if len(mdeltas) != 1:
                    error_in_line(filename, lineno,
                                  'Failed to interpret "%s, %d" as either an arithmetic or geometric progression' %
                                  (", ".join([str(r) for r in rhs[:-2]]), last_val))
                delta = list(mdeltas)[0]
                if delta == 0:
                    error_in_line(filename, lineno, "Decreasing geometric progressions are not currently supported")
                self.increment = lambda x: x * delta

        # Determine when to stop.
        second_val = self.increment(self.first_val)
        if self.first_val < second_val:
            self.finished = lambda x: x > last_val
        else:
            self.finished = lambda x: x < last_val

    def __iter__(self):
        self.next_val = self.first_val
        return self

    def __next__(self):
        if self.finished(self.next_val):
            raise StopIteration
        head = self.next_val
        self.next_val = self.increment(head)
        return head

    def next(self):
        # Python 2 wrapper for __next__
        return self.__next__()

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

    def update_qmi(self, prefix, next_prefix, problem):
        num = qmasm.symbol_to_number(prefix + self.sym, prefix, next_prefix)
        problem.weights[num] += self.weight

class Chain(Statement):
    "Chain between qubits."
    def __init__(self, filename, lineno, sym1, sym2):
        super(Chain, self).__init__(filename, lineno)
        self.sym1 = sym1
        self.sym2 = sym2

    def as_str(self, prefix=""):
        return "%s%s = %s%s" % (prefix, self.sym1, prefix, self.sym2)

    def update_qmi(self, prefix, next_prefix, problem):
        num1 = qmasm.symbol_to_number(prefix + self.sym1, prefix, next_prefix)
        num2 = qmasm.symbol_to_number(prefix + self.sym2, prefix, next_prefix)
        if num1 == num2:
            self.error_in_line("A chain cannot connect a spin to itself")
        elif num1 > num2:
            num1, num2 = num2, num1
        problem.chains.add((num1, num2))

class AntiChain(Statement):
    "AntiChain between qubits."
    def __init__(self, filename, lineno, sym1, sym2):
        super(AntiChain, self).__init__(filename, lineno)
        self.sym1 = sym1
        self.sym2 = sym2

    def as_str(self, prefix=""):
        return "%s%s /= %s%s" % (prefix, self.sym1, prefix, self.sym2)

    def update_qmi(self, prefix, next_prefix, problem):
        num1 = qmasm.symbol_to_number(prefix + self.sym1, prefix, next_prefix)
        num2 = qmasm.symbol_to_number(prefix + self.sym2, prefix, next_prefix)
        if num1 == num2:
            self.error_in_line("An anti-chain cannot connect a spin to itself")
        elif num1 > num2:
            num1, num2 = num2, num1
        problem.antichains.add((num1, num2))

class Pin(Statement):
    "Pinning of a qubit to true or false."
    def __init__(self, filename, lineno, sym, goal):
        super(Pin, self).__init__(filename, lineno)
        self.sym = sym
        self.goal = goal

    def as_str(self, prefix=""):
        return "%s%s := %s" % (prefix, self.sym, self.goal)

    def update_qmi(self, prefix, next_prefix, problem):
        num = qmasm.symbol_to_number(prefix + self.sym, prefix, next_prefix)
        problem.pinned.append((num, self.goal))

class Alias(Statement):
    "Alias of one symbol to another."
    def __init__(self, filename, lineno, sym1, sym2):
        super(Alias, self).__init__(filename, lineno)
        self.sym1 = sym1
        self.sym2 = sym2

    def as_str(self, prefix=""):
        return "%s%s <-> %s%s" % (prefix, self.sym1, prefix, self.sym2)

    def update_qmi(self, prefix, next_prefix, problem):
        sym1 = prefix + self.sym1
        sym2 = prefix + self.sym2
        if next_prefix != None:
            sym1 = sym1.replace(prefix + "!next.", next_prefix)
            sym2 = sym2.replace(prefix + "!next.", next_prefix)
        qmasm.sym_map.alias(sym1, sym2)

class Strength(Statement):
    "Coupler strength between two qubits."
    def __init__(self, filename, lineno, sym1, sym2, strength):
        super(Strength, self).__init__(filename, lineno)
        self.sym1 = sym1
        self.sym2 = sym2
        self.strength = strength

    def as_str(self, prefix=""):
        return "%s%s %s%s %s" % (prefix, self.sym1, prefix, self.sym2, self.strength)

    def update_qmi(self, prefix, next_prefix, problem):
        num1 = qmasm.symbol_to_number(prefix + self.sym1, prefix, next_prefix)
        num2 = qmasm.symbol_to_number(prefix + self.sym2, prefix, next_prefix)
        if num1 == num2:
            self.error_in_line("A coupler cannot connect a spin to itself")
        elif num1 > num2:
            num1, num2 = num2, num1
        problem.strengths[(num1, num2)] += self.strength

class Assert(Statement):
    "Instantiation of a run-time assertion."
    parser = qmasm.AssertParser()

    def __init__(self, filename, lineno, expr):
        super(Assert, self).__init__(filename, lineno)
        self.expr = expr
        self.ast = self.parser.parse(expr)
        self.ast.compile()

    def as_str(self, prefix=""):
        if prefix == "":
            ast = self.ast
        else:
            ast = copy.deepcopy(self.ast)
            ast.apply_prefix(prefix, None)
            ast.compile()
        return "!assert " + str(ast)

    def update_qmi(self, prefix, next_prefix, problem):
        if prefix == "":
            ast = self.ast
        else:
            ast = copy.deepcopy(self.ast)
            ast.apply_prefix(prefix, next_prefix)
            ast.compile()
        problem.assertions.append(ast)

class MacroUse(Statement):
    "Instantiation of a macro definition."
    def __init__(self, filename, lineno, name, body, prefixes):
        super(MacroUse, self).__init__(filename, lineno)
        self.name = name
        self.body = body
        self.prefixes = prefixes

    def as_str(self, prefix=""):
        stmt_strs = []
        nprefixes = len(self.prefixes)
        if nprefixes == 0:
            # No prefixes -- display the macro body in the current scope (i.e.,
            # using the given prefix).
            for stmt in self.body:
                sstr = stmt.as_str(prefix)
                if "!next." not in sstr:
                    stmt_strs.append(sstr)
        else:
            # At least one prefix -- disaply the macro body in a new scope
            # (i.e., by augmenting the given prefix with each new prefix in
            # turn).
            for p in range(nprefixes):
                pfx = self.prefixes[p]
                for stmt in self.body:
                    sstr = stmt.as_str(prefix + pfx)
                    if "!next." in sstr:
                        if p == nprefixes - 1:
                            # Drop statements that use "!next." if there's
                            # no next prefix.
                            continue
                        next_pfx = self.prefixes[p + 1]
                        sstr = sstr.replace(prefix + pfx + "!next.", prefix + next_pfx)
                    stmt_strs.append(sstr)
        return "\n".join(stmt_strs)

    def update_qmi(self, prefix, next_prefix, problem):
        nprefixes = len(self.prefixes)
        if nprefixes == 0:
            # No prefixes -- import the macro body into the current scope
            # (i.e., using the given prefix).
            for stmt in self.body:
                try:
                    stmt.update_qmi(prefix, None, problem)
                except qmasm.utils.RemainingNextException:
                    pass
        else:
            # At least one prefix -- import the macro body into a new scope
            # (i.e., by augmenting the given prefix with each new prefix in
            # turn).
            for p in range(nprefixes):
                pfx = prefix + self.prefixes[p]
                if p == nprefixes - 1:
                    next_pfx = None
                else:
                    next_pfx = prefix + self.prefixes[p + 1]
                for stmt in self.body:
                    try:
                        stmt.update_qmi(pfx, next_pfx, problem)
                    except qmasm.utils.RemainingNextException:
                        pass

class FileParser(object):
    "Parse a QMASM file."

    def __init__(self):
        self.macros = {}        # Map from a macro name to a list of Statement objects
        self.current_macro = (None, [])   # Macro currently being defined (name and statements)
        self.target = qmasm.program   # Reference to either the program or the current macro
        self.env = Environment()      # Stack of maps from compile-time variable names to values
        self.expr_parser = ExprParser()     # Expression parser
        self.rel_parser = RelationParser()  # Relation parser

        # Establish a mapping from a first-field directive to a parsing function.
        self.dir_to_func = {
            "!include":     self.parse_line_include,
            "!assert":      self.parse_line_assert,
            "!let":         self.parse_line_let,
            "!begin_macro": self.parse_line_begin_macro,
            "!end_macro":   self.parse_line_end_macro,
            "!use_macro":   self.parse_line_use_macro,
            "!alias":       self.parse_line_sym_alias
        }

    def is_float(self, str):
        "Return True if a string can be treated as a float."
        try:
            float(str)
            return True
        except ValueError:
            return False

    def split_line_into_fields(self, line):
        "Discard comments, then split the line into fields."
        try:
            line = line[:line.index("#")]
        except ValueError:
            pass
        return line.split()

    def find_file_in_path(self, pathnames, filename):
        "Search a list of directories for a file."
        for pname in pathnames:
            fname = os.path.join(pname, filename)
            if os.path.exists(fname):
                return fname
            fname_qmasm = fname + ".qmasm"
            if os.path.exists(fname_qmasm):
                return fname_qmasm
        return None

    def parse_line_include(self, filename, lineno, fields):
        "Parse an !include directive."
        # "!include" "<filename>" -- process a named auxiliary file.
        if len(fields) != 2:
            error_in_line(filename, lineno, "Expected a filename to follow !include")
        incname = " ".join(fields[1:])
        if len(incname) >= 2 and incname[0] == "<" and incname[-1] == ">":
            # Search QMASMPATH for the filename.
            incname = incname[1:-1]
            try:
                qmasmpath = os.environ["QMASMPATH"].split(":")
                qmasmpath.append(".")
            except KeyError:
                qmasmpath = ["."]
            found_incname = self.find_file_in_path(qmasmpath, incname)
            if found_incname != None:
                incname = found_incname
        elif len(incname) >= 2 and incname[0] == '"' and incname[-1] == '"':
            # Search only the current directory for the filename.
            incname = incname[1:-1]
            found_incname = self.find_file_in_path(["."], incname)
            if found_incname != None:
                incname = found_incname
        else:
            error_in_line(filename, lineno, 'Failed to parse "%s"' % (" ".join(fields)))
        try:
            incfile = open(incname)
        except IOError:
            error_in_line(filename, lineno, 'Failed to open %s for input' % incname)
        self.process_file(incname, incfile)
        incfile.close()

    def parse_line_assert(self, filename, lineno, fields):
        "Parse an !assert directive."
        # "!assert" <expr> -- assert a property that must be true at run time.
        if len(fields) < 2:
            error_in_line(filename, lineno, "Expected an expression to follow !assert")
        self.target.append(Assert(filename, lineno, " ".join(self.env.sub_syms(fields[1:]))))

    def parse_line_let(self, filename, lineno, fields):
        "Parse a !let directive."
        # "!let" <name> := <expr> -- evaluate <expr> and assign the result to <name>.
        if len(fields) < 4 or fields[2] != ":=":
            error_in_line(filename, lineno, 'Expected a variable name, ":=", and an expression to follow !let')
        lhs = fields[1]
        if len(fields) == 4:
            # Handle the case of "!let" <name> := <symbol>.
            is_sym = self.env.toks_re.match(fields[3])
            if is_sym != None and is_sym.group(0) == fields[3]:
                self.env[lhs] = self.env.sub_syms(fields[3])
                return
        ast = self.expr_parser.parse(filename, lineno, " ".join(self.env.sub_syms(fields[3:])))
        ast.compile()
        rhs = ast.evaluate(dict(self.env))
        self.env[lhs] = rhs

    def parse_line_begin_macro(self, filename, lineno, fields):
        "Parse a !begin_macro directive."
        # "!begin_macro" <name> -- begin a macro definition.
        if len(fields) != 2:
            error_in_line(filename, lineno, "Expected a macro name to follow !begin_macro")
        name = fields[1]
        if name in self.macros:
            error_in_line(self, filename, lineno, "Macro %s is multiply defined" % name)
        if self.current_macro[0] != None:
            error_in_line(filename, lineno, "Nested macros are not supported")
        self.current_macro = (name, [])
        self.target = self.current_macro[1]
        self.env.push()

    def parse_line_end_macro(self, filename, lineno, fields):
        "Parse an !end_macro directive."
        # "!end_macro" <name> -- end a macro definition.
        if len(fields) != 2:
            error_in_line(filename, lineno, "Expected a macro name to follow !end_macro")
        name = fields[1]
        if self.current_macro[0] == None:
            error_in_line(filename, lineno, "Ended macro %s with no corresponding begin" % name)
        if self.current_macro[0] != name:
            error_in_line(filename, lineno, "Ended macro %s after beginning macro %s" % (name, self.current_macro[0]))
        self.macros[name] = self.current_macro[1]
        self.target = qmasm.program
        self.current_macro = (None, [])
        self.env.pop()

    def parse_line_weight(self, filename, lineno, fields):
        "Parse a qubit weight."
        # <symbol> <weight> -- increment a symbol's point weight.
        if len(fields) != 2:
            error_in_line(filename, lineno, "Internal error in parse_line_weight")
        try:
            val = float(self.env.sub_syms(fields[1]))
        except ValueError:
            error_in_line(filename, lineno, 'Failed to parse "%s %s" as a symbol followed by a numerical weight' % (fields[0], fields[1]))
        self.target.append(Weight(filename, lineno, self.env.sub_syms(fields[0]), val))

    def parse_line_chain(self, filename, lineno, fields):
        "Parse a qubit chain."
        # <symbol_1> = <symbol_2> -- create a chain between <symbol_1>
        # and <symbol_2>.
        if len(fields) != 3 or fields[1] != "=":
            error_in_line(filename, lineno, "Internal error in parse_line_chain")
        code = "%s = %s" % (self.env.sub_syms(fields[0]), self.env.sub_syms(fields[2]))
        self.target.extend(self.process_chain(filename, lineno, code))

    def parse_line_antichain(self, filename, lineno, fields):
        "Parse a qubit anti-chain."
        # <symbol_1> /= <symbol_2> -- create an anti-chain between <symbol_1>
        # and <symbol_2>.
        if len(fields) != 3 or fields[1] != "/=":
            error_in_line(filename, lineno, "Internal error in parse_line_antichain")
        code = "%s /= %s" % (self.env.sub_syms(fields[0]), self.env.sub_syms(fields[2]))
        self.target.extend(self.process_antichain(filename, lineno, code))

    def parse_line_pin(self, filename, lineno, fields):
        "Parse a qubit pin."
        # <symbol> := <value> -- force symbol <symbol> to have value <value>.
        if len(fields) != 3 or fields[1] != ":=":
            error_in_line(filename, lineno, "Internal error in parse_line_pin")
        code = "%s := %s" % (self.env.sub_syms(fields[0]), self.env.sub_syms(fields[2]))
        self.target.extend(self.process_pin(filename, lineno, code))

    def parse_line_alias(self, filename, lineno, fields):
        "Parse a qubit alias."
        # <symbol_1> <-> <symbol_2> -- make <symbol_1> an alias of <symbol_2>.
        if len(fields) != 3 or fields[1] != "<->":
            error_in_line(filename, lineno, "Internal error in parse_line_alias")
        code = "%s <-> %s" % (self.env.sub_syms(fields[0]), self.env.sub_syms(fields[2]))
        self.target.extend(self.process_alias(filename, lineno, code))

    def parse_line_strength(self, filename, lineno, fields):
        "Parse a coupler strength."
        # <symbol_1> <symbol_2> <strength> -- increment a coupler strength.
        if len(fields) != 3:
            error_in_line(filename, lineno, "Internal error in parse_line_strength")
        try:
            strength = float(self.env.sub_syms(fields[2]))
        except ValueError:
            error_in_line(filename, lineno, 'Failed to parse "%s" as a number' % fields[2])
        self.target.append(Strength(filename, lineno, self.env.sub_syms(fields[0]), self.env.sub_syms(fields[1]), strength))

    def parse_line_use_macro(self, filename, lineno, fields):
        "Parse a !use_macro directive."
        # "!use_macro" <macro_name> [<instance_name> ...] -- instantiate a
        # macro using <instance_name> as each variable's prefix.
        if len(fields) < 2:
            error_in_line(filename, lineno, "Expected a macro name to follow !use_macro")
        name = self.env.sub_syms(fields[1])
        prefixes = [self.env.sub_syms(p) + "." for p in fields[2:]]
        try:
            self.target.append(MacroUse(filename, lineno, name, self.macros[name], prefixes))
        except KeyError:
            error_in_line(filename, lineno, "Unknown macro %s" % name)

    def parse_line_sym_alias(self, filename, lineno, fields):
        "Parse an !alias directive."
        sys.stderr.write('%s:%d: warning: !alias is deprecated; use "!let %s := %s" instead\n' % (filename, lineno, fields[1], fields[2]))
        if len(fields) != 3:
            error_in_line(filename, lineno, "Expected a symbol name and replacement to follow !alias")
        self.env[fields[1]] = self.env.sub_syms(fields[2])

    def process_if(self, filename, lineno, fields, all_lines):
        """Parse and process an !if directive.  Recursively parse the remaining
        file contents."""
        if len(fields) < 2:
            error_in_line(filename, lineno, "Expected a relational expression to follow !if")

        # Scan ahead for the matching !end_if.
        ends_needed = 1
        else_idx = -1
        for idx in range(1, len(all_lines)):
            # Split the line into fields and apply text aliases.
            lno, line = all_lines[idx]
            if line == "":
                continue
            flds = self.split_line_into_fields(line)
            if flds[0] == "!if":
                ends_needed += 1
            elif flds[0] == "!else":
                if len(flds) != 1:
                    error_in_line(filename, lno, "Unexpected text after !else")
                if else_idx != -1:
                    error_in_line(filename, lno, "An !else matching line %d's !if already appeared in line %d" % (lineno, all_lines[else_idx][0]))
                else_idx = idx
            elif flds[0] == "!end_if":
                if len(flds) != 1:
                    error_in_line(filename, lno, "Unexpected text after !end_if")
                ends_needed -= 1
                if ends_needed == 0:
                    break

        # If the condition is true, parse the body.  Otherwise do nothing.
        if ends_needed != 0:
            error_in_line(filename, lineno, "Failed to find a matching !end_if directive")
        end_idx = idx
        ast = self.rel_parser.parse(filename, lineno, " ".join(self.env.sub_syms(fields[1:])))
        ast.compile()
        rhs = ast.evaluate(dict(self.env))
        if rhs:
            # Evaluate the then clause in a new scope.
            self.env.push()
            if else_idx == -1:
                # No else clause
                self.process_file_contents(filename, all_lines[1:end_idx])
            else:
                # else clause
                self.process_file_contents(filename, all_lines[1:else_idx])
            self.env.pop()
        elif else_idx != -1:
            # Evaluate the else clause in a new scope.
            self.env.push()
            self.process_file_contents(filename, all_lines[else_idx+1:end_idx])
            self.env.pop()

        # Process the rest of the file.
        self.process_file_contents(filename, all_lines[end_idx+1:])

    def process_for(self, filename, lineno, fields, all_lines):
        """Parse and process a !for directive.  Recursively parse the remaining
        file contents."""
        # Parse the !for line.
        if len(fields) < 4 or fields[2] != ":=":
            error_in_line(filename, lineno, 'Expected a variable name, an ":=", and a comma-separated list to follow !for')
        seq = [s.strip() for s in " ".join(self.env.sub_syms(fields[3:])).split(",")]
        for i in range(len(seq)):
            if seq[i] == "..." and i != len(seq) - 2:
                error_in_line(filename, lineno, '"..." can appear only in the penultimate position in a sequence')

        # Construct an iterator based on the given sequence.
        if "..." in seq:
            seq_ints = []
            for i in range(len(seq)):
                if seq[i] == "...":
                    seq_ints.append("...")
                else:
                    ast = self.expr_parser.parse(filename, lineno, seq[i])
                    ast.compile()
                    val = ast.evaluate(dict(self.env))
                    seq_ints.append(val)
            iter = LoopIterator(filename, lineno, seq_ints)
        else:
            iter = seq

        # Scan ahead for the matching !end_for.
        ends_needed = 1
        for idx in range(1, len(all_lines)):
            # Split the line into fields and apply text aliases.
            lno, line = all_lines[idx]
            if line == "":
                continue
            flds = self.split_line_into_fields(line)
            if flds[0] == "!for":
                ends_needed += 1
            elif flds[0] == "!end_for":
                if len(flds) != 1:
                    error_in_line(filename, lno, "Unexpected text after !end_for")
                ends_needed -= 1
                if ends_needed == 0:
                    break

        # Parse the body once per iterator element.
        if ends_needed != 0:
            error_in_line(filename, lineno, "Failed to find a matching !end_for directive")
        end_idx = idx
        for val in iter:
            self.env.push()
            self.env[fields[1]] = val
            self.process_file_contents(filename, all_lines[1:end_idx])
            self.env.pop()

        # Process the rest of the file.
        self.process_file_contents(filename, all_lines[end_idx+1:])


    def process_pin(self, filename, lineno, pin_str):
        "Parse a pin statement into one or more Pin objects and add these to the program."
        lhs_rhs = pin_str.split(":=")
        if len(lhs_rhs) != 2:
            qmasm.abend('Failed to parse pin statement "%s"' % pin_str)
        pin_parser = PinParser()
        lhs_list = pin_parser.parse_lhs(self.env.sub_syms(lhs_rhs[0]))
        rhs_list = pin_parser.parse_rhs(self.env.sub_syms(lhs_rhs[1]))
        if len(lhs_list) != len(rhs_list):
            qmasm.abend('Different number of left- and right-hand-side values in "%s" (%d vs. %d)' % (pin_str, len(lhs_list), len(rhs_list)))
        return [Pin(filename, lineno, l, r) for l, r in zip(lhs_list, rhs_list)]

    def process_chain(self, filename, lineno, chain_str):
        "Parse a chain statement into one or more Chain objects and add these to the program."
        # We use the LHS parser from PinParser to parse both sides of the chain.
        lhs_rhs = chain_str.split("=")
        if len(lhs_rhs) != 2:
            qmasm.abend('Failed to parse chain statement "%s"' % chain_str)
        pin_parser = PinParser()
        lhs_list = pin_parser.parse_lhs(self.env.sub_syms(lhs_rhs[0]))
        rhs_list = pin_parser.parse_lhs(self.env.sub_syms(lhs_rhs[1]))  # Note use of parse_lhs to parse the RHS.
        if len(lhs_list) != len(rhs_list):
            qmasm.abend('Different number of left- and right-hand-side values in "%s" (%d vs. %d)' % (chain_str, len(lhs_list), len(rhs_list)))
        return [Chain(filename, lineno, l, r) for l, r in zip(lhs_list, rhs_list)]

    def process_antichain(self, filename, lineno, antichain_str):
        "Parse an anti-chain statement into one or more AntiChain objects and add these to the program."
        # We use the LHS parser from PinParser to parse both sides of the anti-chain.
        lhs_rhs = antichain_str.split("/=")
        if len(lhs_rhs) != 2:
            qmasm.abend('Failed to parse anti-chain statement "%s"' % antichain_str)
        pin_parser = PinParser()
        lhs_list = pin_parser.parse_lhs(self.env.sub_syms(lhs_rhs[0]))
        rhs_list = pin_parser.parse_lhs(self.env.sub_syms(lhs_rhs[1]))  # Note use of parse_lhs to parse the RHS.
        if len(lhs_list) != len(rhs_list):
            qmasm.abend('Different number of left- and right-hand-side values in "%s" (%d vs. %d)' % (antichain_str, len(lhs_list), len(rhs_list)))
        return [AntiChain(filename, lineno, l, r) for l, r in zip(lhs_list, rhs_list)]

    def process_alias(self, filename, lineno, alias_str):
        "Parse an alias statement into one or more Alias objects and add these to the program."
        # We use the LHS parser from PinParser to parse both sides of the alias.
        lhs_rhs = alias_str.split("<->")
        if len(lhs_rhs) != 2:
            qmasm.abend('Failed to parse alias statement "%s"' % alias_str)
        pin_parser = PinParser()
        lhs_list = pin_parser.parse_lhs(self.env.sub_syms(lhs_rhs[0]))
        rhs_list = pin_parser.parse_lhs(self.env.sub_syms(lhs_rhs[1]))  # Note use of parse_lhs to parse the RHS.
        if len(lhs_list) != len(rhs_list):
            qmasm.abend('Different number of left- and right-hand-side values in "%s" (%d vs. %d)' % (alias_str, len(lhs_list), len(rhs_list)))
        return [Alias(filename, lineno, l, r) for l, r in zip(lhs_list, rhs_list)]

    def process_file_contents(self, filename, all_lines):
        """Parse the contents of a file.  Contents are passed as a list plus an
        initial line number."""
        for idx in range(len(all_lines)):
            # Split the line into fields and apply text aliases.
            lineno, line = all_lines[idx]
            if line == "":
                continue
            fields = self.split_line_into_fields(line)
            nfields = len(fields)
            if nfields == 0:
                continue

            # Process the line.
            if fields[0] == "!if":
                # Special case for !if directives
                return self.process_if(filename, lineno, fields, all_lines[idx:])
            elif fields[0] == "!for":
                # Special case for !for directives
                return self.process_for(filename, lineno, fields, all_lines[idx:])
            try:
                # Parse first-field directives.
                func = self.dir_to_func[fields[0]]
            except KeyError:
                # Reject unrecognized directives.
                if fields[0][0] == "!":
                    if fields[0] in ["!else", "!end_if"]:
                        error_in_line(filename, lineno, "Encountered an %s without a matching !if" % fields[0])
                    elif fields[0] == "!end_for":
                        error_in_line(filename, lineno, "Encountered an !end_for without a matching !for")
                    else:
                        error_in_line(filename, lineno, "Unrecognized directive %s" % fields[0])

                # Prohibit "!next." outside of macros.
                if self.current_macro[0] == None:
                    for f in fields:
                        if "!next." in f:
                            error_in_line(filename, lineno, '"!next." is allowed only within !begin_macro...!end_macro blocks')

                # Parse all lines not containing a directive in the first field.
                if nfields == 2:
                    func = self.parse_line_weight
                elif nfields == 3 and fields[1] == "=":
                    func = self.parse_line_chain
                elif nfields == 3 and fields[1] == "/=":
                    func = self.parse_line_antichain
                elif nfields == 3 and fields[1] == ":=":
                    func = self.parse_line_pin
                elif nfields == 3 and fields[1] == "<->":
                    func = self.parse_line_alias
                elif nfields == 3 and self.is_float(fields[2]):
                    func = self.parse_line_strength
                else:
                    # None of the above
                    error_in_line(filename, lineno, 'Failed to parse "%s"' % line)
            func(filename, lineno, fields)

    def process_file(self, filename, infile):
        """Define a function that parses an input file into an internal
        representation.  This function can be called recursively (due to !include
        directives)."""

        # Read the entire file into a list.
        all_lines = []
        lineno = 1
        for line in infile:
            all_lines.append((lineno, line.strip()))
            lineno += 1
        self.process_file_contents(filename, all_lines)

    def process_files(self, file_list):
        "Parse a list of file(s) into an internal representation."
        if file_list == []:
            # No files were specified: Read from standard input.
            self.process_file("<stdin>", sys.stdin)
            if self.current_macro[0] != None:
                error_in_line(filename, lineno, "Unterminated definition of macro %s" % self.current_macro[0])
        else:
            # Files were specified: Process each in turn.
            for infilename in file_list:
                try:
                    infile = open(infilename)
                except IOError:
                    qmasm.abend('Failed to open %s for input' % infilename)
                self.process_file(infilename, infile)
                if self.current_macro[0] != None:
                    error_in_line(filename, lineno, "Unterminated definition of macro %s" % self.current_macro[0])
                infile.close()

class PinParser(object):
    "Provide methods for parsing a pin statement."

    def __init__(self):
        self.bracket_re = re.compile(r'^\s*(\d+)(\s*(?:\.\.|:)\s*(\d+))?\s*$')
        self.bool_re = re.compile(r'TRUE|FALSE|T|F|0|[-+]?1', re.IGNORECASE)

        # Define synonyms for "true" and "false".
        self.str2bool = {s: True for s in ["1", "+1", "T", "TRUE"]}
        self.str2bool.update({s: False for s in ["0", "-1", "F", "FALSE"]})

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
        return [self.str2bool[t.upper()] for t in self.bool_re.findall(rhs)]
