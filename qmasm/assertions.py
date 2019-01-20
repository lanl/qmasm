###################################
# Parse an !assert directive      #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import qmasm
import re
import sys

class AssertAST(object):
    "Represent an abstract syntax tree."

    def __init__(self, type, value, kids=[]):
        self.type = type
        self.value = value
        self.kids = kids
        self._str = None   # Memoized string representation

    def _str_helper(self):
        "Do most of the work for the __str__ method."
        nkids = len(self.kids)
        if nkids == 0:
            return str(self.value)
        if nkids == 1:
            if self.type == "unary" and self.value == "id":
                return str(self.kids[0])
            if self.type == "unary":
                return "%s%s" % (self.value, str(self.kids[0]))
            if self.type == "factor" and self.kids[0].type == "expr":
                return "(%s)" % str(self.kids[0])
            return str(self.kids[0])
        if nkids == 2:
            if self.value in ["*", "/", "%", "&", "<<", ">>", "**"]:
                return "%s%s%s" % (str(self.kids[0]), self.value, str(self.kids[1]))
            elif self.type == "conn":
                return "(%s) %s (%s)" % (str(self.kids[0]), self.value, str(self.kids[1]))
            else:
                return "%s %s %s" % (str(self.kids[0]), self.value, str(self.kids[1]))
        if nkids == 3:
            if self.type == "if_expr":
                return "if %s then %s else %s endif" % (str(self.kids[0]), str(self.kids[1]), str(self.kids[2]))
        raise Exception("Internal error parsing (%s, %s)" % (repr(self.type), repr(self.value)))

    def __str__(self):
        if self._str == None:
            self._str = self._str_helper()
        return self._str

    def apply_prefix(self, prefix, next_prefix):
        "Prefix every identifier with a given string."
        if self.type == "ident":
            self.value = qmasm.apply_prefix(self.value, prefix, next_prefix)
        else:
            for k in self.kids:
                k.apply_prefix(prefix, next_prefix)

    class EvaluationError(Exception):
        "Represent an exception thrown during AST evaluation."
        pass

    def _evaluate_ident(self, i2b):
        "Evaluate a variable."
        try:
            bit = i2b[self.value]
            if bit == None:
                raise self.EvaluationError("Unused variable %s" % self.value)
            return bit
        except KeyError:
            raise self.EvaluationError("Undefined variable %s" % self.value)

    def _evaluate_unary(self, i2b, kvals):
        "Evaluate a unary expression."
        if self.value == "-":
            return -kvals[0]
        elif self.value == "~":
            return ~kvals[0]
        elif self.value == "!":
            if kvals[0] == 0:
                return 1
            else:
                return 0
        elif self.value in ["+", "id"]:
            return kvals[0]
        else:
            raise self.EvaluationError('Internal error evaluating unary "%s"' % self.value)

    def _evaluate_arith(self, i2b, kvals):
        "Evaluate an arithmetic expression."
        if self.value == "+":
            return kvals[0] + kvals[1]
        elif self.value == "-":
            return kvals[0] - kvals[1]
        elif self.value == "*":
            return kvals[0] * kvals[1]
        elif self.value == "/":
            return kvals[0] // kvals[1]
        elif self.value == "%":
            return kvals[0] % kvals[1]
        elif self.value == "&":
            return kvals[0] & kvals[1]
        elif self.value == "|":
            return kvals[0] | kvals[1]
        elif self.value == "^":
            return kvals[0] ^ kvals[1]
        elif self.value == "<<":
            return kvals[0] << kvals[1]
        elif self.value == ">>":
            return kvals[0] >> kvals[1]
        elif self.value == "**":
            if kvals[1] < 0:
                raise self.EvaluationError("Negative powers (%d) are not allowed" % kvals[1])
            return kvals[0] ** kvals[1]
        else:
            raise self.EvaluationError("Internal error evaluating arithmetic operator %s" % self.value)

    def _evaluate_rel(self, i2b, kvals):
        "Evaluate a relational expression."
        if self.value == "=":
            return kvals[0] == kvals[1]
        elif self.value == "<>":
            return kvals[0] != kvals[1]
        elif self.value == "<":
            return kvals[0] < kvals[1]
        elif self.value == "<=":
            return kvals[0] <= kvals[1]
        elif self.value == ">":
            return kvals[0] > kvals[1]
        elif self.value == ">=":
            return kvals[0] >= kvals[1]
        else:
            raise self.EvaluationError("Internal error evaluating relational operator %s" % self.value)

    def _evaluate_conn(self, i2b, kvals):
        "Evaluate a logical connective."
        if self.value == "&&":
            return kvals[0] and kvals[1]
        elif self.value == "||":
            return kvals[0] or kvals[1]
        else:
            raise self.EvaluationError("Internal error evaluating logical connective %s" % self.value)

    def _evaluate_if_expr(self, i2b, kvals):
        "Evaluate an if...then...else expression."
        if kvals[0]:
            return kvals[1]
        else:
            return kvals[2]

    def _evaluate_node(self, i2b):
        """Evaluate the AST to either True or False given a mapping from
        identifiers to bits."""
        kvals = [k._evaluate_node(i2b) for k in self.kids]
        if self.type == "ident":
            # Variable
            return self._evaluate_ident(i2b)
        elif self.type == "int":
            # Constant
            return self.value
        elif self.type == "unary":
            # Unary expression
            return self._evaluate_unary(i2b, kvals)
        elif len(kvals) == 1:
            # All other single-child nodes return their child unmodified.
            return kvals[0]
        elif self.type in ["power", "term", "expr"]:
            return self._evaluate_arith(i2b, kvals)
        elif self.type == "rel":
            return self._evaluate_rel(i2b, kvals)
        elif self.type == "conn":
            return self._evaluate_conn(i2b, kvals)
        elif self.type == "if_expr":
            return self._evaluate_if_expr(i2b, kvals)
        else:
            raise self.EvaluationError("Internal error evaluating AST node of type %s, value %s" % (repr(self.type), repr(self.value)))

    def evaluate(self, i2b):
        """Evaluate the AST to either True or False given a mapping from
        identifiers to bits."""
        try:
            return self._evaluate_node(i2b)
        except self.EvaluationError as e:
            qmasm.abend("%s in assertion %s" % (e, self))

class AssertParser(object):
    int_re = re.compile(r'\d+')
    conn_re = re.compile(r'\|\||&&')
    rel_re = re.compile(r'<[=>]?|>=?|=')
    arith_re = re.compile(r'[-+/%&\|^~!]|>>|<<|\*\*?')
    ident_re = re.compile(r'[^-+*/%&\|^~!()<=>\s]+')
    keyword_re = re.compile(r'\b(if|then|else|endif)\b')

    class ParseError(Exception):
        pass

    def lex(self, s):
        "Split a string into tokens (tuples of type and value)."
        tokens = []
        s = s.lstrip()
        while len(s) > 0:
            # Match parentheses.
            if s[0] == "(":
                tokens.append(("lparen", "("))
                s = s[1:].lstrip()
                continue
            if s[0] == ")":
                tokens.append(("rparen", ")"))
                s = s[1:].lstrip()
                continue

            # Match keywords.
            mo = self.keyword_re.match(s)
            if mo != None:
                match = mo.group(0)
                tokens.append((match, match))
                s = s[len(match):].lstrip()
                continue

            # Match positive integers.
            mo = self.int_re.match(s)
            if mo != None:
                match = mo.group(0)
                tokens.append(("int", int(match)))
                s = s[len(match):].lstrip()
                continue

            # Match connectives.
            mo = self.conn_re.match(s)
            if mo != None:
                match = mo.group(0)
                tokens.append(("conn", match))
                s = s[len(match):].lstrip()
                continue

            # Match "<<" and ">>" before we match "<" and ">".
            if len(s) >= 2 and (s[:2] == "<<" or s[:2] == ">>"):
                tokens.append(("arith", s[:2]))
                s = s[2:].lstrip()
                continue

            # Match relational operators.
            mo = self.rel_re.match(s)
            if mo != None:
                match = mo.group(0)
                tokens.append(("rel", match))
                s = s[len(match):].lstrip()
                continue

            # Match "**" before we match "*".
            if len(s) >= 2 and s[:2] == "**":
                tokens.append(("power", s[:2]))
                s = s[2:].lstrip()
                continue

            # Match arithmetic operators.
            mo = self.arith_re.match(s)
            if mo != None:
                match = mo.group(0)
                tokens.append(("arith", match))
                s = s[len(match):].lstrip()
                continue

            # Everything else is an identifier.
            mo = self.ident_re.match(s)
            if mo != None:
                match = mo.group(0)
                tokens.append(("ident", match))
                s = s[len(match):].lstrip()
                continue
            raise self.ParseError("Failed to parse %s" % s)
        tokens.append(("EOF", "EOF"))
        return tokens

    def advance(self):
        "Advance to the next symbol."
        self.tokidx += 1
        self.sym = self.tokens[self.tokidx]

    def accept(self, ty):
        """Advance to the next token if the current token matches a given
        token type and return True.  Otherwise, return False."""
        if self.sym[0] == ty:
            self.advance()
            return True
        return False

    def expect(self, ty):
        """Advance to the next token if the current token matches a given
        token.  Otherwise, fail."""
        if not self.accept(ty):
            raise self.ParseError("Expected %s but saw %s" % (ty, repr(self.sym)))

    def if_expr(self):
        "Return an if...then...else expression."
        self.expect("if")
        cond = self.conjunction()
        self.expect("then")
        then_expr = self.expression()
        self.expect("else")
        else_expr = self.expression()
        self.expect("endif")
        return AssertAST("if_expr", None, [cond, then_expr, else_expr])

    def factor(self):
        "Return a factor (variable, integer, or expression)."
        val = self.sym[1]
        if self.accept("ident"):
            child = AssertAST("ident", val)
        elif self.accept("int"):
            child = AssertAST("int", val)
        elif self.accept("lparen"):
            child = self.expression()
            self.expect("rparen")
        elif self.sym[0] == "arith":
            child = self.unary()
        elif self.sym[0] == "if":
            child = self.if_expr()
        elif val == "EOF":
            raise self.ParseError("Parse error at end of expression")
        else:
            raise self.ParseError('Parse error at "%s"' % val)
        return AssertAST("factor", None, [child])

    def power(self):
        "Return a factor or a factor raised to the power of a second factor."
        f1 = self.factor()
        op = self.sym[1]
        if self.sym[0] == "power" and op == "**":
            self.advance()
            f2 = self.power()
            return AssertAST("power", op, [f1, f2])
        return AssertAST("power", op, [f1])

    def unary(self):
        "Return a unary operator applied to a power."
        op = self.sym[1]
        if op in ["+", "-", "~", "!"]:
            self.advance()
        else:
            op = "id"
        return AssertAST("unary", op, [self.power()])

    def term(self):
        "Return a term (product of an optional term and a unary)."
        # Produce a list of ASTs representing unaries.
        u = self.unary()
        ops = [self.sym[1]]
        asts = [u]
        while self.sym[0] == "arith" and ops[-1] in ["*", "/", "%", "&", "<<", ">>"]:
            self.advance()
            u = self.unary()
            ops.append(self.sym[1])
            asts.append(u)

        # Handle the trivial case of the identity operation.
        if len(asts) == 1:
            return AssertAST("term", ops[0], asts)

        # Merge the ASTs in a left-associative fashion into a single AST.
        ops.pop()
        while len(asts) > 1:
            asts = [AssertAST("term", ops[0], [asts[0], asts[1]])] + asts[2:]
            ops.pop(0)
        return asts[0]

    def expression(self):
        "Return an expression (sum of an optional expression and a term)."
        # Produce a list of ASTs representing terms.
        t = self.term()
        ops = [self.sym[1]]
        asts = [t]
        while self.sym[0] == "arith" and ops[-1] in ["+", "-", "|", "^"]:
            self.advance()
            t = self.term()
            ops.append(self.sym[1])
            asts.append(t)

        # Handle the trivial case of the identity operation.
        if len(asts) == 1:
            return AssertAST("expr", ops[0], asts)

        # Merge the ASTs in a left-associative fashion into a single AST.
        ops.pop()
        while len(asts) > 1:
            asts = [AssertAST("expr", ops[0], [asts[0], asts[1]])] + asts[2:]
            ops.pop(0)
        return asts[0]

    def comparison(self):
        "Return a comparison of exactly two expressions."
        e1 = self.expression()
        op = self.sym[1]
        if self.sym[0] != "rel":
            raise self.ParseError('Expected a relational operator but saw "%s"' % op)
        self.advance()
        e2 = self.expression()
        return AssertAST("rel", op, [e1, e2])

    def disjunction(self):
        "Return an disjunction of one or two comparisons."
        c1 = self.comparison()
        op = self.sym[1]
        if op != "&&":
            return AssertAST("conn", None, [c1])
        self.advance()
        c2 = self.disjunction()
        return AssertAST("conn", op, [c1, c2])

    def conjunction(self):
        "Return a conjunction of one or two disjunctions."
        d1 = self.disjunction()
        op = self.sym[1]
        if op != "||":
            return AssertAST("conn", None, [d1])
        self.advance()
        d2 = self.conjunction()
        return AssertAST("conn", op, [d1, d2])

    def parse(self, s):
        "Parse a relational expression into an AST"
        self.tokens = self.lex(s)
        self.tokidx = -1
        self.advance()
        try:
            ast = self.conjunction()
            if self.sym[0] != "EOF":
                raise self.ParseError('Parse error at "%s"' % self.sym[1])
        except self.ParseError as e:
            qmasm.abend('%s in "%s"' % (e, s))
        return ast
