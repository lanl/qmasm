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
            if type == "unary":
                return "%s%s" % (self.value, str(self.kids[0]))
            if type == "factor" and self.kids[0].type == "expr":
                return "(%s)" % str(self.kids[0])
            return str(self.kids[0])
        if nkids == 2:
            if self.value in ["*", "/", "%", "&"]:
                return "%s%s%s" % (str(self.kids[0]), self.value, str(self.kids[1]))
            elif self.type == "conn":
                return "(%s) %s (%s)" % (str(self.kids[0]), self.value, str(self.kids[1]))
            else:
                return "%s %s %s" % (str(self.kids[0]), self.value, str(self.kids[1]))
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
        elif self.value == "+":
            return kvals[0]
        else:
            raise self.EvaluationError("Internal error evaluating unary %s" % self.value)

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
            return kvals[0] ^ kvals[1]
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
        elif self.type in ["term", "expr"]:
            return self._evaluate_arith(i2b, kvals)
        elif self.type == "rel":
            return self._evaluate_rel(i2b, kvals)
        elif self.type == "conn":
            return self._evaluate_conn(i2b, kvals)
        else:
            raise self.EvaluationError("Internal error evaluating AST node type %s" % self.value)

    def evaluate(self, i2b):
        """Evaluate the AST to either True or False given a mapping from
        identifiers to bits."""
        try:
            return self._evaluate_node(i2b)
        except self.EvaluationError, e:
            qmasm.abend("%s in assertion %s" % (e, self))

class AssertParser(object):
    int_re = re.compile(r'\d+')
    conn_re = re.compile(r'\|\||&&')
    rel_re = re.compile(r'<[=>]?|>=?|=')
    arith_re = re.compile(r'[-+*/%&\|^~]')
    ident_re = re.compile(r'[^-+*/%&\|^~()<=>\s]+')

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

            # Match relational operators.
            mo = self.rel_re.match(s)
            if mo != None:
                match = mo.group(0)
                tokens.append(("rel", match))
                s = s[len(match):].lstrip()
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
        else:
            raise self.ParseError("Parse error at %s" % val)
        return AssertAST("factor", None, [child])

    def unary(self):
        "Return a unary operator applied to a factor."
        op = self.sym[1]
        if op not in ["+", "-", "~"]:
            raise self.ParseError('unexpected unary operator "%s"' % op)
        self.advance()
        return AssertAST("unary", op, self.factor())

    def term(self):
        "Return a term (product of one or two factors)."
        f1 = self.factor()
        op = self.sym[1]
        if self.sym[0] == "arith" and op in ["*", "/", "%", "&"]:
            self.advance()
            f2 = self.term()
            return AssertAST("term", op, [f1, f2])
        return AssertAST("term", op, [f1])

    def expression(self):
        "Return an expression (sum of one or two terms)."
        t1 = self.term()
        op = self.sym[1]
        if self.sym[0] != "arith" or op not in ["+", "-", "|", "^"]:
            return AssertAST("expr", None, [t1])
        self.advance()
        t2 = self.expression()
        return AssertAST("expr", op, [t1, t2])

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
        "Parse an relational expression into an AST"
        self.tokens = self.lex(s)
        self.tokidx = -1
        self.advance()
        ast = self.conjunction()
        if self.sym[0] != "EOF":
            raise self.ParseError('Parse error at "%s"' % self.sym[1])
        return ast
