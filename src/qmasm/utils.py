###################################
# QMASM utility functions         #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from collections import defaultdict
import copy
import math
import qmasm
import sys

class RemainingNextException(Exception):
    'This exception is thrown if a "!next." directive can\'t be replaced.'
    pass

class Utilities(object):
    "Provide various utility functions as mixins for QMASM."

    def abend(self, str):
        "Abort the program on an error."
        sys.stderr.write("%s: %s\n" % (sys.argv[0], str))
        sys.exit(1)

    def warn(self, str):
        "Issue a warning message but continue execution."
        sys.stderr.write("%s: Warning: %s\n" % (sys.argv[0], str))

    def symbol_to_number(self, sym, prefix=None, next_prefix=None):
        "Map from a symbol to a number, creating a new association if necessary."
        # Replace "!next." by substituting prefixes in the name.
        if "!next." in sym:
            if prefix == None or next_prefix == None:
                raise RemainingNextException
            sym = sym.replace(prefix + "!next.", next_prefix)

        # Return the symbol's logical qubit number or allocate a new one.
        try:
            return self.sym_map.to_number(sym)
        except KeyError:
            return self.sym_map.new_symbol(sym)

    def apply_prefix(self, sym, prefix=None, next_prefix=None):
        "Apply a prefix to a symbol name, replacing !next. with the next prefix."
        if prefix != None:
            sym = prefix + sym
        if "!next." in sym:
            if prefix == None or next_prefix == None:
                raise RemainingNextException
            sym = sym.replace(prefix + "!next.", next_prefix)
        return sym

    def canonicalize_strengths(self, strs):
        "Combine edges (A, B) and (B, A) into (A, B) with A < B."
        new_strs = defaultdict(lambda: 0.0)
        for (q1, q2), wt in strs.items():
            if q1 == q2:
                continue          # Discard vertex weights.
            if wt == 0.0:
                continue          # Discard zero weights.
            if q1 > q2:
                q1, q2 = q2, q1   # Canonicalize vertex order.
            new_strs[(q1, q2)] += wt
        return new_strs

class SymbolMapping:
    "Map between symbols and numbers."

    def __init__(self):
        self.sym2num = {}
        self.num2syms = {}
        self.next_sym_num = 0

    def new_symbol(self, sym):
        "Assign the next available number to a symbol."
        if sym in self.sym2num:
            raise Exception("Internal error: Symbol %s is already defined" % sym)
        self.sym2num[sym] = self.next_sym_num
        self.num2syms[self.next_sym_num] = set([sym])
        self.next_sym_num += 1
        return self.sym2num[sym]

    def to_number(self, sym):
        "Map a symbol to a single number."
        return self.sym2num[sym]

    def to_symbols(self, num):
        "Map a number to a set of one or more symbols."
        return self.num2syms[num]

    def max_number(self):
        "Return the maximum symbol number used so far."
        return self.next_sym_num - 1

    def all_numbers(self):
        "Return an unordered list of all numbers used."
        return self.num2syms.keys()

    def all_symbols(self):
        "Return an unordered list of all symbols used."
        return self.sym2num.keys()

    def symbol_number_items(self):
        "Return a list of {symbol, number} pairs."
        return self.sym2num.items()

    def alias(self, sym1, sym2):
        "Make two symbols point to the same number."
        # Ensure that both symbols are defined.
        n_new_defs = 0   # Number of new definitions made
        try:
            num1 = self.sym2num[sym1]
        except KeyError:
            num1 = self.new_symbol(sym1)
            n_new_defs += 1
        try:
            num2 = self.sym2num[sym2]
        except KeyError:
            num2 = self.new_symbol(sym2)
            n_new_defs += 1

        # Do nothing if the symbols are already aliased.
        if num1 == num2:
            return num1

        # Abort if both symbols were previously defined (and do not alias each
        # other).  In this case, we'd need to merge the associated weights and
        # strengths, and we don't currently do that.
        if n_new_defs == 0:
            abend("Unable to alias pre-existing variables %s and %s" % (sym1, sym2))

        # Replace all occurrences of the larger number with the smaller in
        # num2syms.
        new_num = min(num1, num2)
        old_num = max(num1, num2)
        self.num2syms[new_num].update(self.num2syms[old_num])
        del self.num2syms[old_num]

        # Regenerate sym2num from num2syms.
        self.sym2num = {}
        for n, ss in self.num2syms.items():
            for s in ss:
                self.sym2num[s] = n
        return new_num

    def overwrite_with(self, sym2num):
        "Overwrite the map's contents with a given symbol-to-number map."
        self.sym2num = sym2num
        self.num2syms = {}
        for s, n in self.sym2num.items():
            try:
                self.num2syms[n].add(s)
            except KeyError:
                self.num2syms[n] = set([s])
        if len(self.num2syms.keys()) == 0:
            self.next_sym_num = 0
        else:
            self.next_sym_num = max(self.num2syms.keys()) + 1

    def replace_all(self, before_syms, after_syms):
        "Replace each symbol in before_syms with its peer in after_syms."
        sym2num = copy.deepcopy(self.sym2num)
        before_nums = []
        for s in before_syms:
            try:
                before_nums.append(sym2num[s])
            except KeyError:
                abend('Failed to rename nonexistent symbol "%s"' % s)
            del sym2num[s]
        for s, n in zip(after_syms, before_nums):
            sym2num[s] = n
        self.overwrite_with(sym2num)
