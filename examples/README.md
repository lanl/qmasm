QMASM examples
=============

This directory contains examples of QMASM code.

One of five
-----------

* Main file: [`1of5.qmasm`](1of5.qmasm)

* Command line: `qmasm --run 1of5.qmasm`

This is a trivial demonstration of QMASM.  The program defines variables *A*, *B*, *C*, *D*, and *E* and outputs all (Boolean) values of those in which exactly one variable is *true*:
```
Solution #1 (energy = -21.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          -1  False
    B          -1  False
    C          -1  False
    D          -1  False
    E          +1  True

Solution #2 (energy = -21.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          -1  False
    B          -1  False
    C          -1  False
    D          +1  True
    E          -1  False

Solution #3 (energy = -21.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          -1  False
    B          -1  False
    C          +1  True
    D          -1  False
    E          -1  False

Solution #4 (energy = -21.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          -1  False
    B          +1  True
    C          -1  False
    D          -1  False
    E          -1  False

Solution #5 (energy = -21.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          +1  True
    B          -1  False
    C          -1  False
    D          -1  False
    E          -1  False
```

Try experimenting with the `--pin` option.  If a variable is pinned to *true*, QMASM will output the single solution that honors that constraint.  If a variable is pinned to *false*, QMASM will output the four solutions.  If *two* variables are pinned to *true*, a situation `1of5.qmasm` does not include in its ground state, QMASM may return no solutions or it may return one or more *incorrect* solutions—whatever exhibits the lowest total energy and doesn't break any chains or pins.

Circuit satisfiability
----------------------

* Main file: [`circsat.qmasm`](circsat.qmasm)

* Helper file: [`gates.qmasm`](gates.qmasm)

* Command line: `qmasm --run --pin="x10 := true" circsat.qmasm`

Given a Boolean expression with *n* inputs and one output, determine the sets of inputs that make the output *true*.  This is a classic NP-complete problem.  `circsat.qmasm` represents a particular 3-input Boolean expression borrowed from the [*Introduction to Algorithms*](https://mitpress.mit.edu/books/introduction-algorithms) textbook's discussion of NP-completeness.

`gates.qmasm` defines macros for various Boolean operators (NOT, 2-input OR, 3-input AND), which are then used by the top-level program, `circsat.qmasm`.  `circsat.qmasm` maps inputs *x1*, *x2*, and *x3* to output *x10*.  By default, the program returns all sets of inputs and the corresponding output.  Pinning *x10* to *true* returns only the solutions to the circuit-satisfiability problem.  In this case, there is only one solution:
```
Solution #1 (energy = -80.00):

    Name(s)  Spin  Boolean
    -------  ----  -------
    x1         +1  True
    x10        +1  True
    x2         +1  True
    x3         -1  False
```

Sorting
-------

* Main file: [`sort4.qmasm`](sort4.qmasm)

* Helper file: [`comparator.qmasm`](comparator.qmasm)

* Command line: `qmasm --run --pin="i1 i2 i3 i4 := 0 1 1 0" sort4.qmasm`

Sort a list of four 1-bit numbers.  `sort4.qmasm` implements a 4-element sorting network from [*Sorting and Searching*](http://www.informit.com/store/art-of-computer-programming-volume-3-sorting-and-searching-9780201896855).  Specify values for inputs *i1*, *i2*, *i3*, and *i4*, and the program will sort these into *o1*, *o2*, *o3*, and *o4*:
```
Solution #1 (energy = -89.50):

    Name(s)  Spin  Boolean
    -------  ----  -------
    i1         -1  False
    i2         +1  True
    i3         +1  True
    i4         -1  False
    o1         -1  False
    o2         -1  False
    o3         +1  True
    o4         +1  True
```
Because there is no clear distinction between inputs and outputs, one can also specify the outputs and receive a list of inputs that would sort to those outputs.  For example, `--pin="o1 o2 o3 o4 := 0 0 1 1"` produces the following 6 solutions:
```
Solution #1 (energy = -89.50):

    Name(s)  Spin  Boolean
    -------  ----  -------
    i1         +1  True
    i2         -1  False
    i3         -1  False
    i4         +1  True
    o1         -1  False
    o2         -1  False
    o3         +1  True
    o4         +1  True

Solution #2 (energy = -89.50):

    Name(s)  Spin  Boolean
    -------  ----  -------
    i1         -1  False
    i2         -1  False
    i3         +1  True
    i4         +1  True
    o1         -1  False
    o2         -1  False
    o3         +1  True
    o4         +1  True

Solution #3 (energy = -89.50):

    Name(s)  Spin  Boolean
    -------  ----  -------
    i1         +1  True
    i2         -1  False
    i3         +1  True
    i4         -1  False
    o1         -1  False
    o2         -1  False
    o3         +1  True
    o4         +1  True

Solution #4 (energy = -89.50):

    Name(s)  Spin  Boolean
    -------  ----  -------
    i1         -1  False
    i2         +1  True
    i3         -1  False
    i4         +1  True
    o1         -1  False
    o2         -1  False
    o3         +1  True
    o4         +1  True

Solution #5 (energy = -89.50):

    Name(s)  Spin  Boolean
    -------  ----  -------
    i1         +1  True
    i2         +1  True
    i3         -1  False
    i4         -1  False
    o1         -1  False
    o2         -1  False
    o3         +1  True
    o4         +1  True

Solution #6 (energy = -89.50):

    Name(s)  Spin  Boolean
    -------  ----  -------
    i1         -1  False
    i2         +1  True
    i3         +1  True
    i4         -1  False
    o1         -1  False
    o2         -1  False
    o3         +1  True
    o4         +1  True
```
One can even specify combinations of inputs and outputs.  As an exercise, see what solutions `--pin="i1 i2 i3 o2 := 0 1 1 0"` leads to.
