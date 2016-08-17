QASM examples
=============

This directory contains examples of QASM code.

One of five
-----------

* Main file: [`1of5.qasm`](1of5.qasm)

* Command line: `qasm --run --chain=-8 1of5.qasm`

This is a trivial demonstration of QASM.  The program defines variables *A*, *B*, *C*, *D*, and *E* and outputs all (Boolean) values of those in which exactly one variable is *true*.  What makes this program interesting is that the default chain strength—in this case, computed to be -4—is insufficiently strong to produce all five solutions with equal probability.  We therefore need to double the chain strength on the command line to get the desired result:
```
Solution #1 (energy = -18.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          -1  False
    B          -1  False
    C          -1  False
    D          -1  False
    E          +1  True

Solution #2 (energy = -18.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          -1  False
    B          -1  False
    C          -1  False
    D          +1  True
    E          -1  False

Solution #3 (energy = -18.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          -1  False
    B          -1  False
    C          +1  True
    D          -1  False
    E          -1  False

Solution #4 (energy = -18.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          -1  False
    B          +1  True
    C          -1  False
    D          -1  False
    E          -1  False

Solution #5 (energy = -18.12):

    Name(s)  Spin  Boolean
    -------  ----  -------
    A          +1  True
    B          -1  False
    C          -1  False
    D          -1  False
    E          -1  False
```

Try experimenting with the `--pin` option.  If a variable is pinned to *true*, QASM will output the single solution that honors that constraint.  If a variable is pinned to *false*, QASM will output the four solutions.  If *two* variables are pinned to *true*, a situation `1of5.qasm` does not include in its ground state, QASM may return no solutions or it may return one or more *incorrect* solutions—whatever exhibits the lowest total energy and doesn't break any chains or pins.

Circuit satisfiability
----------------------

* Main file: [`circsat.qasm`](circsat.qasm)

* Helper file: [`gates.qasm`](gates.qasm)

* Command line: `qasm --run --pin="x10 := true" circsat.qasm`

Given a Boolean expression with *n* inputs and one output, determine the sets of inputs that make the output *true*.  This is a classic NP-complete problem.  `circsat.qasm` represents a particular 3-input Boolean expression borrowed from the [*Introduction to Algorithms*](https://mitpress.mit.edu/books/introduction-algorithms) textbook's discussion of NP-completeness.

`gates.qasm` defines macros for various Boolean operators (NOT, 2-input OR, 3-input AND), which are then used by the top-level program, `circsat.qasm`.  `circsat.qasm` maps inputs *x1*, *x2*, and *x3* to output *x10*.  By default, the program returns all sets of inputs and the corresponding output.  Pinning *x10* to *true* returns only the solutions to the circuit-satisfiability problem.  In this case, there is only one solution:
```
Solution #1 (energy = -48.50):

    Name(s)  Spin  Boolean
    -------  ----  -------
    x1         +1  True
    x10        +1  True
    x2         +1  True
    x3         -1  False
```
