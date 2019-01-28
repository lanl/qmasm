QMASM examples
=============

This directory contains examples of QMASM code.

Feature test
------------

* Main file: [`feature-test.qmasm`](feature-test.qmasm)

* Command line: `qmasm --format=qmasm feature-test.qmasm`

This program does nothing useful.  Rather, it provides a showcase of [all of QMASM's language constructs](//github.com/lanl/qmasm/wiki/File-format).  Read the [`feature-test.qmasm` source code](feature-test.qmasm) to see examples of the various mechanisms that a QMASM program can use.

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

Macro chaining
--------------

* Main file: [`and4.qmasm`](and4.qmasm)

* Command line: `qmasm --run and4.qmasm`

This is a demonstration of QMASM's multi-instantiation form of `!use_macro`.  The program defines a 2-input AND gate (*Y* = *A* AND *B*) called `and_chain` that chains its output (`Y`) to the next macro instance's `A` input.  The program then defines a macro called `big_and` that instantiates `and_chain` three times, as `$and1`, `$and2`, and `$and3`, then names its overall inputs (`A`, `B`, `C`, and `D`) and output (`Y`) in terms of its constituent macros' inputs and outputs.  The result is a 4-input AND gate.

```
Solution #1 (energy = -21.25, tally = 236):

    Name(s)       Spin  Boolean
    ------------  ----  --------
    big_and.A       -1  False
    big_and.B       -1  False
    big_and.C       -1  False
    big_and.D       -1  False
    big_and.Y       -1  False

Solution #2 (energy = -21.25, tally = 138):

    Name(s)       Spin  Boolean
    ------------  ----  --------
    big_and.A       -1  False
    big_and.B       -1  False
    big_and.C       -1  False
    big_and.D       +1  True
    big_and.Y       -1  False

Solution #3 (energy = -21.25, tally = 131):

    Name(s)       Spin  Boolean
    ------------  ----  --------
    big_and.A       -1  False
    big_and.B       -1  False
    big_and.C       +1  True
    big_and.D       -1  False
    big_and.Y       -1  False
```
<p style="text-align: center">⋮</p>

```
Solution #15 (energy = -21.25, tally = 7):

    Name(s)       Spin  Boolean
    ------------  ----  --------
    big_and.A       +1  True
    big_and.B       +1  True
    big_and.C       +1  True
    big_and.D       -1  False
    big_and.Y       -1  False

Solution #16 (energy = -21.25, tally = 3):

    Name(s)       Spin  Boolean
    ------------  ----  --------
    big_and.A       +1  True
    big_and.B       +1  True
    big_and.C       +1  True
    big_and.D       +1  True
    big_and.Y       +1  True
```

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

* Command line: `qmasm --run --pin="in[1:4] := 0110" sort4.qmasm`

Sort a list of four 1-bit numbers.  `sort4.qmasm` implements a 4-element sorting network from [*Sorting and Searching*](http://www.informit.com/store/art-of-computer-programming-volume-3-sorting-and-searching-9780201896855).  Specify values for inputs *in[1]*, *in[2]*, *in[3]*, and *in[4]*, and the program will sort these into *out[1]*, *out[2]*, *out[3]*, and *out[4]*:
```
Solution #1 (energy = -122.00, tally = 3):

    Name(s)    Spin  Boolean
    ---------  ----  --------
    in[1]        -1  False
    in[2]        +1  True
    in[3]        +1  True
    in[4]        -1  False
    out[1]       -1  False
    out[2]       -1  False
    out[3]       -1  False
    out[4]       +1  True
```

Oops, that didn't work.  The D-Wave incorrectly sorted `0 1 1 0` into `0 0 0 1`, due to hardware artifacts such as limited precision or cross-qubit interference.  This is a case where optimization postprocessing is useful for transforming nearly correct solutions into truly correct solutions:
```
$ qmasm --postproc=opt --run --pin="in[1:4] := 0110" sort4.qmasm
Solution #1 (energy = -93.50, tally = 1000):

    Name(s)    Spin  Boolean
    ---------  ----  --------
    in[1]        -1  False
    in[2]        +1  True
    in[3]        +1  True
    in[4]        -1  False
    out[1]       -1  False
    out[2]       -1  False
    out[3]       +1  True
    out[4]       +1  True
```
Because there is no clear distinction between inputs and outputs, one can also specify the outputs and receive a list of inputs that would sort to those outputs.  For example, `qmasm --run --pin="out[1:4] := 0011" sort4.qmasm` produces the following six solutions:
```
Solution #1 (energy = -92.50, tally = 63):

    Name(s)    Spin  Boolean
    ---------  ----  --------
    in[1]        +1  True
    in[2]        -1  False
    in[3]        -1  False
    in[4]        +1  True
    out[1]       -1  False
    out[2]       -1  False
    out[3]       +1  True
    out[4]       +1  True

Solution #2 (energy = -92.50, tally = 34):

    Name(s)    Spin  Boolean
    ---------  ----  --------
    in[1]        -1  False
    in[2]        -1  False
    in[3]        +1  True
    in[4]        +1  True
    out[1]       -1  False
    out[2]       -1  False
    out[3]       +1  True
    out[4]       +1  True

Solution #3 (energy = -92.50, tally = 128):

    Name(s)    Spin  Boolean
    ---------  ----  --------
    in[1]        +1  True
    in[2]        -1  False
    in[3]        +1  True
    in[4]        -1  False
    out[1]       -1  False
    out[2]       -1  False
    out[3]       +1  True
    out[4]       +1  True

Solution #4 (energy = -92.50, tally = 16):

    Name(s)    Spin  Boolean
    ---------  ----  --------
    in[1]        -1  False
    in[2]        +1  True
    in[3]        -1  False
    in[4]        +1  True
    out[1]       -1  False
    out[2]       -1  False
    out[3]       +1  True
    out[4]       +1  True

Solution #5 (energy = -92.50, tally = 56):

    Name(s)    Spin  Boolean
    ---------  ----  --------
    in[1]        +1  True
    in[2]        +1  True
    in[3]        -1  False
    in[4]        -1  False
    out[1]       -1  False
    out[2]       -1  False
    out[3]       +1  True
    out[4]       +1  True

Solution #6 (energy = -92.50, tally = 25):

    Name(s)    Spin  Boolean
    ---------  ----  --------
    in[1]        -1  False
    in[2]        +1  True
    in[3]        +1  True
    in[4]        -1  False
    out[1]       -1  False
    out[2]       -1  False
    out[3]       +1  True
    out[4]       +1  True
```
That is, it permuted {0, 0, 1, 1} into {1, 0, 0, 1}, {0, 0, 1, 1}, {1, 0, 1, 0}, {0, 1, 0, 1}, {1, 1, 0, 0}, and {0, 1, 1, 0}.

One can even specify combinations of inputs and outputs.  As an exercise, see what solutions `qmasm --run --pin="in[1:3] out[2] := 0110" sort4.qmasm` leads to.

Shortest path through a maze
----------------------------

* Main file: [`maze6x6.qmasm`](maze6x6.qmasm)

* Command line: `qmasm --postproc=opt --run maze6x6.qmasm | egrep 'Solution|True'`

Find the shortest path through a 6×6 maze.
```
   A  B  C  D  E  F
  +--+--+--+--+--+  +
1 |           |     |
  +  +  +--+  +  +  +
2 |  |     |  |  |  |
  +--+--+  +  +  +--+
3 |  |     |        |
  +  +--+--+  +  +  +
4 |  |     |  |  |  |
  +  +--+  +--+  +--+
5 |                 |
  +  +  +  +--+  +--+
6 |  |  |  |        |
  +--+  +--+--+--+--+
```
The cleverness in the implementation is that each room of the maze (macro `room`) constrains the shortest path to traversing either zero or two of the four compass directions.  The former case implies that the shortest path does not pass through that room.  The latter case implies that the shortest path enters and exits the room exactly once.  All other options require more energy.  Pinning the ingress and egress to *true* (done within `maze6x6.qmasm` itself) causes the minimal-energy solution to represent the shortest path between the corresponding two rooms.

The program takes a long time to embed in the Chimera graph so be patient.  (`maze6x6.qmasm` can therefore serve as a good test for new embedding algorithms.)  When it runs, it finds a single valid solution:
```
Solution #1 (energy = -664.00, tally = 1000):
    B5.E       +1  True
    B5.S       +1  True
    B6.N       +1  True
    B6.S       +1  True
    C5.E       +1  True
    C5.W       +1  True
    D5.E       +1  True
    D5.W       +1  True
    E1.E       +1  True
    E1.S       +1  True
    E2.N       +1  True
    E2.S       +1  True
    E3.N       +1  True
    E3.S       +1  True
    E4.N       +1  True
    E4.S       +1  True
    E5.N       +1  True
    E5.W       +1  True
    F1.N       +1  True
    F1.W       +1  True
```
This can be validated by beginning at the ingress at `F1.N`.  The only exit from F1 other than the ingress is `F1.W`, which takes us to `E1`.  Because we came from `E1.E` the only other exit is `E1.S`, which takes us to `E2`.  The process continues until we reach the egress at `B6.S`.

If you want to experiment with other mazes, [`qmasm-maze.go`](qmasm-maze.go) is the source code for a maze-generating program written in [Go](https://golang.org/).  Once you've installed a Go compiler, build `qmasm-maze` with
```bash
go get github.com/spakin/disjoint
go build qmasm-maze.go
```
In case you're unfamiliar with Go, the first line of the above installs a dependency, and the second line compiles the program.  You'll need to set your `GOPATH` environment variable so `go` knows where to install to (e.g., `export GOPATH=$HOME/go`).

Once compiled, `qmasm-maze` accepts a pair of dimensions (the number of rooms wide and tall) to use for the maze and outputs QMASM code:
```bash
./qmasm-maze gen 6 6 > my-maze.qmasm
```

`qmasm-maze` can also validate solutions produced by running the code.  For example, the solution shown above can be seen to be correct:
```
$ qmasm --postproc=opt --run maze6x6.qmasm | ./qmasm-maze validate
Solution 1: F1 -- E1 -- E2 -- E3 -- E4 -- E5 -- D5 -- C5 -- B5 -- B6
```
Without optimization postprocessing, QMASM frequently returns invalid paths through the maze:
```
$ ./qmasm-maze gen 6 6 | qmasm --run 2>&1 | ./qmasm-maze validate
Solution 1: [Room F1 contains 1 exit(s) but should have 0 or 2]
Solution 2: [Room E1 contains 1 exit(s) but should have 0 or 2]
```
