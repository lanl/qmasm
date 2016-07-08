QASM: A Quantum Macro Assembler
===============================

Description
-----------

[D-Wave Systems, Inc.](http://www.dwavesys.com/) is a vendor of adiabatic quantum computers or _quantum annealers_.  These are special-purpose devices that rapidly find the spins in an [Ising spin system](https://en.wikipedia.org/wiki/Ising_model) that minimize total energy.  D-Wave provides various [software interfaces](http://www.dwavesys.com/software) for programming their hardware, but, to date, these are either too high-level (requiring a lot of classical pre- and/or post-processing) or too low-level (requiring manual layout of a problem on the physical topology).

QASM fills a gap in the D-Wave's software ecosystem by shielding the programmer from having to know system-specific hardware details while still enabling programs to be expressed at a fairly low level of abstraction.  It is therefore analogous to a conventional macro assembler and can be used in much the same way: as a target either for programmers who want a great deal of control over the hardware or for compilers that implement higher-level languages.

Despite having the same name, Los Alamos's QASM has nothing to do with [MIT's QASM](http://www.media.mit.edu/quanta/quanta-web/projects/qasm-tools/), which is used to describe quantum circuits (a different model of quantum computation from quantum annealing).

Installation
------------

QASM is written in Python and uses [Setuptools](https://pythonhosted.org/an_example_pypi_project/setuptools.html) for installation.  Use
```Python
python setup.py install
```
to install in the default location and
```Python
python setup.py install --prefix=/my/install/directory
```
to install elsewhere.

License
-------

QASM is provided under a BSD-ish license with a "modifications must be indicated" clause.  See [the LICENSE file](http://github.com/losalamos/qasm/blob/master/LICENSE.md) for the full text.

This package is part of the Hybrid Quantum-Classical Computing suite, known internally as LA-CC-16-032.

Author
------

Scott Pakin, <pakin@lanl.gov>
