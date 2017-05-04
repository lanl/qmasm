QMASM: A Quantum Macro Assembler
================================

[![Build Status](https://travis-ci.org/lanl/qmasm.svg?branch=master)](https://travis-ci.org/lanl/qmasm)

Description
-----------

QMASM fills a gap in the software ecosystem for [D-Wave's adiabatic quantum computers](http://www.dwavesys.com/) by shielding the programmer from having to know system-specific hardware details while still enabling programs to be expressed at a fairly low level of abstraction.  It is therefore analogous to a conventional macro assembler and can be used in much the same way: as a target either for programmers who want a great deal of control over the hardware or for compilers that implement higher-level languages.

N.B. This tool used to be called "QASM" but was renamed to avoid confusion with [MIT's QASM](http://www.media.mit.edu/quanta/quanta-web/projects/qasm-tools/), which is used to describe quantum circuits (a different model of quantum computation from what the D-Wave uses) and the [IBM Quantum Experience](http://www.research.ibm.com/quantum/)'s [QASM language](https://github.com/IBMQuantum/QASM), also used for describing quantum circuits.

Installation
------------

QMASM is written in Python and uses [Setuptools](https://pythonhosted.org/an_example_pypi_project/setuptools.html) for installation.  Use
```bash
python setup.py install
```
to install in the default location and
```bash
python setup.py install --prefix=/my/install/directory
```
to install elsewhere.

Documentation
-------------

Documentation for QMASM can be found on the [QMASM wiki](https://github.com/lanl/qmasm/wiki).

QMASM (then known as QASM) is discussed in the following publication:

> Scott Pakin. "A Quantum Macro Assembler". In _Proceedings of the 20th Annual IEEE High Performance Extreme Computing Conference_ (HPEC 2016), Waltham, Massachusetts, USA, 13–15 September 2016.  DOI: [10.1109/HPEC.2016.7761637](http://dx.doi.org/10.1109/HPEC.2016.7761637).


License
-------

QMASM is provided under a BSD-ish license with a "modifications must be indicated" clause.  See [the LICENSE file](http://github.com/lanl/qmasm/blob/master/LICENSE.md) for the full text.

This package is part of the Hybrid Quantum-Classical Computing suite, known internally as LA-CC-16-032.

Author
------

Scott Pakin, <pakin@lanl.gov>
