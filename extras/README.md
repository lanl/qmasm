QASM extras
=============

This directory contains additional tools that may be of use to QASM programmers.

QASM stylesheets
----------------

[`qasm.ssh`](qasm.ssh) is a stylsheet for the [a2ps](https://www.gnu.org/software/a2ps/) "Any to PostScript" pretty printer.  One can print QASM source code with a command like the following:
```bash
a2ps -1 --prologue=color -E/path/to/qasm.ssh my-program.qasm | lpr
```
See the a2ps documentation for instructions on how to install `qasm.ssh` and use that stylesheet implicitly for all source files ending in `.qasm`.  (In short, put `qasm.ssh` in the `a2ps/sheets` directory, and edit `a2ps/sheets/sheets.map` to associate the `.qasm` file extension with the `qasm.ssh` stylesheet.)
