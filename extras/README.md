QMASM extras
=============

This directory contains additional tools that may be of use to QMASM programmers.

QMASM stylesheets
----------------

[`qmasm.ssh`](qmasm.ssh) is a stylsheet for the [a2ps](https://www.gnu.org/software/a2ps/) "Any to PostScript" pretty printer.  One can print QMASM source code with a command like the following:
```bash
a2ps -1 --prologue=color -E/path/to/qmasm.ssh my-program.qmasm | lpr
```
See the a2ps documentation for instructions on how to install `qmasm.ssh` and use that stylesheet implicitly for all source files ending in `.qmasm`.  (In short, put `qmasm.ssh` in the `a2ps/sheets` directory, and edit `a2ps/sheets/sheets.map` to associate the `.qmasm` file extension with the `qmasm.ssh` stylesheet.)
