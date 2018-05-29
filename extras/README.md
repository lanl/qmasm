QMASM extras
============

This directory contains additional tools that may be of use to QMASM programmers.

QMASM stylesheets
-----------------

[`qmasm-mode.el`](qmasm-mode.el) is an Emacs major mode for editing QMASM source code.  It is currently a very simple mode that provides only syntax highlighting.  Load the mode manually with `M-x load-library` then `M-x qmasm-mode` or automatically by including statements like the following in your `.emacs` file:
```Emacs Lisp
(load-library "qmasm-mode")
(setq auto-mode-alist (cons '("\\.qmasm$" . qmasm-mode) auto-mode-alist))
```

[`qmasm.ssh`](qmasm.ssh) is a stylsheet for the [a2ps](https://www.gnu.org/software/a2ps/) "Any to PostScript" pretty printer.  One can print QMASM source code with a command like the following:
```bash
a2ps -1 --prologue=color -E/path/to/qmasm.ssh my-program.qmasm | lpr
```
See the a2ps documentation for instructions on how to install `qmasm.ssh` and use that stylesheet implicitly for all source files ending in `.qmasm`.  (In short, put `qmasm.ssh` in the `a2ps/sheets` directory, and edit `a2ps/sheets/sheets.map` to associate the `.qmasm` file extension with the `qmasm.ssh` stylesheet.)

Topology generation
-------------------

QMASM's `--topology-file` option lets the user define a graph topology to target in place of the D-Wave hardware's actual topology.  The format is a list of space-separated vertex pairs, one pair per line.  Comments, which go from the first `#` character to the end of the line, can also be included in the file.

The following scripts can be used to construct files that can be used as an argument to `--topology-file`:

* [`qmasm-gen-chimera`](qmasm-gen-chimera) generates a complete Chimera graph of arbitrary size.  It takes three arguments: the width of the Chimera graph in unit cells, the height of the Chimera graph in unit cells, and the number of vertices in each of a unit cell's two partitions.  For example, a complete D-Wave 2000Q could be generated with `qmasm-gen-chimera 16 16 4`.

* [`qmasm-gen-current`](qmasm-gen-current) outputs the current topology.  It takes no arguments but expects the various `DW_INTERNAL__*` environment variables to be set properly.

* [`qmasm-gen-all-to-all`](qmasm-gen-all-to-all) outputs a complete graph of a given number of vertices.  The script takes one argument, which is the number of vertices.
