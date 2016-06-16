#! /usr/bin/env python

###################################
# Install QASM                    #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from setuptools import setup

setup(name = "QASM",
      version = "1.0",
      description = "Quantum Macro Assembler",
      author = "Scott Pakin",
      author_email = "pakin@lanl.gov",
      url = "https://github.com/losalamos/qasm",
      scripts = ["qasm"],
      license = "BSD",
      keywords = "quantum assembler d-wave",
)
