#! /usr/bin/env python

###################################
# Install QMASM                   #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import os.path
from setuptools import setup, find_packages
from setuptools.command.install import install as _install

script_list = ["qmasm", "qb2qmasm", "qmasm-ground-state", "qmasm-qbsolv"]

class install(_install):
    def run(self):
        _install.run(self)

setup(name = "QMASM",
      version = "3.0",
      description = "Quantum Macro Assembler",
      author = "Scott Pakin",
      author_email = "pakin@lanl.gov",
      classifiers = ["Topic :: Software Development :: Compilers"],
      url = "https://github.com/lanl/qmasm",
      license = "BSD",
      keywords = "quantum assembler d-wave",
      packages = find_packages(),
      scripts = ["scripts/"+s for s in script_list],
      cmdclass = {"install": install}
)
