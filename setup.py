#! /usr/bin/env python

###################################
# Install QMASM                   #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from setuptools import setup, find_packages

setup(name = "qmasm",
      version = "4.0",
      description = "Quantum Macro Assembler",
      author = "Scott Pakin",
      author_email = "pakin@lanl.gov",
      classifiers = ["Topic :: Software Development :: Compilers"],
      url = "https://github.com/lanl/qmasm",
      license = "BSD",
      keywords = "quantum assembler d-wave",
      entry_points = {
          "console_scripts": ["qmasm = qmasm.__main__:main"]
      },
      scripts = ["scripts/qb2qmasm",
                 "extras/qmasm-gen-all-to-all",
                 "extras/qmasm-gen-chimera",
                 "extras/qmasm-gen-current"],
      packages = find_packages("src"),
      package_dir = {"": "src"},
      python_requires = ">= 3.8",
      install_requires = ["dwave-ocean-sdk"]
)
