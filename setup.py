#! /usr/bin/env python

###################################
# Install QASM                    #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import os.path
from setuptools import setup, find_packages
from setuptools.command.install import install as _install

class install(_install):
    def run(self):
        # Install, then remove qasm.py, keeping only qasm.
        _install.run(self)
        pyscript = os.path.join(self.install_scripts, "qasm.py")
        script = os.path.join(self.install_scripts, "qasm")
        os.rename(pyscript, script)

setup(name = "QASM",
      version = "1.0",
      description = "Quantum Macro Assembler",
      author = "Scott Pakin",
      author_email = "pakin@lanl.gov",
      classifiers = ["Topic :: Software Development :: Compilers"],
      url = "https://github.com/losalamos/qasm",
      license = "BSD",
      keywords = "quantum assembler d-wave",
      packages = find_packages(),
      scripts = ["qasm.py"],
      cmdclass = {"install": install}
)
