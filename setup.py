#! /usr/bin/env python

###################################
# Install QMASM                   #
# By Scott Pakin <pakin@lanl.gov> #
###################################

import os.path
from setuptools import setup, find_packages
from setuptools.command.install import install as _install

class install(_install):
    def run(self):
        # Install, then remove qmasm.py, keeping only qmasm.
        _install.run(self)
        for scr in ["qmasm", "qb2qmasm"]:
            pyscript = os.path.join(self.install_scripts, scr + ".py")
            script = os.path.join(self.install_scripts, scr)
            os.rename(pyscript, script)

setup(name = "QMASM",
      version = "1.2",
      description = "Quantum Macro Assembler",
      author = "Scott Pakin",
      author_email = "pakin@lanl.gov",
      classifiers = ["Topic :: Software Development :: Compilers"],
      url = "https://github.com/losalamos/qmasm",
      license = "BSD",
      keywords = "quantum assembler d-wave",
      packages = find_packages(),
      scripts = ["qmasm.py", "qb2qmasm.py"],
      cmdclass = {"install": install}
)
