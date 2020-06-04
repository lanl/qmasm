#! /usr/bin/env python

###################################
# Install QMASM                   #
# By Scott Pakin <pakin@lanl.gov> #
###################################

from setuptools import setup, find_packages

long_description = '''QMASM fills a gap in the software ecosystem for [D-Wave's adiabatic quantum computers](http://www.dwavesys.com/) by shielding the programmer from having to know system-specific hardware details while still enabling programs to be expressed at a fairly low level of abstraction.  It is therefore analogous to a conventional macro assembler and can be used in much the same way: as a target either for programmers who want a great deal of control over the hardware or for compilers that implement higher-level languages.'''

setup(name = 'qmasm',
      version = '4.0.2',
      description = 'Quantum Macro Assembler',
      long_description = long_description,
      long_description_content_type = 'text/markdown',
      author = 'Scott Pakin',
      author_email = 'pakin@lanl.gov',
      classifiers = [
          'Topic :: Software Development :: Compilers',
          'Development Status :: 5 - Production/Stable',
          'Environment :: Console',
          'Programming Language :: Python :: 3',
          'Intended Audience :: Developers'],
      url = 'https://github.com/lanl/qmasm',
      download_url = 'https://github.com/lanl/qmasm/archive/v4.0.2.tar.gz',
      license = 'BSD-ish',
      keywords = ['quantum', 'annealing', 'macro', 'assembler', 'd-wave'],
      entry_points = {
          'console_scripts': ['qmasm = qmasm.__main__:main']
      },
      scripts = ['scripts/qb2qmasm',
                 'extras/qmasm-gen-all-to-all',
                 'extras/qmasm-gen-chimera',
                 'extras/qmasm-gen-current'],
      packages = find_packages('src'),
      package_dir = {'': 'src'},
      python_requires = '>= 3.8',
      install_requires = ['dwave-ocean-sdk']
)
