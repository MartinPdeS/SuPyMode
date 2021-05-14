#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys

try:
    from setuptools import setup, find_packages
except ImportError:
    print('Please install or upgrade setuptools or pip to continue')
    sys.exit(1)


requirements = [ 'numpy', 'scipy', 'matplotlib' ]

setup_dict = dict(
      description='Coupled mode modlisation for fiber optic coupler',
      name = 'SuPyMode',
      version = '1.1',
      author = 'Martin Poinsinet de Sivry',
      author_email = 'Martin.poinsinet.de.sivry@gmail.com',
      packages=find_packages(),
      data_file=['*.json'],
      py_modules = [],
      install_requires = ['pandas', 'scipy', 'numpy', 'matplotlib','shapely','descartes'],#requirements,
      license = 'Full private, no reproduction authorized',
      url='https://gitlab.com/PolyMtlLFO/SuPyModes',
      long_description=open('README.md').read(),
      platforms = ['Linux', 'Max OSX']
)

setup(**setup_dict)
