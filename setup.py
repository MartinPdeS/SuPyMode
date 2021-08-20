#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
import numpy
from shutil import rmtree
import pathlib
from setuptools import setup, Extension
import subprocess

from setuptools import find_packages, setup, Command



# Package meta-data.
NAME            = 'SuPyMode'
DESCRIPTION     = 'A package for light propagation in fiber optics.'
URL             = 'https://github.com/MartinPdeS/SuPyModes'
EMAIL           = 'Martin.poinsinet.de.sivry@gmail.com'
AUTHOR          = 'Martin Poinsinet de Sivry',
REQUIRES_PYTHON = '>3.8.0'
VERSION         = '0.0.4'

# What packages are required for this module to be executed?
REQUIRED = ['scipy',
            'matplotlib',
            'cython',
            'pybind11',
            'vtk',
            'pandas',
            'mayavi',
            'coverage',
            'vtk',
            'numpy']


EXTRAS = {}

here = os.path.abspath(os.path.dirname(__file__))


macro = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')]


try:
    with io.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION


# Where the magic happens:
setup(
    name                          = NAME,
    version                       = about['__version__'],
    description                   = DESCRIPTION,
    long_description              = long_description,
    long_description_content_type = 'text/markdown',
    author                        = AUTHOR,
    author_email                  = EMAIL,
    setup_requires                = ['numpy'],
    python_requires               = '>=3.6',#REQUIRES_PYTHON,
    url                           = URL,
    packages                      = find_packages(exclude=["tests", "*.tests", "*.tests.*", "tests.*"]),
    install_requires              = REQUIRED,
    extras_require                = EXTRAS,
    dependency_links              = [],
    include_package_data          = True,
    ext_modules                   = None,
    license                       = 'MIT',
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: MIT License',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
    ],
    # $ setup.py publish support.
    #cmdclass={'upload': UploadCommand},
)
