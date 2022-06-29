#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
import numpy
import pathlib
import subprocess
import pkg_resources

from shutil     import rmtree
from setuptools import setup, Extension
from setuptools import find_packages, setup, Command



# Package meta-data.
NAME            = 'SuPyMode'
DESCRIPTION     = 'A package for light propagation in fiber optics.'
URL             = 'https://github.com/MartinPdeS/SuPyMode'
EMAIL           = 'Martin.poinsinet.de.sivry@gmail.com'
AUTHOR          = 'Martin Poinsinet de Sivry',
REQUIRES_PYTHON = '>3.8.0'


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

with open(os.path.join(__location__, 'Version.py'), "r+") as f:
    Version = f.read().rstrip("\n").split(".")
    Major, Mid, Minor = int(Version[0]), int(Version[1]), int(Version[2])

if '--NewMajor' in sys.argv:
    Major += 1
    sys.argv.remove('--NewMajor')
if '--NewMid' in sys.argv:
    Mid += 1
    sys.argv.remove('--NewMidr')
if '--NewMinor' in sys.argv:
    Minor += 1
    sys.argv.remove('--NewMinor')

Version = f'{Major}.{Mid}.{Minor}'

print(f"SuPyMode Version: {Version}")

with open(os.path.join(__location__, 'Version.py'), "w+") as f:
    f.writelines(Version)


# What packages are required for this module to be executed?
requirementPath = os.path.join(os.path.dirname(__file__), 'requirements.txt')

with open(requirementPath,'r') as requirements_txt:
    REQUIRED = [
        str(requirement)
        for requirement
        in pkg_resources.parse_requirements(requirements_txt)
    ]


EXTRAS = {}

here = os.path.abspath(os.path.dirname(__file__))


macro = [('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_Version')]


try:
    with io.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


about = {}
if not Version:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(os.path.join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = Version


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
)
