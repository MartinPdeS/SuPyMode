SuPyMode
========

|python|
|docs|
|Citation|
|Unittest|
|PyPi|
|PyPi_download|
|colab|


..  figure:: https://github.com/MartinPdeS/SuPyMode/blob/master/docs/images/mode_propagation.gif?raw=true
   :alt: some image
   :class: with-shadow float-left
   :width: 800px

   Propagation of mode in an adiabatic 2x1 modally-specific photonic lantern.




This project aims to develop an useful tool design and optimize fiber optic tapered component.
SuPyMode is a Python library linked to a c++ core allowing for a flexible interface and fast computing core.
The library also aims to offer the end-user a great vizual tools for data analysis.
To this day, SuPyMode as been proven a useful tool to develop very-short 2x1 and 3x1 modally specific photonic lantern with very low loss and cross-talk.

----

Documentation
**************
All the latest available documentation is available `here <https://supymodes.readthedocs.io/en/latest/>`_ or you can click the following badge:

|docs|


----


Installation
------------


Pip installation
****************

The package have been uploaded as wheel for a few OS (Linux, MacOS) and need Python 3.10.
As such, with the adequate configuration one can simply do

.. code-block:: python

   >>> pip3 install SuPyMode



Manual installation
*******************

To install manually (os independent) you will need to install:

1. cmake (3.16+)

In order to use the SuPyMode Simulator Library, one must have installed the python dependencies:

.. code-block:: python

    >>> pip3 install -r requirements.txt

Then, download and install the SuPyMode package:

.. code-block:: python

    >>> git clone https://github.com/MartinPdeS/SuPyMode.git
    >>> cd SuPyMode && mkdir build && cd build
    >>> cmake ..
    >>> make install                    (Linux, MacOs)
    >>> cd ..
    >>> pip3 install .

----


Contact Information
*******************

As of 2021 the project is still under development if you want to collaborate it would be a pleasure. I encourage you to contact me.

PyMieSim was written by `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_  .

Email:`martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=SuPyMode>`_ .


.. |python| image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
   :target: https://www.python.org/

.. |docs| image:: https://readthedocs.org/projects/supymodes/badge/?version=latest
   :target: https://supymodes.readthedocs.io/en/latest/
   :alt: Documentation Status

.. |Citation| image:: https://zenodo.org/badge/366930899.svg
   :target: https://zenodo.org/badge/latestdoi/366930899

.. |Unittest| image:: https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/MartinPdeS/2cac8ebb51cd9ce07dc7a955648301d5/raw/dbe91085da6e5b228bc9a0dcaeb4f50ae5abbe2d/coverage_badge.json

.. |PyPi| image:: https://badge.fury.io/py/SuPyMode.svg
   :target: https://pypi.org/project/SuPyMode/

.. |PyPi_download| image:: https://img.shields.io/pypi/dm/SuPyMode.svg
   :target: https://pypi.org/project/SuPyMode/

.. |colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :target: https://colab.research.google.com/github/MartinPdeS/SuPyMode/blob/master/SuPyModes.ipynb



