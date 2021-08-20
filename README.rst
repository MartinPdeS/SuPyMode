SuPyMode
========

|python|
|docs|

This project aims to produce an useful tool (python library) to simulate propagation mode in a very wide range of parameters and form.
It also offer the possiblity to compute the coupling coefficient between those mode and hence is a pratical tool to simulate couplers.


----

Documentation
**************
All the latest available documentation is available `here <https://supymode.readthedocs.io/en/latest/>`_ or you can click the following badge:

|docs|


----


Installation
------------


Pip installation
****************

The package have been uploaded as wheel for a few OS (Linux, Windows) and a few Python version (3.6, 3.8).
As such, with the adequate configuration one can simply do

.. code-block:: python

   >>> pip3 install SuPyMode



Manual installation
*******************

To install manually (os independent) you will need to install:

1. cmake (3.0+)
2. Boost (1.58+)
3. Boost component: iostream filesystem

In order to use the SuPyMode Simulator Library, one must have installed the python dependencies:

.. code-block:: python

    >>> pip3 install Numpy
    >>> pip3 install Scipy
    >>> pip3 install Pandas
    >>> pip3 install Matplotlib
    >>> pip3 install Shapely
    >>> pip3 install Descartes

Then, download and install the SuPyMode package:

.. code-block:: python

    >>> git clone https://github.com/MartinPdeS/SuPyModes.git
    >>> cd SuPyModes && mkdir build && cd build
    >>> cmake ..
    >>> make install                    (Linux)
    >>> msbuild INSTALL.vcxproj         (from visual Studio powershell)
    >>> cd ..
    >>> pip3 install .


----


Sellmeier
*********

Refractive index for some materials are computed via Sellmeier equation:

.. code-block:: python

    Fused_silica(wavelength = 1.55)
    Ambiant_air(wavelength = 1.55)
    BK7_glass(wavelength = 1.55)


IMPORTANT NOTICE: All units in the simulator are micrometers!
The 1.55 stand for a wavelength of 1550 nm.


----


Contact Information
************************
As of 2021 the project is still under development if you want to collaborate it would be a pleasure. I encourage you to contact me.

PyMieSim was written by `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_  .

Email:`martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=PyMieSim>`_ .


.. |python| image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
   :target: https://www.python.org/

.. |docs| image:: https://readthedocs.org/projects/supymode/badge/?version=latest
   :target: https://supymode.readthedocs.io/en/latest/?badge=latest
