SuPyMode
========


This project aims to produce an useful tool (python library) to simulate propagation mode in a very wide range of parameters and form.
It also offer the possiblity to compute the coupling coefficient between those mode and hence is a pratical tool to simulate couplers.


----

Documentation
**************

All the latest available documentation is available in the docs/build/index.html file


----


Installation
------------

Manual installation
*******************

To install manually (os independent) you will need to install:
    - cmake (3.0+)
    - Boost (1.58+) 
    - Boost component: iostream filesystem

Once it's done, do the following:

.. code-block::
    >>> git clone https://github.com/MartinPdeS/SuPyModes.git
    >>> cd SuPyModes && mkdir build && cd build
    >>> cmake ..
    >>> make install                    (Linux)
    >>> msbuild INSTALL.vcxproj         (from visual Studio powershell)
.

----

Packages depedencies
********************

In order to use the SuPyMode Simulator Library, one must have installed the following packages:


1. Numpy
2. Scipy
4. Matplotlib
5. Shapely
6. Descartes

----

Or to install individually the packages:

.. code-block:: python
    >>> pip3 install Numpy
    >>> pip3 install Scipy
    >>> pip3 install Pandas
    >>> pip3 install Matplotlib
    >>> pip3 install Shapely
    >>> pip3 install Descartes
    >>> git clone https://github.com/MartinPdeS/SuPyMode.git

    >>> sudo apt-get install gnuplot
    >>> sudo apt-get install libgnuplot-iostream-dev

    >>> cd SuPyMode && git submodule init && git submodule update
    >>> cd extern/eigen && mkdir build && cd build && cmake .. && make install && cd ..
    >>> cd extern/spectra && mkdir build && cd build && cmake .. && make install && cd ..
    >>> cmake .
    >>> make EigenSolver
    >>> python3 setup.py install -v




IMPORTANT NOTICE: All units in the simulator are micrometers!

----


Sellmeier
*********

Refractive index for some materials are computed via Sellmeier equation:

.. code-block:: python

    Fused_silica(wavelength = 1.55)
    Ambiant_air(wavelength = 1.55)
    BK7_glass(wavelength = 1.55)


The 1.55 stand for a wavelength of 1550 nm.

----


Contact Information
************************
As of 2021 the project is still under development if you want to collaborate it would be a pleasure. I encourage you to contact me.

PyMieSim was written by `Martin Poinsinet de Sivry-Houle <https://github.com/MartinPdS>`_  .

Email:`martin.poinsinet-de-sivry@polymtl.ca <mailto:martin.poinsinet-de-sivry@polymtl.ca?subject=PyMieSim>`_ .
