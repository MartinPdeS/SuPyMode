SuPyModes
==========


This project aims to produce an useful tool (python library) to simulate propagation mode in a very wide range of parameters and form.
It also offer the possiblity to compute the coupling coefficient between those mode and hence is a pratical tool to simulate couplers.


----

Documentation
**************

All the latest available documentation is available in the docs/build/index.html file


----


Packages depedencies
********************

In order to use the SuPyMode Simulator Library, one must have installed the following packages:


1. Numpy
2. Scipy
3. Pandas
4. Matplotlib
5. Shapely
6. Descartes

----

Using pip3 one can use the following commands:

.. code-block:: python

    pip3 install -r requirement.txt
    pip3 install -e ../SuPyModes


Or to install individually the packages:

.. code-block:: python

    >>> pip3 install Numpy
    >>> pip3 install Scipy
    >>> pip3 install Pandas
    >>> pip3 install Matplotlib
    >>> pip3 install Shapely
    >>> pip3 install Descartes
    >>> apt-get install python-sphinx (for Unix OS)
    >>> sudo port install py27-sphinx (for Mac OS)




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