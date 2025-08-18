.. figure:: https://github.com/MartinPdeS/SuPyMode/blob/master/docs/images/mode_propagation.gif?raw=true
   :alt: Propagation of mode in an adiabatic 2x1 modally-specific photonic lantern.
   :width: 800px
   class: with-shadow float-left

.. list-table::
   :widths: 10 25 25 25
   :header-rows: 0

   * - Meta
     - |python|
     - |docs|
     - |zenodo|
   * - Testing
     - |ci/cd|
     - |coverage|
     - |colab|
   * - PyPI
     - |PyPI|
     - |PyPI_download|
     -
   * - Anaconda
     - |anaconda|
     - |anaconda_download|
     - |anaconda_date|

SuPyMode
========

This project aims to develop an useful tool design and optimize fiber optic tapered component.
SuPyMode is a Python library linked to a c++ core allowing for a flexible interface and fast computing core.
The library also aims to offer the end-user a great vizual tools for data analysis.
To this day, SuPyMode as been proven a useful tool to develop very-short 2x1 and 3x1 modally specific photonic lantern with very low loss and cross-talk.

Features
--------
- Fast and efficient simulation of fiber optic tapered components.
- User-friendly interface for easy integration into existing workflows.
- Comprehensive visualization tools for data analysis and interpretation.

Installation
------------
**SuPyMode** is available on PyPI and Anaconda.  Install it with:

.. code-block:: bash

   pip install SuPyMode
   conda install SuPyMode

See the `online documentation <https://supymode.readthedocs.io/>`_ for detailed
usage and additional examples.


Quick example
-------------
Below is a short example computing the mode propgation in a simple fiber.

.. code-block:: python

   from SuPyMode.workflow import Workflow, fiber_loader, Boundaries, BoundaryValue, DomainAlignment

   wavelength = 1550e-9

   fiber = fiber_loader.load_fiber('SMF28', clad_refractive_index=1.4444, remove_cladding=False)

   boundaries = [
      Boundaries(right=BoundaryValue.SYMMETRIC, top=BoundaryValue.SYMMETRIC),
      Boundaries(right=BoundaryValue.SYMMETRIC, top=BoundaryValue.ANTI_SYMMETRIC)
   ]


   workflow = Workflow(
      fiber_list=[fiber],             # List of fiber to be added in the mesh, the order matters.
      wavelength=wavelength,          # Wavelength used for the mode computation.
      resolution=80,                  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
      x_bounds=DomainAlignment.LEFT,  # Mesh x-boundary structure.
      y_bounds=DomainAlignment.BOTTOM,# Mesh y-boundary structure.
      boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
      n_sorted_mode=3,                # Total computed and sorted mode.
      n_added_mode=2,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
      plot_geometry=True,             # Plot the geometry mesh before computation.
      auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
      itr_final=0.05,                 # Final value of inverse taper ratio to simulate
      index_scrambling=0              # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
   )

   workflow.superset.plot(plot_type='field', itr_list=[1.0, 0.1])

   workflow.superset.plot(plot_type='index')

   workflow.superset.plot(plot_type='normalized-coupling')

   workflow.superset.plot(plot_type='adiabatic')


Building from source
--------------------
For development or manual compilation, clone the repository and run:

.. code-block:: bash

   git submodule update --init
   mkdir build && cd build
   cmake ../ -G"Unix Makefiles"
   sudo make install
   cd ..
   python -m pip install .

Testing
-------
Run the unit tests with:

.. code-block:: bash

   pip install SuPyMode[testing]
   pytest

Citing SuPyMode
---------------
If you use SuPyMode in academic work, please cite:

.. code-block:: none

   @article{de2024supymode,
      title={SuPyMode: an open-source library for design and optimization of fiber optic components},
      author={de Sivry-Houle, Martin Poinsinet and Becerra Deana, Rodrigo Itzamna and Virally, St{\'e}phane and Godbout, Nicolas and Boudoux, Caroline},
      journal={Optics Continuum},
      volume={3},
      number={2},
      pages={242--255},
      year={2024},
      publisher={Optica Publishing Group}
   }



Contact
-------
For questions or contributions, contact `martin.poinsinet.de.sivry@gmail.com <mailto:martin.poinsinet.de.sivry@gmail.com>`_.

.. |python| image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
    :alt: Python
    :target: https://www.python.org/
.. |zenodo| image:: https://zenodo.org/badge/366930899.svg
   :target: https://zenodo.org/badge/latestdoi/366930899
   :alt: Scientific article
.. |colab| image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Google Colab
   :target: https://colab.research.google.com/github/MartinPdeS/SuPyMode/blob/master/notebook.ipynb
.. |docs| image:: https://github.com/martinpdes/supymode/actions/workflows/deploy_documentation.yml/badge.svg
   :target: https://martinpdes.github.io/SuPyMode/
   :alt: Documentation Status
.. |PyPi| image:: https://badge.fury.io/py/SuPyMode.svg
   :alt: PyPI version
   :target: https://pypi.org/project/SuPyMode/
.. |PyPi_download| image:: https://img.shields.io/pypi/dm/supymode.svg
   :alt: PyPI downloads
   :target: https://pypistats.org/packages/supymode
.. |coverage| image:: https://raw.githubusercontent.com/MartinPdeS/SuPyMode/python-coverage-comment-action-data/badge.svg
   :alt: Unittest coverage
   :target: https://htmlpreview.github.io/?https://github.com/MartinPdeS/SuPyMode/blob/python-coverage-comment-action-data/htmlcov/index.html
.. |ci/cd| image:: https://github.com/martinpdes/supymode/actions/workflows/deploy_coverage.yml/badge.svg
    :alt: Unittest Status
.. |anaconda| image:: https://anaconda.org/martinpdes/supymode/badges/version.svg
    :alt: Anaconda version
    :target: https://anaconda.org/martinpdes/supymode
.. |anaconda_download| image:: https://anaconda.org/martinpdes/supymode/badges/downloads.svg
    :alt: Anaconda downloads
    :target: https://anaconda.org/martinpdes/supymode
.. |anaconda_date| image:: https://anaconda.org/martinpdes/supymode/badges/latest_release_relative_date.svg
    :alt: Latest release date
    :target: https://anaconda.org/martinpdes/supymode




