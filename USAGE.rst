SuPyMode usage guide
====================

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

Citing SuPyMode
---------------

If you use SuPyMode in academic work, please cite:

.. code-block:: text

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
