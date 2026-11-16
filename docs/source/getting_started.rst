Getting started
===============

Install and verify
------------------

.. code-block:: console

   $ python -m pip install SuPyMode
   $ python -c "import SuPyMode; print(SuPyMode.__version__)"

Units and validation
--------------------

SuPyMode shares TypedUnit's registry, available as ``SuPyMode.ureg``.
Wavelengths, coordinate vectors, grid spacings, taper radii, and segment lengths
accept explicit compatible quantities. Existing bare values continue to mean
SI metres. Workflow rotation angles accept explicit angular quantities or
legacy bare degrees.

TypedUnit validates inputs before conversion. Native model data, taper arrays,
and numerical outputs retain their existing SI conventions. Supplied quantity
objects and their display units are preserved; the conversion happens at the
numerical or external geometry boundary.

Refractive indices, taper ratios, alpha parameters, and amplitudes require plain
numbers. ``NoUnit`` contracts reject Pint quantities for these values, even if
the quantity is dimensionless. Shape, finiteness, and positivity checks use
TypedUnit; geometry, interpolation, and solver relationships stay in SuPyMode.
The generic interpolators continue to operate on numerical SI arrays.

.. code-block:: python

   from SuPyMode import ureg
   from SuPyMode.taper import AlphaProfile

   profile = AlphaProfile(initial_radius=62.5 * ureg.micrometer)

   profile.add_taper_segment(
       alpha=0.0,
       initial_heating_length=10 * ureg.millimeter,
       stretching_length=20 * ureg.millimeter,
   )

   profile.initialize()

Run a first fiber model
-----------------------

Lengths in this workflow are expressed in metres. The small grid and taper
interval below are for learning the API; check convergence for scientific use.

.. code-block:: python

   from SuPyMode.workflow import (
       Workflow, fiber_loader, Boundaries, BoundaryValue, DomainAlignment,
   )

   fiber = fiber_loader.load_fiber(
       fiber_name="SMF28", clad_refractive_index=1.4444, remove_cladding=False,
   )
   workflow = Workflow(
       fiber_list=[fiber], wavelength=1550e-9,
       resolution=40, n_step=50, n_sorted_mode=2, n_added_mode=2,
       itr_initial=1.0, itr_final=0.8,
       x_bounds=DomainAlignment.LEFT, y_bounds=DomainAlignment.BOTTOM,
       boundaries=[Boundaries(
           right=BoundaryValue.SYMMETRIC, top=BoundaryValue.SYMMETRIC,
       )],
   )
   workflow.initialize_geometry(plot=False)
   workflow.run_solver()
   workflow.plot(plot_type="index")

Inspect fields as well as effective indices when checking mode order. Continue
with :doc:`guide` or the :doc:`gallery/workflow/plot_workflow_01` example.

If installation or rendering fails, see :doc:`troubleshooting`.
