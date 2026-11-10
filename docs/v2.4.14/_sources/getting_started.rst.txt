Getting started
===============

Install and verify
------------------

.. code-block:: console

   $ python -m pip install SuPyMode
   $ python -c "import SuPyMode; print(SuPyMode.__version__)"

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
