#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import numpy
from SuPyMode.workflow import Workflow, fiber_loader, Boundaries, DomainAlignment, BoundaryValue
from PyFiberModes.__future__ import get_normalized_LP_coupling
import PyFiberModes


def test_normalized_coupling(
        wavelength: float = 1550e-9,
        resolution: int = 140,
        fiber_name: str = 'test_multimode_fiber',
        n_step: int = 80,
        itr_final: float = 0.8,
        **kwargs):
    """
    Test the normalized coupling between LP01 and LP02 modes using both analytical and simulation methods.

    This test compares the analytical normalized coupling derived from PyFiberModes with the simulation output
    from SuPyMode. If significant discrepancies exist beyond the acceptable threshold, it raises an assertion error.

    Args:
        wavelength (float): Operating wavelength of the fiber.
        resolution (int): Spatial resolution of the simulation.
        n_step (int): Number of steps in the simulation.
        itr_final (float): Final iteration value.
        plot_results (bool): Whether to plot results using MPSPlots.
        x_bounds (list | str): Horizontal bounds for the simulation.
        y_bounds (list | str): Vertical bounds for the simulation.
        **kwargs: Additional keyword arguments to pass to the Workflow.
    """
    fiber = fiber_loader.load_fiber(fiber_name, clad_refractive_index=1.4444, remove_cladding=False)

    # Set up the workflow with specified parameters and boundaries
    workflow = Workflow(
        fiber_list=[fiber],
        wavelength=wavelength,
        resolution=resolution,
        x_bounds=DomainAlignment.LEFT,
        y_bounds=DomainAlignment.BOTTOM,
        boundaries=[Boundaries(right=BoundaryValue.SYMMETRIC, top=BoundaryValue.SYMMETRIC)],
        n_sorted_mode=4,
        n_added_mode=8,
        debug_mode=1,
        auto_label=True,
        itr_final=itr_final,
        n_step=n_step,
        **kwargs
    )

    superset = workflow.get_superset()

    itr_list = superset.model_parameters.itr_list

    # Load and scale the analytical fiber
    smf28 = PyFiberModes.fiber.load_fiber(fiber_name, wavelength=wavelength, add_air_layer=False)

    # Compute the analytical normalized coupling
    analytical = numpy.empty(itr_list.shape)
    for idx, itr in enumerate(itr_list):
        scaled_fiber = smf28.scale(factor=itr)
        analytical[idx] = get_normalized_LP_coupling(
            fiber=scaled_fiber,
            mode_0=PyFiberModes.LP01,
            mode_1=PyFiberModes.LP02
        )

    analytical = numpy.abs(analytical)
    # Extract the simulation data
    simulation = numpy.abs(superset.LP01.normalized_coupling.get_values(superset.LP02))

    # Calculate error metrics and validate against acceptable threshold
    error = numpy.abs(analytical - simulation)
    relative_error = error / numpy.abs(analytical)
    mean_relative_error = relative_error.mean()

    if mean_relative_error > 0.1:
        raise AssertionError(f"Discrepancy between computed and analytical normalized coupling: [Mean Error: {error.mean()}, Mean Relative Error: {mean_relative_error}]")


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
