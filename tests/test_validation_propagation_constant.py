#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries
import PyFiberModes


def test_propagation_constant(
        wavelength: float = 1550e-9,
        fiber_name: str = 'test_multimode_fiber',
        resolution: int = 140,
        n_step: int = 80,
        itr_final: float = 0.5,
        **kwargs):
    """
    Tests the consistency between analytical and simulated propagation constants over a range of fiber scaling iterations.

    This test scales a single fiber model and computes the propagation constant for the LP01 mode both analytically and via simulation,
    then compares these values to ensure they fall within a specified tolerance.

    Args:
        wavelength (float): Operating wavelength of the fiber (in meters).
        resolution (int): Spatial resolution of the simulation.
        n_step (int): Number of steps in the simulation.
        itr_final (float): Final iteration value for scaling the fiber.
        plot_results (bool): Flag to determine if results should be plotted.
        **kwargs: Additional keyword arguments to pass to the Workflow.
    """
    # Load and scale the reference fiber
    fiber = fiber_catalogue.load_fiber(fiber_name, wavelength=wavelength)

    # Initialize the simulation workflow
    workflow = Workflow(
        fiber_list=[fiber],
        wavelength=wavelength,
        resolution=resolution,
        x_bounds='left',
        y_bounds='bottom',
        boundaries=[Boundaries(right='symmetric', top='symmetric')],
        n_sorted_mode=4,
        n_added_mode=8,
        debug_mode=0,
        auto_label=True,
        itr_final=itr_final,
        n_step=n_step,
        **kwargs
    )

    superset = workflow.get_superset()
    itr_list = superset.model_parameters.itr_list

    # Load and scale the analytical fiber model
    pfm_fiber = PyFiberModes.fiber.load_fiber(
        fiber_name=fiber_name,
        wavelength=wavelength,
        add_air_layer=False
    )

    # Compute the analytical propagation constant for each scaling iteration
    analytical = numpy.array([pfm_fiber.scale(factor=itr).get_propagation_constant(mode=PyFiberModes.LP01) for itr in itr_list])

    # Retrieve simulation results
    simulation = superset.LP01.beta.data

    # Check discrepancies
    discrepancies = numpy.isclose(analytical, simulation, rtol=1e-2)
    if numpy.mean(discrepancies) < 0.9:
        error = numpy.abs(analytical - simulation)
        relative_error = error / numpy.abs(analytical)
        raise AssertionError(f"Discrepancy between computed and analytical propagation constants. Mean Error: {error.mean()}, Mean Relative Error: {relative_error.mean()}")


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
