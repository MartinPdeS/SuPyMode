#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries
from PyFiberModes.fiber import load_fiber
from PyFiberModes import LP01
from MPSPlots.render2D import SceneList


def test_propagation_constant(
        wavelength: float = 1550e-9,
        resolution: int = 120,
        n_step: int = 200,
        itr_final: float = 0.2,
        plot_results: bool = False, **kwargs):
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
    fiber = fiber_catalogue.load_fiber('SMF28', wavelength=wavelength)
    fiber = fiber.scale(factor=10)

    # Initialize the simulation workflow
    workflow = Workflow(
        fiber_list=[fiber],
        wavelength=wavelength,
        resolution=resolution,
        x_bounds=[-200e-6, 0],
        y_bounds=[-200e-6, 0],
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
    itr_list = superset.itr_list

    # Load and scale the analytical fiber model
    smf28 = load_fiber(
        fiber_name='SMF28',
        wavelength=wavelength,
        add_air_layer=False
    )
    smf28 = smf28.scale(factor=10)

    # Compute the analytical propagation constant for each scaling iteration
    analytical = numpy.array([smf28.scale(factor=itr).get_propagation_constant(mode=LP01) for itr in itr_list])

    # Retrieve simulation results
    simulation = superset.LP01.beta.get_values()

    # Optionally plot the results
    if plot_results:
        figure = SceneList()
        ax = figure.append_ax(x_label='Inverse taper ratio', y_label='Propagation constant', show_legend=True)
        ax.add_line(x=itr_list, y=analytical, label='Analytical')
        ax.add_scatter(x=itr_list, y=simulation, label='Simulation')
        figure.show()

    # Check discrepancies
    discrepancies = numpy.isclose(analytical, simulation, rtol=1e-2)
    if numpy.mean(discrepancies) < 0.9:
        error = numpy.abs(analytical - simulation)
        relative_error = error / numpy.abs(analytical)
        raise AssertionError(f"Discrepancy between computed and analytical propagation constants. "
                             f"Mean Error: {error.mean()}, Mean Relative Error: {relative_error.mean()}")

