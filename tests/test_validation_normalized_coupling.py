#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries
from PyFiberModes.fiber import load_fiber
from PyFiberModes.__future__ import get_normalized_LP_coupling
from PyFiberModes import LP01, LP02
from MPSPlots.render2D import SceneList


def test_normalized_coupling(
        wavelength: float = 1550e-9,
        resolution: int = 200,
        fiber_name: str = 'SMF28',
        n_step: int = 400,
        itr_final: float = 0.6,
        plot_results: bool = False,
        x_bounds: list | str = [-200e-6, 0],
        y_bounds: list | str = [-200e-6, 0],
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
    # Load and scale the fiber
    fiber = fiber_catalogue.load_fiber(fiber_name, wavelength=wavelength)
    fiber = fiber.scale(factor=10)

    # Set up the workflow with specified parameters and boundaries
    workflow = Workflow(
        fiber_list=[fiber],
        wavelength=wavelength,
        resolution=resolution,
        x_bounds=x_bounds,
        y_bounds=y_bounds,
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

    # Load and scale the analytical fiber
    smf28 = load_fiber(fiber_name, wavelength=wavelength, add_air_layer=False)
    smf28 = smf28.scale(factor=10)

    # Compute the analytical normalized coupling
    analytical = numpy.empty(itr_list.shape)
    for idx, itr in enumerate(itr_list):
        scaled_fiber = smf28.scale(factor=itr)
        analytical[idx] = get_normalized_LP_coupling(
            fiber=scaled_fiber,
            mode_0=LP01,
            mode_1=LP02
        )

    # Extract the simulation data
    simulation = -numpy.abs(superset.LP01.normalized_coupling.get_values(superset.LP02))

    # Optionally plot the results for visual comparison
    if plot_results:
        figure = SceneList()
        ax = figure.append_ax(x_label='Inverse taper ratio', y_label='Normalized coupling', show_legend=True)
        ax.add_line(x=itr_list, y=analytical, label='Analytical')
        ax.add_scatter(x=itr_list, y=simulation, label='SuPyMode Simulation')
        figure.show()

    # Calculate error metrics and validate against acceptable threshold
    error = numpy.abs(analytical - simulation)
    relative_error = error / numpy.abs(analytical)
    mean_relative_error = relative_error.mean()

    if mean_relative_error > 0.1:
        raise AssertionError(f"Discrepancy between computed and analytical normalized coupling: "
                             f"[Mean Error: {error.mean()}, Mean Relative Error: {mean_relative_error}]")
