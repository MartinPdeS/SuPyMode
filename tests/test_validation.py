#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy
from SuPyMode.tools.analytics.data_visualizer import DataVisualizer
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries2D, configuration
from SuPyMode.superset import SuperSet


def get_superset_from_preset(
        fiber_type,
        wavelength: float = 1550e-9,
        resolution: int = 160,
        n_step: int = 100) -> SuperSet:
    """
    Gets the superset from preset.

    :param      fiber_type:  The fiber type
    :type       fiber_type:  { type_description }
    :param      wavelength:  The wavelength
    :type       wavelength:  float
    :param      resolution:  The resolution
    :type       resolution:  int 
    :param      n_step:      The n step
    :type       n_step:      int

    :returns:   The superset from preset.
    :rtype:     Superset
    """
    fiber = fiber_type(wavelength=wavelength)

    clad_structure = configuration.ring.FusedProfile_01x01

    boundaries = [
        Boundaries2D(right='symmetric', bottom='symmetric'),
        Boundaries2D(right='symmetric', bottom='anti-symmetric')
    ]

    workflow = Workflow(
        fiber_list=[fiber],
        clad_structure=clad_structure,
        wavelength=wavelength,
        resolution=resolution,
        n_step=n_step,
        x_bounds="centering-left",
        y_bounds="centering-top",
        boundaries=boundaries,
        n_sorted_mode=6,
        n_added_mode=4,
        debug_mode=True,
        auto_label=True,
    )

    return workflow.get_superset()


def test_propagation_constant(
        fiber_type: object = fiber_catalogue.SMF28,
        mode_numbers: list = ['LP01', 'LP11_b', 'LP02'],
        wavelength: float = 1550e-9,
        resolution: int = 160,
        n_step: int = 100) -> None:
    """
    Validate the propagation constant computation with analytical prediction
    for a specific fiber type.

    :param      fiber_type:  The fiber type
    :type       fiber_type:  { type_description }
    :param      wavelength:  The wavelength
    :type       wavelength:  float
    :param      resolution:  The resolution
    :type       resolution:  int
    :param      n_step:      The n step
    :type       n_step:      int

    :returns:   The superset from preset.
    :rtype:     Superset
    """
    superset = get_superset_from_preset(
        fiber_type=fiber_type,
        wavelength=wavelength,
        resolution=resolution,
        n_step=n_step
    )

    fibermode_solver = DataVisualizer(wavelength=wavelength)

    fibermodes_data_sets = fibermode_solver.get_beta_vs_itr(
        mode_numbers=[m[:4] for m in mode_numbers],
        itr_list=superset.itr_list,
        debug_mode=False
    )

    for mode_number, data_set in zip(mode_numbers, fibermodes_data_sets):

        analytic_data = data_set.y
        simulation_data = getattr(superset, mode_number).beta.get_values()

        not_nan_idx = numpy.where(~numpy.isnan(analytic_data))

        analytic_data = analytic_data[not_nan_idx]
        simulation_data = simulation_data[not_nan_idx]

        discrepency = numpy.isclose(
            analytic_data,
            simulation_data,
            atol=1e-4,
            rtol=1e-4
        )

        if not discrepency.mean() > 0.9:
            raise ValueError(f'Discrepencies detected from analytical and simulated predictions [mean error: {discrepency.mean()}]')


if __name__ == '__main__':
    test_propagation_constant(
        fiber_type=fiber_catalogue.SMF28,
        mode_numbers=['LP01', 'LP11_b', 'LP02'],
        wavelength=1550e-9,
        resolution=20,
        n_step=10
    )

# -
