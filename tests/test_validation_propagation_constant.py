#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries2D
from PyFiberModes.fiber import load_fiber
from PyFiberModes import LP01
from MPSPlots.render2D import SceneList


def test_propagation_constant(wavelength: float = 1550e-8, resolution: int = 120, n_step: int = 200, itr_final: float = 0.2, plot_results: bool = False, **kwargs):

    fiber = fiber_catalogue.load_fiber('SMF28', wavelength=wavelength)

    fiber = fiber.scale(factor=10)

    workflow = Workflow(
        fiber_list=[fiber],
        wavelength=wavelength,
        resolution=resolution,
        x_bounds=[-200e-6, 0],
        y_bounds=[-200e-6, 0],
        boundaries=[Boundaries2D(right='symmetric', top='symmetric')],
        n_sorted_mode=4,
        n_added_mode=8,
        auto_label=True,
        itr_final=itr_final,
        n_step=n_step,
        **kwargs
    )

    superset = workflow.get_superset()
    itr_list = superset.itr_list

    smf28 = load_fiber(
        fiber_name='SMF28',
        wavelength=wavelength,
        add_air_layer=False
    )

    smf28 = smf28.scale(factor=10)

    analytical = numpy.empty(itr_list.shape)
    for idx, itr in enumerate(itr_list):
        _fiber = smf28.scale(factor=itr)
        analytical[idx] = _fiber.get_propagation_constant(mode=LP01)

    simulation = superset.LP01.beta.get_values()

    if plot_results:
        figure = SceneList()

        ax = figure.append_ax(
            x_label='Inverse taper ratio',
            y_label='Propagation constant',
            show_legend=True
        )

        ax.add_line(x=itr_list, y=analytical, label='Analytical')
        ax.add_scatter(x=itr_list, y=simulation, label='SuPyMode')

        figure.show()

    discrepency_bool = numpy.isclose(analytical, simulation, rtol=1e-2)
    if discrepency_bool.mean() < 0.9:
        error = abs(analytical - simulation)
        relative_error = error / abs(analytical)
        raise AssertionError(f"Discrepency between the computed and analytical normalized coupling, discrepency: [mean error: {error.mean()}, mean relative error: {relative_error.mean()}]")


if __name__ == '__main__':
    test_normalized_coupling(
        wavelength=1550e-9,
        resolution=120,
        itr_final=0.6,
        n_step=100,
        plot_field=True,
        debug_mode=True,
        plot_results=True
    )

# -
