#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries2D
from PyFiberModes.fiber import load_fiber
from PyFiberModes.__future__ import get_normalized_LP_coupling
from PyFiberModes import LP01, LP02
from MPSPlots.render2D import SceneList


def test_normalized_coupling(
        wavelength: float = 1550e-9,
        resolution: int = 160,
        n_step: int = 200,
        itr_final: float = 0.6,
        plot_results: bool = False,
        **kwargs):

    fiber = fiber_catalogue.SMF28(wavelength=wavelength)

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
        analytical[idx] = get_normalized_LP_coupling(fiber=_fiber, mode_0=LP01, mode_1=LP02)

    simulation = -abs(superset.LP01.normalized_coupling.get_values(superset.LP02))

    if plot_results:
        figure = SceneList()

        ax = figure.append_ax(
            x_label='Inverse taper ratio',
            y_label='Normalized coupling',
            show_legend=True
        )

        ax.add_line(x=itr_list, y=analytical, label='Analytical')
        ax.add_scatter(x=itr_list, y=simulation, label='SuPyMode')

        figure.show()

    error = abs(analytical - simulation)
    relative_error = error / abs(analytical)
    if relative_error.mean() > 0.1:
        raise AssertionError(f"Discrepency between the computed and analytical normalized coupling, discrepency: [mean error: {error.mean()}, mean relative error: {relative_error.mean()}]")


if __name__ == '__main__':
    test_normalized_coupling(
        debug_mode=True,
        plot_results=True
    )

# -
