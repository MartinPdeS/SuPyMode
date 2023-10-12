#!/usr/bin/env python
# -*- coding: utf-8 -*-

from SuPyMode.workflow import Workflow, fiber_catalogue, configuration, Boundaries2D
from SuPyMode.profiles import AlphaProfile

fiber = [
    fiber_catalogue.DCF1300S_33(wavelength=1550e-9),
    fiber_catalogue.DCF1300S_20(wavelength=1550e-9),
    fiber_catalogue.DCF1300S_33(wavelength=1550e-9)
]


def test_propagation():
    profile = AlphaProfile(
        initial_radius=62.5e-6,
        symmetric=False,
        label='test profile'
    )

    profile.add_taper_segment(
        alpha=0,
        initial_heating_length=5e-3,
        stretching_length=2e-3,
        n_point=200
    )

    profile.initialize()

    workflow = Workflow(
        fiber_list=fiber,
        clad_structure=configuration.ring.FusedProfile_03x03,
        fusion_degree=0.3,
        wavelength=1550e-9,
        resolution=40,
        x_bounds="centering-left",
        y_bounds="centering",
        boundaries=[Boundaries2D(right='symmetric')],
        n_sorted_mode=2,
        n_added_mode=2,
        debug_mode=False
    )

    superset = workflow.get_superset()

    figure, _ = superset.plot_propagation(
        profile=profile,
        initial_amplitude=[0, 1],
        coupling='normalized',
        method='RK45',
        show_energy=True,
        show_amplitude=False,
        max_step=1550e-9 / 5000
    )

    # figure.show()


if __name__ == '__main__':
    test_propagation()

# -
