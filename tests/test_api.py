#!/usr/bin/env python
# -*- coding: utf-8 -*-

from SuPyMode.workflow import configuration, Workflow, fiber_catalogue, Boundaries


def test_superset_plot():
    fibers = [
        fiber_catalogue.load_fiber('SMF28', wavelength=1550e-9),
        fiber_catalogue.load_fiber('SMF28', wavelength=1550e-9),
    ]

    fused_structure = configuration.ring.FusedProfile_02x02

    workflow = Workflow(
        fiber_list=fibers,
        clad_structure=fused_structure,
        wavelength=1550e-9,
        resolution=20,
        x_bounds="left",
        y_bounds="centering",
        debug_mode=0,
        auto_label=True,
        boundaries=[Boundaries(right='symmetric')],
        n_sorted_mode=2,
        n_added_mode=2,
    )

    solver = workflow.solver

    _ = solver.eigen_value_to_index(3e6)

    _ = solver.coordinate_system

    mode = workflow.superset[0]

    _ = mode.geometry

    _ = mode.coordinate_system

    _ = mode.itr_list

    _ = mode.model_parameters

    _ = mode.binding_number

    _ = mode.get_field_interpolation(itr=1.0)

    _ = mode.get_field_interpolation(slice_number=3)

# -
