#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from SuPyMode.workflow import configuration, Workflow, fiber_catalogue, Boundaries2D

fused_structure_list = [
    configuration.ring.FusedProfile_01x01,
    configuration.ring.FusedProfile_02x02,
    configuration.ring.FusedProfile_03x03,
    configuration.ring.FusedProfile_04x04,
    configuration.ring.FusedProfile_07x07
]


@pytest.mark.parametrize('fused_structure', fused_structure_list, ids=fused_structure_list)
def test_fused_structure_workflow(fused_structure, debug_mode: bool = False):
    fibers = [
        fiber_catalogue.DCF1300S_33(wavelength=1550e-9) for _ in range(fused_structure.number_of_fibers)
    ]

    _ = Workflow(
        fiber_list=fibers,
        clad_structure=fused_structure,
        wavelength=1550e-9,
        resolution=30,
        x_bounds="centering-left",
        y_bounds="centering",
        boundaries=[Boundaries2D(right='symmetric')],
        n_sorted_mode=2,
        n_added_mode=2,
        debug_mode=debug_mode
    )


if __name__ == '__main__':
    for fused_structure in fused_structure_list:
        test_fused_structure_workflow(fused_structure=fused_structure, debug_mode=True)
# -