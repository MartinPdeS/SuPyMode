#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from SuPyMode.workflow import configuration, Workflow, fiber_catalogue, Boundaries

# Define a list of fused fiber structures from the configuration module
fused_structure_list = [
    configuration.ring.FusedProfile_01x01,
    configuration.ring.FusedProfile_02x02,
    configuration.ring.FusedProfile_03x03,
    configuration.ring.FusedProfile_04x04,
    configuration.ring.FusedProfile_07x07
]

# Test IDs for better readability in pytest output
test_ids = [f"Profile_{structure.number_of_fibers}x{structure.number_of_fibers}" for structure in fused_structure_list]


@pytest.mark.parametrize('fused_structure', fused_structure_list, ids=test_ids)
def test_fused_structure_workflow(fused_structure):
    """
    Test the Workflow initialization with various configurations of fused fiber structures.

    This test verifies that the Workflow can be instantiated with different numbers of fibers
    specified by the fused structures without raising any exceptions.

    Args:
        fused_structure (FusedProfile): A fused fiber structure configuration object.
    """
    # Load fibers based on the number required by the fused structure, all with the same specified wavelength
    fibers = [
        fiber_catalogue.load_fiber('DCF1300S_33', wavelength=1550e-9)
        for _ in range(fused_structure.number_of_fibers)
    ]

    # Instantiate the Workflow with the given configuration
    workflow = Workflow(
        fiber_list=fibers,
        clad_structure=fused_structure,
        wavelength=1550e-9,
        resolution=30,
        x_bounds="left",
        y_bounds="centering",
        boundaries=[Boundaries(right='symmetric')],
        n_sorted_mode=2,
        n_added_mode=2,
        debug_mode=0,
        fusion_degree='auto'
    )

    # Assert that the workflow instance has been successfully created (basic check)
    assert workflow is not None, "Workflow should be successfully instantiated with the given configurations."
