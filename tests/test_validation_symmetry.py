#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from SuPyMode.workflow import Workflow, configuration, fiber_catalogue, Boundaries


def test_symmetry(fiber_name: str = 'DCF1300S_33', wavelength: float = 1.55e-6, resolution: int = 160):
    """
    Tests the effect of symmetric and asymmetric boundary conditions on the computed indices in a fiber modeling workflow.

    This function sets up two workflows with identical settings except for the boundary conditions: one symmetric and one asymmetric.
    It then compares the index data produced by these workflows to ensure that symmetry conditions are implemented correctly,
    and that they have the expected impact on the results.

    Raises:
        ValueError: If the mean discrepancy between symmetric and asymmetric configurations exceeds the tolerance threshold.
    """
    # Setup for symmetric boundary conditions
    symmetric_workflow = Workflow(
        fiber_list=[fiber_catalogue.load_fiber(fiber_name, wavelength=wavelength)],
        clad_structure=configuration.ring.FusedProfile_01x01,
        wavelength=wavelength,
        resolution=resolution,
        n_step=10,
        itr_final=0.9,
        x_bounds="left",  # Assuming symmetrically influenced by the 'left' boundary
        y_bounds="centering",
        boundaries=[Boundaries(right='symmetric')],
        n_sorted_mode=1,
        n_added_mode=2,
        debug_mode=0,
    )

    # Setup for asymmetric boundary conditions
    asymmetric_workflow = Workflow(
        fiber_list=[fiber_catalogue.load_fiber(fiber_name, wavelength=wavelength)],
        clad_structure=configuration.ring.FusedProfile_01x01,
        wavelength=wavelength,
        resolution=resolution,
        n_step=10,
        itr_final=0.9,
        x_bounds="centering",  # Centered, implying no special treatment to symmetry
        y_bounds="centering",
        boundaries=[Boundaries()],  # No special symmetry boundaries
        n_sorted_mode=1,
        n_added_mode=2,
        debug_mode=0,
    )

    # Compare the effective index data between the two configurations
    symmetric_index_data = symmetric_workflow.superset[0].index.data
    asymmetric_index_data = asymmetric_workflow.superset[0].index.data
    discrepancy = numpy.isclose(
        symmetric_index_data,
        asymmetric_index_data,
        atol=1e-5,
        rtol=1e-5
    )

    if numpy.mean(discrepancy) <= 0.9:
        raise ValueError(f"Symmetric and non-symmetric index data do not match sufficiently; mean error: {numpy.mean(discrepancy)}")
# -
