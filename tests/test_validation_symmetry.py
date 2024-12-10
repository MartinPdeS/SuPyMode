#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import numpy
from SuPyMode.workflow import Workflow, configuration, fiber_catalogue, Boundaries

BOUNDARIES_LIST = [
    dict(x_bounds='left', boundaries=[Boundaries(right='symmetric')]),
    dict(x_bounds='right', boundaries=[Boundaries(left='symmetric')]),
    dict(y_bounds='top', boundaries=[Boundaries(bottom='symmetric')]),
    dict(y_bounds='bottom', boundaries=[Boundaries(top='symmetric')]),
]


@pytest.mark.parametrize('boundaries', BOUNDARIES_LIST)
def test_symmetry(boundaries, fiber_name: str = 'DCF1300S_33', wavelength: float = 1.55e-6, resolution: int = 21):
    """
    Tests the effect of symmetric and asymmetric boundary conditions on the computed indices in a fiber modeling workflow.

    This function sets up two workflows with identical settings except for the boundary conditions: one symmetric and one asymmetric.
    It then compares the index data produced by these workflows to ensure that symmetry conditions are implemented correctly,
    and that they have the expected impact on the results.

    Raises:
        ValueError: If the mean discrepancy between symmetric and asymmetric configurations exceeds the tolerance threshold.
    """
    kwargs = dict(
        wavelength=wavelength,
        resolution=resolution,
        n_step=10,
        itr_final=0.9,
        debug_mode=2,
        n_sorted_mode=1,
        n_added_mode=2,
    )

    fiber_list = [fiber_catalogue.load_fiber(fiber_name, wavelength=wavelength)]
    clad_structure = configuration.ring.FusedProfile_01x01

    # Setup for asymmetric boundary conditions
    reference_workflow = Workflow(
        fiber_list=fiber_list,
        clad_structure=clad_structure,
        **kwargs,
        x_bounds="centering",
        y_bounds="centering",
        boundaries=[Boundaries()],
    )

    # Setup for symmetric boundary conditions
    left_workflow = Workflow(
        fiber_list=fiber_list,
        clad_structure=clad_structure,
        **boundaries,
        **kwargs
    )

    # Compare the effective index data between the two configurations
    discrepancy = numpy.isclose(
        reference_workflow.superset[0].index.data,
        left_workflow.superset[0].index.data,
        atol=1e-10,
        rtol=1e-10
    )

    difference = abs(reference_workflow.superset[0].index.data - left_workflow.superset[0].index.data)

    if numpy.mean(discrepancy) <= 0.9:
        raise ValueError(f"Mismatch [{numpy.mean(difference):.5e}] between: non-symmetric and symmetric symmetry-based formulation of the numerical problem.")


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
