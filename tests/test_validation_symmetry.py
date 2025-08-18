#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest
import numpy
from SuPyMode.workflow import Workflow, fiber_loader, Boundaries, Profile, DomainAlignment, BoundaryValue, StructureType

BOUNDARIES_LIST = [
    dict(x_bounds=DomainAlignment.LEFT, boundaries=[Boundaries(right=BoundaryValue.SYMMETRIC)]),
    dict(x_bounds=DomainAlignment.RIGHT, boundaries=[Boundaries(left=BoundaryValue.SYMMETRIC)]),
    dict(y_bounds=DomainAlignment.TOP, boundaries=[Boundaries(bottom=BoundaryValue.SYMMETRIC)]),
    dict(y_bounds=DomainAlignment.BOTTOM, boundaries=[Boundaries(top=BoundaryValue.SYMMETRIC)]),
]


@pytest.fixture
def reference_clad():
    clad_structure = Profile()
    clad_structure.add_structure(
        structure_type=StructureType.CIRCULAR,
        number_of_fibers=4,
        fusion_degree=0.3,
        fiber_radius=62.5e-6,
        compute_fusing=True
    )
    clad_structure.refractive_index = 1.4444
    return clad_structure

@pytest.fixture
def symmetric_clad():
    clad_structure = Profile()
    clad_structure.add_structure(
        structure_type=StructureType.CIRCULAR,
        number_of_fibers=4,
        fusion_degree=0.3,
        fiber_radius=62.5e-6,
        compute_fusing=True
    )
    clad_structure.refractive_index = 1.4444
    return clad_structure

@pytest.mark.parametrize('boundaries', BOUNDARIES_LIST)
def test_symmetry(boundaries, reference_clad, symmetric_clad, fiber_name: str = 'DCF1300S_33', wavelength: float = 1.55e-6, resolution: int = 61):
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

    fiber_list = [
        fiber_loader.load_fiber(fiber_name, clad_refractive_index=1.4444, position=core) for core in reference_clad.cores
    ]

    # Setup for asymmetric boundary conditions
    reference_workflow = Workflow(
        fiber_list=fiber_list,
        clad_structure=reference_clad,
        **kwargs,
        x_bounds=DomainAlignment.CENTERING,
        y_bounds=DomainAlignment.CENTERING,
        boundaries=[Boundaries()],
    )

    reference_workflow.initialize_geometry()
    reference_workflow.run_solver()

    # Setup for symmetric boundary conditions
    left_workflow = Workflow(
        fiber_list=fiber_list,
        clad_structure=symmetric_clad,
        **boundaries,
        **kwargs
    )

    left_workflow.initialize_geometry()
    left_workflow.run_solver()

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
