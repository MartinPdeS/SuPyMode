#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from SuPyMode.workflow import Workflow, fiber_loader, Boundaries, BoundaryValue, Profile, StructureType, DomainAlignment


@pytest.fixture(scope="module")
def precomputed_workflow():
    """
    Fixture to initialize the SuPyMode workflow for reuse across multiple tests.
    """
    fibers = [
        fiber_loader.load_fiber('SMF28', clad_refractive_index=1.4444),
        fiber_loader.load_fiber('SMF28', clad_refractive_index=1.4444),
    ]
    fused_structure = Profile()

    fused_structure.add_structure(
        structure_type=StructureType.CIRCULAR,
        number_of_fibers=2,
        fusion_degree=0.3,
        fiber_radius=62.5e-6,
        compute_fusing=True
    )

    fused_structure.refractive_index = 1.4444

    workflow = Workflow(
        fiber_list=fibers,
        clad_structure=fused_structure,
        wavelength=1550e-9,
        resolution=20,
        x_bounds=DomainAlignment.LEFT,
        y_bounds=DomainAlignment.CENTERING,
        debug_mode=0,
        auto_label=True,
        boundaries=[Boundaries(right=BoundaryValue.SYMMETRIC)],
        n_sorted_mode=2,
        n_added_mode=2,
    )

    workflow.initialize_geometry()
    workflow.run_solver()

    return workflow


def test_load_fibers(precomputed_workflow):
    """
    Test loading fibers from the fiber catalogue.
    """
    fibers = precomputed_workflow.fiber_list
    assert fibers is not None, "Failed to load fibers from the catalogue."
    assert len(fibers) == 2, "Incorrect number of fibers loaded."


def test_initialize_workflow(precomputed_workflow):
    """
    Test initializing a SuPyMode workflow with specific parameters.
    """
    assert precomputed_workflow is not None, "Workflow initialization failed."


def test_solver_properties(precomputed_workflow):
    """
    Test accessing solver properties and performing eigenvalue conversion.
    """
    solver = precomputed_workflow.solver
    assert solver.eigen_value_to_index(3e6) is not None, "Eigenvalue conversion failed."
    assert solver.coordinate_system is not None, "Solver coordinate system not accessible."


def test_mode_properties(precomputed_workflow):
    """
    Test accessing and validating properties of the first mode in the superset.
    """
    mode = precomputed_workflow.superset[0]
    assert mode.geometry is not None, "Mode geometry not accessible."
    assert mode.coordinate_system is not None, "Mode coordinate system not accessible."
    assert mode.itr_list is not None, "Mode ITR list not accessible."
    assert mode.model_parameters is not None, "Mode model parameters not accessible."
    assert mode.binding_number is not None, "Mode binding number not accessible."


def test_field_interpolation(precomputed_workflow):
    """
    Test field interpolation for a mode using ITR and slice number.
    """
    mode = precomputed_workflow.superset[0]
    assert mode.get_field_interpolation(itr=1.0) is not None, "Field interpolation by ITR failed."
    assert mode.get_field_interpolation(slice_number=3) is not None, "Field interpolation by slice number failed."


def test_superset_operations(precomputed_workflow):
    """
    Test labeling, resetting labels, sorting, and exporting supermodes.
    """
    precomputed_workflow.superset.label_supermodes('a', 'b')
    precomputed_workflow.superset.reset_labels()
    precomputed_workflow.superset.sort_modes('beta')
    precomputed_workflow.superset.sort_modes('symmetry+beta')
    precomputed_workflow.superset.export_data(filename='test_data')


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
