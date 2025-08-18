#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from SuPyMode.workflow import Profile, Workflow, fiber_loader, Boundaries, DomainAlignment, BoundaryValue, StructureType
import matplotlib.pyplot as plt
# Define a list of fused fiber structures from the configuration module

clad_structure_0 = Profile()

clad_structure_0.add_structure(
    structure_type=StructureType.CIRCULAR,
    number_of_fibers=2,
    fusion_degree=0.3,
    fiber_radius=62.5e-6,
    compute_fusing=True
)

clad_structure_0.refractive_index = 1.4444

clad_structure_1 = Profile()

clad_structure_1.add_structure(
    structure_type=StructureType.CIRCULAR,
    number_of_fibers=3,
    fusion_degree=0.3,
    fiber_radius=62.5e-6,
    compute_fusing=True
)

clad_structure_1.refractive_index = 1.4444

fused_structure_list = [clad_structure_0, clad_structure_1]


@pytest.mark.parametrize('fused_structure', fused_structure_list, ids=lambda x: f"n_fiber: {len(x.cores)}")
def test_workflow(fused_structure):
    """
    Test the Workflow initialization with various configurations of fused fiber structures.

    This test verifies that the Workflow can be instantiated with different numbers of fibers
    specified by the fused structures without raising any exceptions.

    Args:
        fused_structure (FusedProfile): A fused fiber structure configuration object.
    """
    # Load fibers based on the number required by the fused structure, all with the same specified wavelength
    fibers = [
        fiber_loader.load_fiber('DCF1300S_33', clad_refractive_index=1.4444)
        for _ in range(len(fused_structure.cores))
    ]

    # Instantiate the Workflow with the given configuration
    workflow = Workflow(
        fiber_list=fibers,
        clad_structure=fused_structure,
        wavelength=1550e-9,
        resolution=10,
        x_bounds=DomainAlignment.LEFT,
        y_bounds=DomainAlignment.CENTERING,
        boundaries=[Boundaries(right=BoundaryValue.SYMMETRIC)],
        n_sorted_mode=2,
        n_added_mode=2,
        debug_mode=0,
    )

    workflow.initialize_geometry()

    workflow.run_solver()

    workflow.generate_pdf_report(filename='test_0')

    plt.close()

    # workflow.save_superset_instance(filename='test_0')

    # Assert that the workflow instance has been successfully created (basic check)
    assert workflow is not None, "Workflow should be successfully instantiated with the given configurations."


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
