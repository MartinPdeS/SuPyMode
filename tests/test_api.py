#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from SuPyMode.workflow import configuration, Workflow, fiber_catalogue, Boundaries


def test_superset_plot():
    """
    Test the creation and manipulation of a mode superset in a SuPyMode workflow.

    This function performs the following:
    - Loads fibers from the fiber catalogue.
    - Configures a fused structure for the fiber cladding.
    - Initializes a workflow with specific parameters.
    - Accesses and manipulates various solver and mode properties.
    - Tests the labeling and resetting of supermodes.
    """

    # Load two identical fibers from the catalogue
    fibers = [
        fiber_catalogue.load_fiber('SMF28', wavelength=1550e-9),
        fiber_catalogue.load_fiber('SMF28', wavelength=1550e-9),
    ]

    # Define the cladding structure
    fused_structure = configuration.ring.FusedProfile_02x02

    # Initialize the workflow
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

    # Access the solver
    solver = workflow.solver

    # Example operation: convert eigenvalue to index
    _ = solver.eigen_value_to_index(3e6)

    # Access the solver's coordinate system
    _ = solver.coordinate_system

    # Access the first mode in the superset
    mode = workflow.superset[0]

    # Access various mode properties
    _ = mode.geometry
    _ = mode.coordinate_system
    _ = mode.itr_list
    _ = mode.model_parameters
    _ = mode.binding_number

    # Perform field interpolation for specific iterations and slices
    _ = mode.get_field_interpolation(itr=1.0)
    _ = mode.get_field_interpolation(slice_number=3)

    # Label and reset labels for the supermodes
    workflow.superset.label_supermodes('a', 'b')
    workflow.superset.reset_labels()

    workflow.superset.sort_modes('beta')

    workflow.superset.sort_modes('symmetry+beta')

    workflow.superset.export_data(filename='test_data')


if __name__ == '__main__':
    pytest.main([__file__])
