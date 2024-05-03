#!/usr/bin/env python
# -*- coding: utf-8 -*-

from SuPyMode.workflow import Workflow, fiber_catalogue, configuration, Boundaries
from SuPyMode.profiles import AlphaProfile

# Predefine fibers to be used across tests
fibers = [
    fiber_catalogue.load_fiber('DCF1300S_33', wavelength=1550e-9),
    fiber_catalogue.load_fiber('DCF1300S_20', wavelength=1550e-9),
    fiber_catalogue.load_fiber('DCF1300S_33', wavelength=1550e-9)
]


def test_propagation():
    """
    Test the propagation of modes through a tapered fiber setup using the AlphaProfile and Workflow classes.

    This function verifies that the system can initialize and execute propagation calculations without errors
    and that the propagation method correctly handles the specified inputs under a basic scenario.

    The test involves setting up a fiber taper segment and using the Workflow to propagate modes through this configuration.
    """
    # Set up the taper profile with specific characteristics
    profile = AlphaProfile(
        initial_radius=62.5e-6,
        symmetric=False,
        label='test profile'
    )

    profile.add_taper_segment(
        alpha=0,  # No taper angle, implying a uniform taper
        initial_heating_length=1e-3,  # Start length of the taper
        stretching_length=5e-3,  # Total length of the taper
        n_point=200  # Number of points to calculate within the taper segment
    )

    # Initialize the workflow with the specified configuration and fibers
    workflow = Workflow(
        fiber_list=fibers,
        clad_structure=configuration.ring.FusedProfile_03x03,
        fusion_degree='auto',
        wavelength=1550e-9,
        resolution=40,
        x_bounds="left",
        y_bounds="centering",
        boundaries=[Boundaries(right='symmetric')],
        n_sorted_mode=2,
        n_added_mode=2,
        debug_mode=0
    )

    # Get the superset generated by the workflow, which contains the propagation environment
    superset = workflow.get_superset()

    # Perform the propagation calculation using the specified profile and parameters
    propagation_results = superset.propagate(
        profile=profile,
        initial_amplitude=[0, 1],  # Initial amplitudes for the modes
        add_coupling=True,  # Consider mode coupling in the calculation
        method='RK45',  # Integration method to use
        max_step=1550e-9 / 400  # Maximum step size for the integration
    )

    # Asserts to verify that propagation results are generated
    assert propagation_results is not None, "Propagation should successfully return results"
