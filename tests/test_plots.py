#!/usr/bin/env python
# -*- coding: utf-8 -*-

from unittest.mock import patch
import pytest
import matplotlib.pyplot as plt
from SuPyMode.workflow import configuration, Workflow, fiber_catalogue, Boundaries


PARAMETER_LIST = [
    'index', 'beta', 'eigen-value', 'normalized-coupling', 'adiabatic', 'field'
]


@pytest.fixture(scope="module")
def setup_workflow():
    """
    Set up a common workflow instance for testing purposes.

    This fixture initializes a Workflow instance using standard fibers and a specific fused structure.
    The workflow is set up with common parameters such as resolution, wavelength, boundaries, etc.
    It is used across multiple tests to avoid recomputation and ensure consistent test conditions.

    Returns
    -------
    Workflow
        An instance of the Workflow class configured with predefined parameters for testing.
    """
    fibers = [fiber_catalogue.load_fiber('SMF28', wavelength=1550e-9) for _ in range(2)]
    fused_structure = configuration.ring.FusedProfile_02x02

    return Workflow(
        fiber_list=fibers,
        clad_structure=fused_structure,
        wavelength=1550e-9,
        resolution=30,
        x_bounds="left",
        y_bounds="centering",
        debug_mode=0,
        boundaries=[Boundaries(right='symmetric')],
        n_sorted_mode=2,
        n_added_mode=2,
    )


@pytest.mark.parametrize("plot_type", PARAMETER_LIST)
@patch("matplotlib.pyplot.show")
def test_superset_plot(mock_show, setup_workflow, plot_type):
    """
    Test plotting functionalities for various superset properties.

    Uses a mocked display to verify that plots for different properties of the superset can be
    invoked without actually displaying the GUI, ensuring that the plot functions are called correctly.

    Parameters
    ----------
    mock_show : MagicMock
        Mock for `matplotlib.pyplot.show` to prevent GUI display during testing.
    setup_workflow : Workflow
        A Workflow instance set up via a fixture to standardize test setup.
    plot_type : str
        Type of plot to generate and test (e.g., 'index', 'beta').
    """
    superset = setup_workflow.superset
    superset.plot(plot_type=plot_type, mode_of_interest='fundamental')
    mock_show.assert_called_once()
    mock_show.reset_mock()
    plt.close()


@pytest.mark.parametrize("plot_type", PARAMETER_LIST)
@patch("matplotlib.pyplot.show")
def test_representation_plot(mock_show, setup_workflow, plot_type):
    """
    Test plotting functionalities for individual modes within the superset.

    This function ensures that each mode-specific attribute plot can be called and displayed.
    For properties like 'normalized-coupling' and 'adiabatic', which require comparing two modes,
    the test uses a second mode for completeness.

    Parameters
    ----------
    mock_show : MagicMock
        Mock for `matplotlib.pyplot.show` to prevent GUI display during testing.
    setup_workflow : Workflow
        A Workflow instance set up via a fixture.
    plot_type : str
        Type of mode-specific plot to generate and test (e.g., 'index', 'beta', 'normalized-coupling').
    """
    mode = setup_workflow.superset[0]  # Use the first mode from the superset

    # Handling for 'normalized-coupling' and 'adiabatic' which require another mode as a parameter
    if plot_type in ['normalized-coupling', 'adiabatic']:
        other_mode = setup_workflow.superset[1]  # Assuming at least two modes for these tests
        mode.plot(plot_type=plot_type, other_supermode=other_mode)
    else:
        mode.plot(plot_type=plot_type)

    mock_show.assert_called_once()
    mock_show.reset_mock()
    plt.close()


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
