#!/usr/bin/env python
# -*- coding: utf-8 -*-

from unittest.mock import patch
import pytest
import matplotlib.pyplot as plt
from SuPyMode.workflow import configuration, Workflow, fiber_catalogue, Boundaries


@pytest.fixture
def setup_workflow():
    """ Fixture to set up the workflow with common settings for tests. """
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


parameter_list = [
    'index', 'beta', 'eigen-value', 'normalized-coupling', 'adiabatic', 'field'
]


@pytest.mark.parametrize("plot_type", parameter_list)
@patch("matplotlib.pyplot.show")
def test_superset_plot(mock_show, setup_workflow, plot_type):
    """
    Tests plotting functionalities for various superset properties using mocked display to verify plots are invoked without display.

    Args:
        mock_show (MagicMock): Mock for matplotlib.pyplot.show to prevent GUI display during testing.
        setup_workflow (Workflow): A Workflow instance set up via a fixture to standardize test setup.
        plot_type (str): Type of plot to generate and test.
    """
    superset = setup_workflow.superset
    superset.plot(plot_type=plot_type)
    mock_show.assert_called_once()
    mock_show.reset_mock()


@pytest.mark.parametrize("plot_type", parameter_list)
@patch("matplotlib.pyplot.show")
def test_representation_plot(mock_show, setup_workflow, plot_type):
    """
    Tests individual mode plotting functionalities within the superset to ensure each plot type can be generated.

    This function verifies that plots for each mode's specific attribute can be called and displayed (mocked).

    Args:
        mock_show (MagicMock): Mock for matplotlib.pyplot.show to prevent GUI display during testing.
        setup_workflow (Workflow): A Workflow instance set up via a fixture.
        plot_type (str): Type of mode-specific plot to generate and test.
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
    pytest.main([__file__])
# -
