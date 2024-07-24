import pytest
from unittest.mock import patch
from SuPyMode.profiles import AlphaProfile


@pytest.fixture
def alpha_profile():
    """Fixture to create an AlphaProfile instance."""
    profile = AlphaProfile(initial_radius=1)
    return profile


@patch("matplotlib.pyplot.show")
def test_build_single_segment_profile(mock_show, alpha_profile):
    """
    Test building a profile with a single taper segment.
    """
    alpha_profile.add_taper_segment(
        alpha=0,
        initial_heating_length=10e-3,
        stretching_length=0.2e-3 * 200
    )
    alpha_profile.initialize()
    alpha_profile.plot()

    mock_show.assert_called_once()


@patch("matplotlib.pyplot.show")
def test_build_two_segment_profile(mock_show, alpha_profile):
    """
    Test building a profile with two taper segments.
    """
    alpha_profile.add_taper_segment(
        alpha=0,
        initial_heating_length=10e-3,
        stretching_length=0.2e-3 * 200
    )
    alpha_profile.add_taper_segment(
        alpha=0,
        initial_heating_length=3e-3,
        stretching_length=0.2e-3 * 200
    )
    alpha_profile.initialize()
    alpha_profile.plot()

    mock_show.assert_called_once()


@patch("matplotlib.pyplot.show")
def test_build_asymmetric_profile(mock_show):
    """
    Test building an asymmetric profile with a single taper segment.
    """
    asymmetric_profile = AlphaProfile(initial_radius=1, symmetric=False)
    asymmetric_profile.add_taper_segment(
        alpha=0,
        initial_heating_length=10e-3,
        stretching_length=0.2e-3 * 200
    )
    asymmetric_profile.initialize()
    asymmetric_profile.plot()

    mock_show.assert_called_once()


def test_generate_propagation_gif(alpha_profile):
    """
    Test generating a GIF from the profile data.
    """
    alpha_profile.add_taper_segment(
        alpha=0,
        initial_heating_length=3e-3,
        stretching_length=0.2e-3 * 200
    )
    alpha_profile.initialize()
    alpha_profile.generate_propagation_gif(number_of_frames=10)

    # Assertions can be added here if generate_propagation_gif outputs testable results
