import pytest
from unittest.mock import patch
from SuPyMode.profiles import AlphaProfile


@pytest.fixture(scope="module")
def alpha_profile():
    """
    Fixture to create an AlphaProfile instance with initial parameters.

    This fixture is shared across multiple tests to avoid redundant reinitialization,
    improving test performance.

    Returns
    -------
    AlphaProfile
        An initialized AlphaProfile instance.
    """
    return AlphaProfile(initial_radius=1)


@pytest.fixture(scope="module")
def asymmetric_alpha_profile():
    """
    Fixture to create an asymmetric AlphaProfile instance.

    Returns
    -------
    AlphaProfile
        An initialized asymmetric AlphaProfile instance.
    """
    return AlphaProfile(initial_radius=1, symmetric=False)


@patch("matplotlib.pyplot.show")
def test_build_single_segment_profile(mock_show, alpha_profile):
    """
    Test creating and plotting a profile with a single taper segment.

    Parameters
    ----------
    mock_show : MagicMock
        Mock for `matplotlib.pyplot.show` to prevent actual plot display.
    alpha_profile : AlphaProfile
        The profile fixture to use in the test.
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
    Test creating and plotting a profile with two taper segments.

    Parameters
    ----------
    mock_show : MagicMock
        Mock for `matplotlib.pyplot.show` to prevent actual plot display.
    alpha_profile : AlphaProfile
        The profile fixture to use in the test.
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
def test_build_asymmetric_profile(mock_show, asymmetric_alpha_profile):
    """
    Test creating and plotting an asymmetric profile with a single taper segment.

    Parameters
    ----------
    mock_show : MagicMock
        Mock for `matplotlib.pyplot.show` to prevent actual plot display.
    asymmetric_alpha_profile : AlphaProfile
        The asymmetric profile fixture to use in the test.
    """
    asymmetric_alpha_profile.add_taper_segment(
        alpha=0,
        initial_heating_length=10e-3,
        stretching_length=0.2e-3 * 200
    )
    asymmetric_alpha_profile.initialize()
    asymmetric_alpha_profile.plot()

    mock_show.assert_called_once()


def test_generate_propagation_gif(alpha_profile):
    """
    Test generating a propagation GIF from the profile data.

    Parameters
    ----------
    alpha_profile : AlphaProfile
        The profile fixture to use in the test.
    """
    alpha_profile.add_taper_segment(
        alpha=0,
        initial_heating_length=3e-3,
        stretching_length=0.2e-3 * 200
    )
    alpha_profile.initialize()
    alpha_profile.generate_propagation_gif(number_of_frames=10)

    # Add assertions if generate_propagation_gif outputs any testable results


if __name__ == "__main__":
    pytest.main([__file__])
