import pytest
from SuPyMode.workflow import AlphaProfile
from SuPyMode.plotter import generate_propagation_gif


def test_build_single_segment_profile():
    """
    Test creating and plotting a profile with a single taper segment.

    Parameters
    ----------
    alpha_profile : AlphaProfile
        The profile fixture to use in the test.
    """
    alpha_profile = AlphaProfile(initial_radius=1)

    alpha_profile.add_taper_segment(
        alpha=0, initial_heating_length=10e-3, stretching_length=0.2e-3 * 200
    )
    alpha_profile.initialize()


def test_build_two_segment_profile():
    """
    Test creating and plotting a profile with two taper segments.

    Parameters
    ----------
    alpha_profile : AlphaProfile
        The profile fixture to use in the test.
    """
    alpha_profile = AlphaProfile(initial_radius=1)

    alpha_profile.add_taper_segment(
        alpha=0, initial_heating_length=10e-3, stretching_length=0.2e-3 * 200
    )
    alpha_profile.add_taper_segment(
        alpha=0, initial_heating_length=3e-3, stretching_length=0.2e-3 * 200
    )
    alpha_profile.initialize()


def test_build_asymmetric_profile():
    """
    Test creating and plotting an asymmetric profile with a single taper segment.

    Parameters
    ----------
    asymmetric_alpha_profile : AlphaProfile
        The asymmetric profile fixture to use in the test.
    """
    asymmetric_alpha_profile = AlphaProfile(initial_radius=1, symmetric=False)

    asymmetric_alpha_profile.add_taper_segment(
        alpha=0, initial_heating_length=10e-3, stretching_length=0.2e-3 * 200
    )
    asymmetric_alpha_profile.initialize()


def test_generate_propagation_gif():
    """
    Test generating a propagation GIF from the profile data.

    Parameters
    ----------
    alpha_profile : AlphaProfile
        The profile fixture to use in the test.
    """
    alpha_profile = AlphaProfile(initial_radius=1)

    alpha_profile.add_taper_segment(
        alpha=0, initial_heating_length=3e-3, stretching_length=0.2e-3 * 200
    )
    alpha_profile.initialize()
    generate_propagation_gif(profile=alpha_profile, number_of_frames=10)


if __name__ == "__main__":
    pytest.main([__file__])
