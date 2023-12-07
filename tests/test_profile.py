#!/usr/bin/env python
# -*- coding: utf-8 -*-

from unittest.mock import patch

from SuPyMode.profiles import AlphaProfile


@patch("matplotlib.pyplot.show")
def test_build_profile(patch):
    profile = AlphaProfile(initial_radius=1)

    profile.add_taper_segment(
        alpha=0,
        initial_heating_length=10e-3,
        stretching_length=0.2e-3 * 200
    )

    profile.plot().show()


@patch("matplotlib.pyplot.show")
def test_build_2_segment_profile(patch):
    profile = AlphaProfile(initial_radius=1)

    profile.add_taper_segment(
        alpha=0,
        initial_heating_length=10e-3,
        stretching_length=0.2e-3 * 200
    )

    profile.add_taper_segment(
        alpha=0,
        initial_heating_length=3e-3,
        stretching_length=0.2e-3 * 200
    )

    profile.plot().show()


@patch("matplotlib.pyplot.show")
def test_build_symmmetric_profile(patch):
    profile = AlphaProfile(initial_radius=1, symmetric=False)

    profile.add_taper_segment(
        alpha=0,
        initial_heating_length=10e-3,
        stretching_length=0.2e-3 * 200
    )

    profile.plot().show()


def test_generate_gif():
    profile = AlphaProfile()

    profile.add_taper_segment(
        alpha=0,
        initial_heating_length=3e-3,
        stretching_length=0.2e-3 * 200
    )

    profile.generate_propagation_gif(number_of_frames=10)

# -
