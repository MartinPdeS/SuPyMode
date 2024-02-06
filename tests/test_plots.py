#!/usr/bin/env python
# -*- coding: utf-8 -*-

from unittest.mock import patch

from SuPyMode.workflow import configuration, Workflow, fiber_catalogue, Boundaries


@patch("matplotlib.pyplot.show")
def test_superset_plot(patch):
    fibers = [
        fiber_catalogue.load_fiber('SMF28', wavelength=1550e-9),
        fiber_catalogue.load_fiber('SMF28', wavelength=1550e-9),
    ]

    fused_structure = configuration.ring.FusedProfile_02x02

    workflow = Workflow(
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

    superset = workflow.superset

    superset.plot('index').show().close()

    superset.plot('beta').show().close()

    superset.plot('eigen-value').show().close()

    superset.plot('normalized-coupling').show().close()

    superset.plot('adiabatic').show().close()

    superset.plot('field').show().close()


@patch("matplotlib.pyplot.show")
def test_representation_plot(patch):
    fibers = [
        fiber_catalogue.load_fiber('SMF28', wavelength=1550e-9),
        fiber_catalogue.load_fiber('SMF28', wavelength=1550e-9),
    ]

    fused_structure = configuration.ring.FusedProfile_02x02

    workflow = Workflow(
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

    superset = workflow.superset

    mode = superset[0]

    mode.index.plot().show().close()

    mode.beta.plot().show().close()

    mode.eigen_value.plot().show().close()

    mode.normalized_coupling.plot(other_supermode=superset[1]).show().close()

    mode.adiabatic.plot(other_supermode=superset[1]).show().close()

    mode.field.plot().show().close()

# -
