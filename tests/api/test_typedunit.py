"""Validate unit contracts before native conversions and Pydantic coercion."""

import numpy as np
import pint
import pytest
import TypedUnit
from SuPyMode.model_parameters import ModelParameters
from SuPyMode.taper import AlphaProfile
from SuPyMode import ureg
from SuPyMode.solver import SuPySolver
from SuPyMode.workflow import Workflow


def test_shared_registry_and_workflow_preserve_quantities() -> None:
    """Typed quantities survive configuration without being silently stripped."""
    assert ureg is TypedUnit.ureg
    wavelength = 1550 * ureg.nanometer
    rotation = 15 * ureg.degree
    workflow = Workflow(
        wavelength=wavelength, clad_rotation=rotation, boundaries=[], resolution=20
    )
    assert workflow.wavelength is wavelength
    assert workflow.clad_rotation is rotation
    assert workflow.wavelength.units == ureg.nanometer


@pytest.mark.parametrize(
    "value",
    [True, np.nan, 1.4 * ureg.dimensionless, 1.4 * pint.UnitRegistry().dimensionless],
)
def test_workflow_rejects_nonplain_refractive_indices(value: object) -> None:
    """Pydantic must not coerce dimensionless quantities to optical scalars."""
    with pytest.raises(ValueError, match="background_index"):
        Workflow(
            wavelength=1.55e-6, boundaries=[], resolution=20, background_index=value
        )


@pytest.mark.parametrize(
    "value",
    [0, np.nan, 1 * ureg.second, [1.55e-6], 1550 * pint.UnitRegistry().nanometer],
)
def test_invalid_wavelengths_are_rejected(value: object) -> None:
    """Dimensionality, shape, registry, finite values and positivity are checked."""
    with pytest.raises(ValueError):
        Workflow(wavelength=value, boundaries=[], resolution=20)


def test_native_model_matches_si_inputs() -> None:
    """Converting coordinate and wavelength units preserves native model values."""
    x = np.linspace(-2, 2, 7)
    y = np.linspace(-1, 1, 5)
    mesh = np.full((5, 7), 1.45)
    physical = ModelParameters(
        mesh=mesh,
        x_vector=x * ureg.micrometer,
        y_vector=y * ureg.micrometer,
        dx=(x[1] - x[0]) * ureg.micrometer,
        dy=(y[1] - y[0]) * ureg.micrometer,
    )
    legacy = ModelParameters(
        mesh=mesh,
        x_vector=x * 1e-6,
        y_vector=y * 1e-6,
        dx=(x[1] - x[0]) * 1e-6,
        dy=(y[1] - y[0]) * 1e-6,
    )
    physical.initialize(
        wavelength=1550 * ureg.nanometer, itr_initial=1.0, itr_final=0.8, n_step=3
    )
    legacy.initialize(wavelength=1.55e-6, itr_initial=1.0, itr_final=0.8, n_step=3)
    assert physical.wavenumber == pytest.approx(legacy.wavenumber)
    assert physical.dx == pytest.approx(legacy.dx)
    np.testing.assert_allclose(physical.mesh_gradient, legacy.mesh_gradient)


@pytest.mark.parametrize(
    "coordinates",
    [np.arange(3) * ureg.second, np.arange(3) * pint.UnitRegistry().meter],
)
def test_solver_rejects_invalid_coordinate_quantities(coordinates: object) -> None:
    """Solver validation happens before stripping registry and unit information."""
    with pytest.raises(ValueError):
        SuPySolver(
            mesh=np.full((3, 3), 1.45), x=coordinates, y=np.arange(3) * ureg.micrometer
        )


def test_solver_supports_explicit_coordinate_and_wavelength_units() -> None:
    """The Python solver delegates validated magnitudes to the native model."""
    coordinates = np.linspace(-1, 1, 5) * ureg.micrometer
    solver = SuPySolver(mesh=np.full((5, 5), 1.45), x=coordinates, y=coordinates)
    solver.init_superset(wavelength=1550 * ureg.nanometer, n_step=3)
    assert solver.model_parameters.wavelength == pytest.approx(1.55e-6)
    assert coordinates.units == ureg.micrometer


def test_taper_quantity_inputs_match_legacy_si_profile() -> None:
    """Explicit lengths reproduce the original taper arrays after SI conversion."""
    physical = AlphaProfile(initial_radius=62.5 * ureg.micrometer)
    legacy = AlphaProfile(initial_radius=62.5e-6)
    physical.add_taper_segment(
        alpha=0.0,
        initial_heating_length=10 * ureg.millimeter,
        stretching_length=20 * ureg.millimeter,
    )
    legacy.add_taper_segment(
        alpha=0.0, initial_heating_length=0.01, stretching_length=0.02
    )
    physical.initialize()
    legacy.initialize()
    np.testing.assert_allclose(physical.distance, legacy.distance)
    np.testing.assert_allclose(physical.radius, legacy.radius)
    np.testing.assert_allclose(physical.itr_list, legacy.itr_list)
    with pytest.raises(ValueError, match="alpha"):
        physical.add_taper_segment(
            alpha=0.1 * ureg.dimensionless,
            initial_heating_length=0.01,
            stretching_length=0.02,
        )
