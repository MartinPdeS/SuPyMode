"""Validate scalar LP mode coupling against a radial fiber calculation."""

import numpy
import pytest
import PyFiberModes
from SuPyMode.workflow import (
    BoundaryValue,
    Boundaries,
    DomainAlignment,
    Workflow,
    fiber_loader,
)


def radial_normalized_coupling(*, fiber, modes, radial_samples=400, radial_extent=4):
    r"""Evaluate the step-index interface expression for two LP0m modes.

    This uses the interface expression attributed to Bures, *Optical Fiber
    Theory*, equation 7.39, in the former PyFiberModes implementation
    (commit ff7dd86f^). Unlike that implementation, contributions from all
    interfaces are summed. The scalar radial norm is integrated to four core
    radii; this cutoff is suitable only for the well-confined modes tested here.

    Parameters
    ----------
    fiber : PyFiberModes.Fiber
        Step-index fiber with dimensions and wavelength in meters.
    modes : tuple
        Two axisymmetric scalar LP modes.
    radial_samples : int
        Number of radial integration samples.
    radial_extent : float
        Integration cutoff in units of the core radius.

    Returns
    -------
    float
        Magnitude of the dimensionless normalized coupling.
    """
    radii = numpy.linspace(0.0, radial_extent * fiber.radius, radial_samples)
    norms = []
    propagation_constants = []
    for mode in modes:
        field = numpy.array(
            [
                fiber.get_radial_field(mode=mode, radius=float(radius))[0].rho
                for radius in radii
            ]
        )
        norms.append(2 * numpy.pi * numpy.trapezoid(field**2 * radii, x=radii))
        propagation_constants.append(fiber.get_propagation_constant(mode=mode))

    interface_sum = 0.0
    for inner_layer, outer_layer in fiber.iterate_interfaces():
        radius = inner_layer.radius_out
        fields = [
            fiber.get_radial_field(mode=mode, radius=radius)[0].rho for mode in modes
        ]
        index_difference = (
            inner_layer.refractive_index**2 - outer_layer.refractive_index**2
        )
        interface_sum += (
            2 * numpy.pi * index_difference * radius**2 * numpy.prod(fields)
        )

    beta_0, beta_1 = propagation_constants
    wavenumber = 2 * numpy.pi / fiber.wavelength
    return abs(
        0.5
        * wavenumber**2
        * interface_sum
        / (
            numpy.sqrt(numpy.prod(norms))
            * (beta_0 - beta_1)
            * numpy.sqrt(beta_0 * beta_1)
        )
    )


def test_normalized_coupling():
    """Compare the simulated LP01–LP02 coupling with a radial reference."""
    wavelength = 1550e-9
    fiber_name = "test_multimode_fiber"
    reference_fiber = PyFiberModes.fiber.load_fiber(
        fiber_name=fiber_name, wavelength=wavelength, add_air_layer=False
    )
    cladding_index = reference_fiber.index_list[-1]
    fiber = fiber_loader.load_fiber(
        fiber_name, clad_refractive_index=cladding_index, remove_cladding=False
    )
    workflow = Workflow(
        fiber_list=[fiber],
        wavelength=wavelength,
        resolution=140,
        background_index=cladding_index,
        x_bounds=DomainAlignment.LEFT,
        y_bounds=DomainAlignment.BOTTOM,
        boundaries=[
            Boundaries(right=BoundaryValue.SYMMETRIC, top=BoundaryValue.SYMMETRIC)
        ],
        n_sorted_mode=4,
        n_added_mode=8,
        debug_mode=0,
        auto_label=True,
        itr_final=0.9,
        n_step=12,
    )
    workflow.initialize_geometry()
    workflow.run_solver()

    modes = (PyFiberModes.LP01, PyFiberModes.LP02)
    analytical = numpy.array(
        [
            radial_normalized_coupling(
                fiber=reference_fiber.scale(factor=itr), modes=modes
            )
            for itr in workflow.superset.model_parameters.itr_list
        ]
    )
    refined_reference = radial_normalized_coupling(
        fiber=reference_fiber, modes=modes, radial_samples=1600, radial_extent=8
    )
    numpy.testing.assert_allclose(analytical[0], refined_reference, rtol=1e-3)
    simulation = numpy.abs(
        workflow.superset.LP01.normalized_coupling.get_values(workflow.superset.LP02)
    )
    assert numpy.all(numpy.isfinite(analytical))
    numpy.testing.assert_allclose(simulation, analytical, rtol=0.10, atol=0)


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__]))
