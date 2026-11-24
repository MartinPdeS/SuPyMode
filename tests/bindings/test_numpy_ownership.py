"""Exercise the lifetime of NumPy arrays returned by native code."""

import gc
import numpy
from SuPyMode.model_parameters import ModelParameters


def test_model_gradient_survives_reinitialization_and_owner_deletion():
    """Keep returned arrays alive independently of replaced model data."""
    x_coordinates = numpy.linspace(-2.0, 2.0, 17)
    y_coordinates = numpy.linspace(-1.0, 1.0, 13)
    x_mesh, y_mesh = numpy.meshgrid(x_coordinates, y_coordinates)
    model = ModelParameters(
        mesh=1.4 + 0.01 * x_mesh**2 + 0.02 * y_mesh**2,
        x_vector=x_coordinates,
        y_vector=y_coordinates,
        dx=x_coordinates[1] - x_coordinates[0],
        dy=y_coordinates[1] - y_coordinates[0],
    )
    original_gradient = model.mesh_gradient
    expected_gradient = original_gradient.copy()
    gradient_view = original_gradient[::2, ::2]
    for _ in range(100):
        model.initialize(wavelength=1.55e-6, itr_initial=1.0, itr_final=0.9, n_step=3)
    del model, original_gradient
    gc.collect()
    numpy.testing.assert_array_equal(gradient_view, expected_gradient[::2, ::2])
    del gradient_view
    gc.collect()
