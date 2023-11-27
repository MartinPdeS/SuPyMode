
import numpy
from PyFiberModes.__future__ import get_normalized_LP_coupling


def integrate_2d(mesh, x, y):
    dx = abs(x[0] - x[1])
    dy = abs(y[0] - y[1])

    integral = numpy.trapz(mesh, dx=dx)

    integral = numpy.trapz(integral, dx=dy)

    return integral


def get_theoretical_mode_coupling(fiber, mode_0, mode_1, itr_list: numpy.ndarray) -> numpy.ndarray:
    coupling_array = numpy.empty(itr_list.shape)

    for idx, itr in enumerate(itr_list):
        _fiber = fiber.scale(itr)

        coupling = get_normalized_LP_coupling(
            fiber=_fiber,
            mode_0=mode_0,
            mode_1=mode_1
        )

        coupling_array[idx] = coupling

    return coupling_array


def get_mode_coupling(superset, mode_0, mode_1) -> numpy.ndarray:
    coordinate = superset.geometry.coordinate_system
    itr_list = superset.itr_list

    gradient = superset.geometry.n2_gradient * coordinate.rho_mesh

    coupling_array = numpy.empty(itr_list.shape)

    for idx, itr in enumerate(itr_list):

        field_01 = superset.LP01.field._data[idx]
        field_02 = superset.LP02.field._data[idx]

        beta_0 = superset.LP01.beta._data[idx]
        beta_1 = superset.LP02.beta._data[idx]

        mesh = field_01 * field_02 * gradient

        integral = integrate_2d(
            mesh=mesh,
            x=coordinate.x_vector * itr,
            y=coordinate.x_vector * itr
        )

        term_0 = abs(beta_0 - beta_1)
        term_1 = numpy.sqrt(beta_0 * beta_1)
        term_2 = 0.5 * superset.wavenumber**2
        term_3 = term_2 / (term_0 * term_1)

        coupling_array[idx] = integral * term_3

    return coupling_array

# -
