#!/usr/bin/env python
# -*- coding: utf-8 -*-

from SuPyMode.binary.interface_mesh import get_rho_gradient_5p

import numpy
import pytest


@pytest.fixture
def x_vector():
    """Fixture to create an mesh and x and y vector."""
    return numpy.linspace(-10, 10, 100)


@pytest.fixture
def y_vector():
    """Fixture to create an mesh and x and y vector."""
    return numpy.linspace(-10, 10, 100)


@pytest.fixture
def rho_mesh(x_vector, y_vector):
    x_mesh, y_mesh = numpy.meshgrid(x_vector, y_vector)

    return numpy.sqrt(numpy.square(x_mesh) + numpy.square(y_mesh))


def test_constant_gradient(x_vector, y_vector, rho_mesh):
    rho_gradient = get_rho_gradient_5p(
        mesh=rho_mesh,
        x_vector=x_vector,
        y_vector=y_vector
    )

    condition = numpy.isclose(rho_gradient, 1, atol=0.5)
    assert numpy.all(condition), f"Error in constant (expected value: 1.0) gradient computation. Mean gradient value: {condition.mean()}"


def test_null_gradient(x_vector, y_vector, rho_mesh):
    rho_mesh *= 0.0

    rho_gradient = get_rho_gradient_5p(
        mesh=rho_mesh,
        x_vector=x_vector,
        y_vector=y_vector
    )

    condition = numpy.isclose(rho_gradient, 0, atol=0.5)
    assert numpy.all(condition), f"Error in constant (expected value: 0.0) gradient computation. Mean gradient value: {condition.mean()}"


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
