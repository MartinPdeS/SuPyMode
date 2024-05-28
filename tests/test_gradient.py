#!/usr/bin/env python
# -*- coding: utf-8 -*-

from SuPyMode.binary.Example import get_rho_gradient_5p

import numpy
import pytest


@pytest.fixture
def x_y_meshes():
    """Fixture to create an mesh and x and y vector."""
    x_vector = numpy.linspace(-10, 10, 100)
    y_vector = numpy.linspace(-10, 10, 100)

    return x_vector, y_vector


@pytest.fixture
def rho_mesh(x_y_vector):
    x_vector, y_vector = x_y_meshes

    x_mesh, y_mesh = numpy.meshgrid(x_vector, y_vector)

    rho_mesh = numpy.sqrt(numpy.square(x_mesh) + numpy.square(y_mesh))

    return rho_mesh


def test_constant_gradient(x_y_meshes, rho_mesh):
    x_vector, y_vector = x_y_meshes

    rho_gradient = get_rho_gradient_5p(
        mesh=rho_mesh,
        x_vector=x_vector,
        y_vector=y_vector
    )

    condition = numpy.isclose(rho_gradient, 1, atol=0.5)
    assert numpy.all(condition), f"Error in constant (expected value: 1.0) gradient computation. Mean gradient value: {condition.mean()}"


def test_null_gradient(rho_mesh_and_vector):
    rho_mesh, x_vector, y_vector = rho_mesh_and_vector

    rho_mesh *= 0.0

    rho_gradient = get_rho_gradient_5p(
        mesh=rho_mesh,
        x_vector=x_vector,
        y_vector=y_vector
    )

    condition = numpy.isclose(rho_gradient, 0, atol=0.5)
    assert numpy.all(condition), f"Error in constant (expected value: 0.0) gradient computation. Mean gradient value: {condition.mean()}"


# -
