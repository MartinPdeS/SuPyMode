#!/usr/bin/env python
# -*- coding: utf-8 -*-

from SuPyMode.binary.Example import get_rho_gradient_5p
from MPSTools.tools.mathematics import gradientO4
import matplotlib.pyplot as plt
import numpy

x_vector = numpy.linspace(-10, 10, 100)
y_vector = numpy.linspace(-10, 10, 100)
dx = x_vector[1] - x_vector[0]
dy = y_vector[1] - y_vector[0]

x_mesh, y_mesh = numpy.meshgrid(x_vector, y_vector)

rho_mesh = numpy.sqrt(numpy.square(x_mesh) + numpy.square(y_mesh))**2


y_gradient, x_gradient = gradientO4(rho_mesh, dx, dy)

theta_mesh = numpy.arctan2(y_mesh, x_mesh)

gradient = (x_gradient * numpy.cos(theta_mesh) + y_gradient * numpy.sin(theta_mesh))

plt.figure()
plt.pcolormesh(gradient)
plt.colorbar()
plt.show()

rho_gradient = get_rho_gradient_5p(
    mesh=rho_mesh,
    x_vector=x_vector,
    y_vector=y_vector
)

plt.figure()
plt.pcolormesh(rho_gradient)
plt.colorbar()
plt.show()

condition = numpy.isclose(rho_gradient, 1, atol=0.5)

assert numpy.all(condition), f"Error in constant (expected value: 1.0) gradient computation. Mean gradient value: {condition.mean()}"


# -
