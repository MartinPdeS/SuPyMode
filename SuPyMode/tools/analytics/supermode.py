#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
import numpy
from dataclasses import dataclass
from scipy.constants import epsilon_0 as e0, mu_0

# Local imports
from PyFiberModes import Wavelength, FiberFactory, Mode, field as FieldClass


@dataclass
class Supermode():
    mode_number: str
    """ Mode number [starts with LP] """
    fiber: FiberFactory
    """ Fiber type to which compute the mode """
    wavelength: float
    """ Wavelenght to which compute the mode """

    def __post_init__(self):
        l, m = self.mode_number[2:]
        self.l, self.m = int(l), int(m)
        self.mode = Mode('LP', self.l, self.m)
        self.wavenumber = 2 * numpy.pi / self.wavelength

    @property
    def beta(self):
        """ Returns propgation constant of the supermode """
        return self.wavenumber * self.neff

    @property
    def neff(self) -> float:
        """ Returns effective refractive index """
        return self.fiber.neff(self.mode, self.wavelength)

    @property
    def norm_factor(self) -> float:
        """ Returns the norm factor of the supermode """
        factor = 0.5 * numpy.sqrt(e0 / mu_0)

        return factor * (self.beta / self.wavenumber)

    def get_field(self, bound: float, resolution: int = 200) -> numpy.ndarray:
        fields = FieldClass.Field(
            fiber=self.fiber,
            mode=self.mode,
            wl=self.wavelength,
            r=bound,
            np=resolution
        )

        return fields.Ex()

    def get_norm2_field(self, field: numpy.ndarray) -> float:
        """ Returns the L2 norm of the supermode field """
        field = numpy.square(field)
        sum_field = field.sum(axis=0).sum(axis=0)
        return numpy.sqrt(sum_field)

    def get_l2_normalized_field(self, bound: float, resolution: int = 200) -> numpy.ndarray:
        """ Returns a L2 normalized field array """
        field_object = FieldClass.Field(
            fiber=self.fiber,
            mode=self.mode,
            wl=self.wavelength,
            r=bound,
            np=resolution
        )

        field_array = field_object.Ex()

        norm_l2 = self.get_norm2_field(field_array)

        normalized_field = field_array / numpy.sqrt(norm_l2)

        return normalized_field

    def evaluate_field_at_r(self, r_space: numpy.ndarray) -> numpy.ndarray:
        """ Returns array corresponding to the supermode field evaluated at the r_space position """
        wavelength = Wavelength(self.wavelength)

        array = numpy.zeros(r_space.size)

        for idx, r in enumerate(r_space):
            er, hr = self.fiber._rfield(self.mode, wavelength, r)
            array[idx] = er[0]

        return array.squeeze()


# -
