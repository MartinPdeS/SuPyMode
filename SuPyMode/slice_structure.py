#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy


class SliceStructure():
    """
    This class describes a coupler slice structure. It includes
    the propagation constant and field structure for each of the
    computed supermodes.
    """
    def __init__(self, parent_superset, itr: float, supermodes: list, add_symmetries: bool = True):
        self.itr = itr
        self.add_symmetries = add_symmetries
        self.slice_value = parent_superset.itr_to_slice(itr)
        self.x, self.y = supermodes[0].get_axis_vector(add_symmetries=self.add_symmetries)
        self.supermodes = supermodes

        for mode in supermodes:
            setattr(self, mode.name, mode)

    @property
    def fields(self) -> list:
        output = [
            mode.field.get_value_from_itr(
                itr=self.itr,
                add_symmetries=self.add_symmetries
            ) for mode in self.supermodes
        ]

        return numpy.asarray(output)

    def get_field_combination(self, amplitudes: list, Linf_normalization: bool = False) -> numpy.ndarray:
        field = numpy.einsum('i, ijk -> jk', abs(amplitudes)**2, self.fields)
        if Linf_normalization:
            field /= numpy.abs(field).max()

        return field
