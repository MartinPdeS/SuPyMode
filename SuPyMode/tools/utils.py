#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy


def get_close_points(tolerane: float, y0: numpy.ndarray, y1: numpy.ndarray, x: numpy.ndarray = None):
    idx = numpy.argwhere(numpy.abs(y0 - y1) < tolerane)

    if x is not None:
        return x[idx], y0[idx]

    return y0[idx]


def get_intersection_(y0: numpy.ndarray, y1: numpy.ndarray, x: numpy.ndarray = None):
    idx = numpy.argwhere(numpy.diff(numpy.sign(y0 - y1))).flatten()

    if x is not None:
        return x[idx], y0[idx]

    return y0[idx]


def get_intersection(y0: numpy.ndarray, y1: numpy.ndarray, x: numpy.ndarray, average: bool = True):

    idx = numpy.argwhere(numpy.diff(numpy.sign(y0 - y1))).flatten()

    if not average:
        return y0[idx]

    else:
        x_mean = (x[idx + 1] + x[idx]) / 2
        y_mean = (y0[idx + 1] + y0[idx]) / 2

        return x_mean, y_mean


def test_valid_input(user_input, valid_inputs: list, variable_name: str = ''):
    if user_input not in valid_inputs:
        raise ValueError(f"[{variable_name}] user_input: {user_input} argument not valid. Valid choices are: {valid_inputs}")


# -
