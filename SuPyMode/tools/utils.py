#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.superset import SuperSet
    from SuPyMode.supermode import SuperMode

import numpy
import pickle
from pathlib import Path
from SuPyMode.tools.directories import instance_directory


def load_superset(filename: str, directory: str = '.'):
    """
    Saves the superset instance as a serialized pickle file.

    :param      filename:  The filename
    :type       filename:  str
    """
    if directory == 'auto':
        directory = instance_directory

    filename = Path(directory).joinpath(filename).with_suffix('.pickle')

    with open(filename, 'rb') as input_file:
        superset = pickle.load(input_file)

    return superset


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

    if len(idx) == 0:  # No intersection
        return None, None

    if not average:
        return y0[idx]

    else:
        x_mean = (x[idx + 1] + x[idx]) / 2
        y_mean = (y0[idx + 1] + y0[idx]) / 2

        return x_mean, y_mean


def test_valid_input(user_input, valid_inputs: list, variable_name: str = '') -> None:
    if user_input not in valid_inputs:
        raise ValueError(f"[{variable_name}] user_input: {user_input} argument not valid. Valid choices are: {valid_inputs}")


def itr_to_slice(itr_list: numpy.ndarray, itr: float | list) -> int | list:
    """
    Convert itr value to slice number

    :param      itr_list:  The itr list
    :type       itr_list:  numpy.ndarray
    :param      itr:       The itr
    :type       itr:       float | list

    :returns:   The equivalent slice numbers
    :rtype:     float | list
    """
    return_scalar = False

    if numpy.isscalar(itr):
        return_scalar = True
        itr = [itr]

    slice_number = [
        numpy.argmin(abs(itr_list - value)) for value in itr
    ]

    if return_scalar:
        return slice_number[0]

    return slice_number


def slice_to_itr(itr_list: numpy.ndarray, slice_number: int | list) -> float | list:
    """
    Convert itr value to slice number

    :param      itr_list:      The itr list
    :type       itr_list:      numpy.ndarray
    :param      slice_number:  The itr
    :type       slice_number:  int | list

    :returns:   The equivalent slice numbers
    :rtype:     float | list
    """
    return_scalar = False

    if numpy.isscalar(slice_number):
        return_scalar = True
        slice_number = [slice_number]

    itr = itr_list[slice_number]

    if return_scalar:
        return itr[0]

    return itr


def interpret_slice_number_and_itr(
        itr_baseline: numpy.ndarray,
        slice_list: int | list[int] = [],
        itr_list: float | list[float] = [],
        sort_slice_number: bool = True) -> tuple:
    """
    Interpret the itr and slice_list as input and return a tuple containing all the associated values

    :param      itr_baseline:      The itr list
    :type       itr_baseline:      numpy.ndarray
    :param      slice_number:  The slice number
    :type       slice_number:  int | list[int]
    :param      itr:           The itr
    :type       itr:           float | list[float]

    :returns:   All the associated slice numbers and itr values
    :rtype:     tuple
    """
    slice_list = numpy.atleast_1d(slice_list)

    slice_from_itr = itr_to_slice(itr_baseline, itr=itr_list)

    slice_from_itr = numpy.atleast_1d(slice_from_itr)

    total_slice_list = [*slice_list, *slice_from_itr]

    total_itr_list = slice_to_itr(itr_baseline, slice_number=total_slice_list)

    slice_list = numpy.asarray(total_slice_list)

    itr_list = numpy.asarray(total_itr_list)

    if sort_slice_number:
        slice_list = numpy.sort(slice_list)[::-1]

        itr_list = numpy.sort(itr_list)[::-1]

    if len(itr_list) == 1:
        return slice_list[0], itr_list[0]

    return slice_list, itr_list


def interpret_mode_of_interest(superset: SuperSet, mode_of_interest: str | SuperMode | list[SuperMode]) -> list[SuperMode]:
    """
    Interpret and returns the input for the mode_of_intereset argument.

    :param      superset:          The superset
    :type       superset:          SuperSet
    :param      mode_of_interest:  The mode of interest
    :type       mode_of_interest:  str | SuperMode | list[SuperMode]

    :returns:   A list of the mode of interest
    :rtype:     list[SuperMode]
    """
    if isinstance(mode_of_interest, str):
        match mode_of_interest:
            case 'fundamental':
                return superset.fundamental_supermodes
            case 'non-fundamental':
                return superset.non_fundamental_supermodes
            case 'all':
                return superset.supermodes

    if not numpy.iterable(mode_of_interest):
        return [mode_of_interest]

    if isinstance(mode_of_interest, list):
        return mode_of_interest

    raise ValueError(f"Invalid input for {mode_of_interest=}. Valid in put must be string ['fundamental', 'non-fundamental', 'all'] or a list or instance of SuperMode.")

# -
