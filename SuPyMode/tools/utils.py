#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
from typing import Iterable
if TYPE_CHECKING:
    from SuPyMode.superset import SuperSet
    from SuPyMode.supermode import SuperMode

import numpy
import pickle
from pathlib import Path
from SuPyMode.directories import instance_directory


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
        sort_slice_number: bool = False) -> tuple:
    """
    Interprets slice numbers and corresponding inverse taper ratios (ITRs), returning arrays of slice numbers
    and their respective ITRs based on the provided lists of slices and ITRs.

    Parameters:
    - itr_baseline (numpy.ndarray): Array of baseline ITR values, typically equidistant, representing the full range.
    - slice_list (Union[int, List[int]], optional): A single slice index or a list of slice indices. Default is an empty list.
    - itr_list (Union[float, List[float]], optional): A single ITR value or a list of ITR values to be converted to slice indices. Default is an empty list.
    - sort_slice_number (bool, optional): Whether to sort the resulting slice numbers and corresponding ITRs in descending order. Default is False.

    Returns:
    - Tuple[numpy.ndarray, numpy.ndarray]: Two numpy arrays, the first containing slice indices and the second containing the corresponding ITR values.

    Raises:
    - ValueError: If the provided ITRs or slice indices are outside the bounds of the baseline ITR array.

    Note:
    - This function assumes that the input ITRs and slice numbers are within the range covered by the `itr_baseline`.
    - The function will also combine and sort the slice indices derived directly and those converted from the given ITRs if `sort_slice_number` is True.
    """
    return_as_iterable = False

    if isinstance(slice_list, Iterable) or isinstance(itr_list, Iterable):
        return_as_iterable = True

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

    if not return_as_iterable:
        return slice_list[0], itr_list[0]

    return slice_list, itr_list


def interpret_mode_of_interest(superset: SuperSet, mode_of_interest: str | SuperMode | list[SuperMode]) -> list[SuperMode]:
    """
    Resolves the mode of interest from user input to the appropriate list of SuperMode instances
    based on the specified criteria or direct references.

    Parameters:
        - superset (SuperSet): The superset containing all supermodes, including fundamental and non-fundamental modes.
        - mode_of_interest (Union[str, SuperMode, List[SuperMode]]): This parameter can be a string specifying a category
          of modes such as 'fundamental', 'non-fundamental', 'all', a single SuperMode instance, or a list of SuperMode instances.

    Returns:
        - List[SuperMode]: A list of SuperMode instances corresponding to the specified mode of interest.

    Raises:
        - ValueError: If the mode_of_interest is not one of the expected types or if the string input does not match
        any known category.
    """
    if isinstance(mode_of_interest, str):
        match mode_of_interest:
            case 'fundamental':
                return superset.fundamental_supermodes
            case 'non-fundamental':
                return superset.non_fundamental_supermodes
            case 'all':
                return superset.supermodes
            case _:
                raise ValueError(f"Unrecognized mode category '{mode_of_interest}'. Expected 'fundamental', 'non-fundamental', or 'all'.")

    if isinstance(mode_of_interest, SuperMode):
        return [mode_of_interest]

    if isinstance(mode_of_interest, list) and all(isinstance(item, SuperMode) for item in mode_of_interest):
        return mode_of_interest

    raise ValueError("mode_of_interest must be either 'fundamental', 'non-fundamental', 'all', a SuperMode instance, or a list of SuperMode instances.")


def get_symmetrized_vector(vector: numpy.ndarray, symmetry_type: str = 'last') -> numpy.ndarray:
    """
    Generate a symmetric version of the input vector based on the specified symmetry type.

    Parameters:
    -----------
    vector : numpy.ndarray
        A one-dimensional array for which the symmetric version is to be calculated.
    symmetry_type : str, optional
        Type of symmetry to apply. Supported types:
        - 'last': Symmetrize using the last element as reference.
        - 'first': Symmetrize using the first element as reference.
        Default is 'last'.

    Returns:
    --------
    numpy.ndarray
        A new vector that is the symmetrized version of the input vector.

    Raises:
    -------
    ValueError
        If the input vector is not one-dimensional or the symmetry type is unsupported.
    """

    if vector.ndim != 1:
        raise ValueError(f"Expected a 1-dimensional vector, but got a {vector.ndim}-dimensional vector instead.")

    if symmetry_type.lower() not in ['last', 'first']:
        raise ValueError("Symmetry type must be 'last' or 'first'.")

    size = len(vector)
    dx = numpy.diff(vector)[0]  # More robust than assuming vector[1] - vector[0]

    if symmetry_type.lower() == 'last':
        start_value = vector[-1]
        expanded = numpy.arange(0, 2 * size - 1) * dx
        return start_value - expanded[::-1] if dx > 0 else start_value + expanded[::-1]
    else:  # 'first'
        start_value = vector[0]
        expanded = numpy.arange(0, 2 * size - 1) * dx
        return start_value + expanded

# -
