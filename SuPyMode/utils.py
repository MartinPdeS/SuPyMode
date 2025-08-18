#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
from typing import Union, List, Iterable, Callable
if TYPE_CHECKING:
    from SuPyMode.superset import SuperSet
    from SuPyMode.supermode import SuperMode

import numpy
import pickle
from pathlib import Path
from SuPyMode.directories import user_data_directory
from functools import wraps
from pathvalidate import sanitize_filepath
import logging


def get_auto_generated_filename(superset) -> str:
    """
    Generates a filename based on the simulation parameters.

    Returns:
        str: Automatically generated filename.
    """
    fiber_name = "".join(fiber.__class__.__name__ for fiber in superset.fiber_list)
    filename = (
        f"structure={superset.clad_structure.__class__.__name__}_"
        f"{fiber_name}_"
        f"resolution={superset.resolution}_"
        f"wavelength={superset.wavelength}"
    )
    return filename.replace('.', '_')


def parse_filename(save_function: Callable) -> Callable:
    @wraps(save_function)
    def wrapper(superset, filename: str = 'auto', directory: str = 'auto', **kwargs):
        filename = get_auto_generated_filename(superset) if filename == 'auto' else filename

        directory = user_data_directory if directory == 'auto' else directory

        filename = Path(filename)

        filename = sanitize_filepath(filename)

        filename = Path(directory).joinpath(filename)

        save_function(superset, filename=filename, **kwargs)

        logging.info(f"Saving data into: {filename}")

        return filename
    return wrapper


def load_superset(filename: str, directory: str = 'auto'):
    """
    Saves the superset instance as a serialized pickle file.

    :param      filename:  The filename
    :type       filename:  str
    """
    if directory == 'auto':
        directory = user_data_directory

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
    Converts an inverse taper ratio (itr) or a list of itrs to their closest corresponding slice numbers
    in the provided itr_list. If a single itr is provided, returns a single slice number. If a list of
    itrs is provided, returns a list of slice numbers.

    Parameters:
        itr_list (np.ndarray): Array of itrs, typically representing a continuous range of values.
        itr (Union[float, List[float]]): A single itr value or a list of itr values.

    Returns:
        Union[int, List[int]]: The closest slice number or list of slice numbers corresponding to the provided itr or itrs.

    """
    if numpy.isscalar(itr):
        return numpy.argmin(abs(itr_list - itr))

    if isinstance(itr_list, Iterable):
        return [numpy.argmin(abs(itr_list - value)) for value in itr]

    raise TypeError("itr must be a float or a list of floats.")


def slice_to_itr(itr_list: numpy.ndarray, slice_number: int | list) -> float | list:
    """
    Converts a slice number or a list of slice numbers to their corresponding inverse taper ratios (itrs)
    from the provided itr_list. If a single slice number is provided, returns a single itr value.
    If a list of slice numbers is provided, returns a list of itr values.

    Parameters:
        - itr_list (np.ndarray): Array of itrs, typically representing a continuous range of values.
        - slice_number (Union[int, List[int]]): A single slice index or a list of slice indices.

    Returns:
         - Union[float, List[float]]: The itr or list of itrs corresponding to the given slice numbers.

    """
    if numpy.isscalar(slice_number):
        return itr_list[slice_number]

    if isinstance(slice_number, Iterable):
        return numpy.take(itr_list, slice_number)

    raise TypeError("itr_list must be an scalar or a list of scalar.")


def interpret_slice_number_and_itr(
        itr_baseline: numpy.ndarray,
        slice_list: Union[int, List[int], None] = None,
        itr_list: Union[float, List[float], None] = None,
        sort_slice_number: bool = True) -> tuple:
    """
    Interprets slice numbers and corresponding inverse taper ratios (ITRs), returning arrays of slice numbers
    and their respective ITRs based on the provided lists of slices and ITRs.

    Parameters:
        - itr_baseline (numpy.ndarray): Array of baseline ITR values, typically equidistant, representing the full range.
        - slice_list (Union[int, List[int]], optional): A single slice index or a list of slice indices. If None and itr_list is also None, defaults to [0, -1].
        - itr_list (Union[float, List[float]], optional): A single ITR value or a list of ITR values to be converted to slice indices. If provided, slice_list defaults to an empty list unless specified.
        - sort_slice_number (bool, optional): Whether to sort the resulting slice numbers and corresponding ITRs in descending order. Default is False.

    Returns:
        - Tuple[numpy.ndarray, numpy.ndarray]: Two numpy arrays, the first containing slice indices and the second containing the corresponding ITR values.

    Raises:
        - ValueError: If the provided ITRs or slice indices are outside the bounds of the baseline ITR array.
    """

    if slice_list is None and itr_list is None:
        slice_list = [0, -1]
    elif itr_list is not None and slice_list is None:
        slice_list = []

    slice_list = numpy.atleast_1d(slice_list)
    itr_list = numpy.atleast_1d(itr_list) if itr_list is not None else numpy.array([])

    # Logic to convert ITRs to slice indices, and combine lists
    slice_from_itr = itr_to_slice(itr_baseline, itr=itr_list) if itr_list.size > 0 else numpy.array([])
    total_slice_list = numpy.concatenate([slice_list, slice_from_itr]).astype(int)

    total_itr_list = slice_to_itr(itr_baseline, slice_number=total_slice_list)

    # Sort if required
    if sort_slice_number:
        indices = numpy.argsort(total_itr_list)
        total_slice_list = total_slice_list[indices][::-1]  # Descending order
        total_itr_list = total_itr_list[indices][::-1]

    return total_slice_list, total_itr_list


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
    from SuPyMode.supermode import SuperMode
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
