#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFinitDiff.finite_difference_2D import Boundaries
from enum import Enum


class Parity(Enum):
    SYMMETRIC = 'symmetric'
    ANTI_SYMMETRIC = 'anti-symmetric'
    ZERO = 'zero'


mode_dict = []


for azimuthal, radial in [(0, 1), (1, 1), (2, 1), (0, 2), (3, 1), (1, 2), (4, 1), (2, 2), (0, 3), (5, 1), (3, 2), (1, 3), (6, 1), ]:

    if azimuthal == 0:
        parities = [[Parity.SYMMETRIC.value, Parity.SYMMETRIC.value]]
        sublabels = ['']

    elif azimuthal % 2 == 1:
        parities = [
            [Parity.ANTI_SYMMETRIC.value, Parity.SYMMETRIC.value],
            [Parity.SYMMETRIC.value, Parity.ANTI_SYMMETRIC.value],
        ]
        sublabels = ['_a', '_b']

    elif azimuthal % 2 == 0:
        parities = [
            [Parity.SYMMETRIC.value, Parity.SYMMETRIC.value],
            [Parity.ANTI_SYMMETRIC.value, Parity.ANTI_SYMMETRIC.value],
        ]
        sublabels = ['_a', '_b']

    for (x_parity, y_parity), sublabel in zip(parities, sublabels):

        mode = dict(mode=(azimuthal, radial, sublabel), x=x_parity, y=y_parity)

        mode_dict.append(mode)


class ModeLabel:
    """
    Represents the LP mode label of an optical fiber based on boundary conditions and mode number.

    The `ModeLabel` class encapsulates information about the optical mode based on the boundary conditions
    of an optical fiber. It calculates and assigns labels such as `LP01`, `LP11_a`, etc., depending on the
    given boundary conditions and mode number.

    Parameters
    ----------
    boundaries : Boundaries
        The boundary conditions for the mode, indicating symmetries in different directions.
    mode_number : int
        The mode number to label, corresponding to the specific optical mode.

    Attributes
    ----------
    boundaries : Boundaries
        The boundary conditions for the mode.
    mode_number : int
        The mode number to label.
    x_parity : str
        The parity in the x direction (symmetric, anti-symmetric, or zero).
    y_parity : str
        The parity in the y direction (symmetric, anti-symmetric, or zero).
    azimuthal : int or None
        The azimuthal mode number, extracted from the filtered mode list.
    radial : int or None
        The radial mode number, extracted from the filtered mode list.
    sub_label : str or None
        The sub-label for the mode, such as `_a` or `_b`.
    raw_label : str
        The base LP label without the sub-label.
    label : str
        The full LP label for the mode, including the sub-label.
    """

    def __init__(self, boundaries: Boundaries, mode_number: int):
        """
        Initializes the `ModeLabel` instance with given boundary conditions and mode number.

        Parameters
        ----------
        boundaries : Boundaries
            The boundary conditions indicating symmetries in the left, right, top, and bottom directions.
        mode_number : int
            The mode number, used to identify the optical mode.
        """
        self.boundaries = boundaries
        self.mode_number = mode_number

        self.initialize()

    def initialize(self) -> None:
        """
        Initialize and calculate the mode label based on boundary conditions and mode number.

        This method sets the `x_parity` and `y_parity` based on the boundary conditions, filters the
        mode list, and assigns appropriate labels for azimuthal and radial modes.
        """
        self.x_parity = self._get_parity(self.boundaries.left, self.boundaries.right)
        self.y_parity = self._get_parity(self.boundaries.top, self.boundaries.bottom)

        filtered_modes = self.get_filtered_mode_list()

        if self.mode_number >= len(filtered_modes):
            raise ValueError(f"Mode number {self.mode_number} exceeds available modes. Max allowed: {len(filtered_modes) - 1}")

        self.azimuthal, self.radial, self.sub_label = filtered_modes[self.mode_number]
        self.raw_label = f"LP{self.azimuthal}{self.radial}"
        self.label = f"{self.raw_label}{self.sub_label}"

    def get_filtered_mode_list(self) -> list[tuple]:
        """
        Filters the list of available modes based on the x and y parity conditions.

        Returns
        -------
        list of tuple
            A list of filtered modes that match the specified x and y parity conditions.
        """
        return [m['mode'] for m in mode_dict if self.x_parity in [m['x'], Parity.ZERO.value] and self.y_parity in [m['y'], Parity.ZERO.value]]

    def _get_parity(self, boundary_1: str, boundary_2: str) -> str:
        """
        Determines the parity in a direction based on boundary conditions.

        Parameters
        ----------
        boundary_1 : str
            The boundary condition for one side (e.g., left or top).
        boundary_2 : str
            The boundary condition for the opposite side (e.g., right or bottom).

        Returns
        -------
        str
            The parity ('symmetric', 'anti-symmetric', or 'zero').
        """
        if Parity.SYMMETRIC.value in (boundary_1, boundary_2):
            return Parity.SYMMETRIC.value
        elif Parity.ANTI_SYMMETRIC.value in (boundary_1, boundary_2):
            return Parity.ANTI_SYMMETRIC.value
        else:
            return Parity.ZERO.value

    def __repr__(self) -> str:
        """
        Returns the string representation of the `ModeLabel`.

        Returns
        -------
        str
            The LP mode label, including the azimuthal and radial information.
        """
        return self.label

    def __str__(self) -> str:
        """
        Returns the string representation of the `ModeLabel`.

        Returns
        -------
        str
            The LP mode label.
        """
        return self.__repr__()

# -
