#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFinitDiff.finite_difference_2D import Boundaries


mode_dict = [
    {'mode': (0, 1, ""), 'x': 'symmetric', 'y': 'symmetric'},
    {'mode': (1, 1, r"_a"), 'x': 'symmetric', 'y': 'anti-symmetric'},
    {'mode': (1, 1, "_b"), 'x': 'anti-symmetric', 'y': 'symmetric'},
    {'mode': (2, 1, "_a"), 'x': 'symmetric', 'y': 'symmetric'},
    {'mode': (2, 1, "_b"), 'x': 'anti-symmetric', 'y': 'anti-symmetric'},
    {'mode': (0, 2, ""), 'x': 'symmetric', 'y': 'symmetric'},
    {'mode': (3, 1, "_a"), 'x': 'symmetric', 'y': 'anti-symmetric'},
    {'mode': (3, 1, "_b"), 'x': 'anti-symmetric', 'y': 'symmetric'},
    {'mode': (1, 2, "_a"), 'x': 'symmetric', 'y': 'anti-symmetric'},
    {'mode': (1, 2, "_b"), 'x': 'anti-symmetric', 'y': 'symmetric'},
    {'mode': (4, 1, "_a"), 'x': 'symmetric', 'y': 'symmetric'},
    {'mode': (4, 1, "_b"), 'x': 'anti-symmetric', 'y': 'anti-symmetric'},
    {'mode': (2, 2, "_a"), 'x': 'symmetric', 'y': 'symmetric'},
    {'mode': (2, 2, "_b"), 'x': 'anti-symmetric', 'y': 'anti-symmetric'},
    {'mode': (0, 3, ""), 'x': 'symmetric', 'y': 'symmetric'},
    {'mode': (5, 1, "_a"), 'x': 'symmetric', 'y': 'anti-symmetric'},
    {'mode': (5, 1, "_b"), 'x': 'anti-symmetric', 'y': 'symmetric'},
    {'mode': (3, 2, "_a"), 'x': 'symmetric', 'y': 'anti-symmetric'},
    {'mode': (3, 2, "_b"), 'x': 'anti-symmetric', 'y': 'symmetric'},
    {'mode': (1, 3, "_a"), 'x': 'symmetric', 'y': 'anti-symmetric'},
    {'mode': (1, 3, "_b"), 'x': 'anti-symmetric', 'y': 'symmetric'},
    {'mode': (6, 1, "_a"), 'x': 'symmetric', 'y': 'symmetric'},
    {'mode': (6, 1, "_b"), 'x': 'anti-symmetric', 'y': 'anti-symmetric'},
    {'mode': (4, 2, "_a"), 'x': 'symmetric', 'y': 'symmetric'},
    {'mode': (4, 2, "_b"), 'x': 'anti-symmetric', 'y': 'anti-symmetric'},
]


class ModeLabel:
    """
    A class to represent the LP mode label of an optical fiber based on boundary conditions and mode number.

    Attributes:
        boundaries (Boundaries): The boundary conditions for the mode.
        mode_number (int): The mode number.
        x_parity (str): The parity in the x direction (symmetric, anti-symmetric, or zero).
        y_parity (str): The parity in the y direction (symmetric, anti-symmetric, or zero).
        azimuthal (int): The azimuthal mode number.
        radial (int): The radial mode number.
        sub_label (str): The sub-label for the mode.
    """

    def __init__(self, boundaries: Boundaries, mode_number: int):
        """
        Initializes the ModeLabel with given boundary conditions and mode number.

        Args:
            boundaries (Boundaries): The boundary conditions.
            mode_number (int): The mode number.
        """
        self.boundaries = boundaries
        self.mode_number = mode_number

        self.initialize()

    def initialize(self) -> None:
        self.x_parity = self._get_x_parity()
        self.y_parity = self._get_y_parity()

        filtered_modes = self.get_filtered_mode_list()

        if self.mode_number >= len(filtered_modes):
            self.azimuthal, self.radial, self.sub_label = None, None, None
            self.raw_label = f"Mode{self.mode_number}"
            self.label = self.raw_label

        else:
            self.azimuthal, self.radial, self.sub_label = filtered_modes[self.mode_number]
            self.raw_label = f"LP{self.azimuthal}{self.radial}"
            self.label = f"{self.raw_label}{self.sub_label}"

    def get_filtered_mode_list(self) -> list[tuple]:
        return [m['mode'] for m in mode_dict if self.x_parity in [m['x'], 'zero'] and self.y_parity in [m['y'], 'zero']]

    def _get_x_parity(self) -> str:
        """
        Determines the parity in the x direction based on boundary conditions.

        Returns:
            str: The x parity (symmetric, anti-symmetric, or zero).
        """
        if self.boundaries.left == 'symmetric' or self.boundaries.right == 'symmetric':
            return 'symmetric'
        elif self.boundaries.left == 'anti-symmetric' or self.boundaries.right == 'anti-symmetric':
            return 'anti-symmetric'
        else:
            return 'zero'

    def _get_y_parity(self) -> str:
        """
        Determines the parity in the y direction based on boundary conditions.

        Returns:
            str: The y parity (symmetric, anti-symmetric, or zero).
        """
        if self.boundaries.top == 'symmetric' or self.boundaries.bottom == 'symmetric':
            return 'symmetric'
        elif self.boundaries.top == 'anti-symmetric' or self.boundaries.bottom == 'anti-symmetric':
            return 'anti-symmetric'
        else:
            return 'zero'

    def __repr__(self) -> str:
        """
        Returns the string representation of the ModeLabel.

        Returns:
            str: The LP mode label.
        """
        return self.label

    def __str__(self) -> str:
        """
        Returns the string representation of the ModeLabel.

        Returns:
            str: The LP mode label.
        """
        return self.__repr__()

# -
