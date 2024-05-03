#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFinitDiff.finite_difference_2D import Boundaries
from dataclasses import dataclass, field
from typing import Dict, List


@dataclass(frozen=True)
class ModeLabel:
    """
    Class for managing mode labels based on boundary conditions and mode numbers.

    Attributes:
        boundaries (Boundaries): The boundary condition object to determine symmetry.
        n_mode (int): The number of mode labels to return, default is 2.
    """
    boundaries: Boundaries
    n_mode: int = 2
    label_list_str: List[str] = field(default_factory=lambda: ['01', '11', '21', '02', '31', '12', '41', '22', '03', '51', '32', '13'])
    label_dictionary: Dict[str, Dict[str, str]] = field(init=False)

    def __post_init__(self):
        object.__setattr__(self, 'label_dictionary', self._initialize_labels())

    def _initialize_labels(self) -> Dict[str, Dict[str, str]]:
        """
        Initializes the label dictionary based on the mode number and their symmetry properties.

        Returns:
            Dict[str, Dict[str, str]]: A dictionary of mode labels with their corresponding x and y symmetry properties.
        """
        label_dict = {}
        for mode_number in self.label_list_str:
            azimuthal_number = int(mode_number[0])
            base_label = f'LP{mode_number}'

            if azimuthal_number == 0:
                label_dict[base_label] = {'x': 'symmetric', 'y': 'symmetric'}
            elif azimuthal_number % 2 == 1:
                label_dict[f'{base_label}_a'] = {'x': 'anti-symmetric', 'y': 'symmetric'}
                label_dict[f'{base_label}_b'] = {'x': 'symmetric', 'y': 'anti-symmetric'}
            else:
                label_dict[f'{base_label}_a'] = {'x': 'symmetric', 'y': 'symmetric'}
                label_dict[f'{base_label}_b'] = {'x': 'anti-symmetric', 'y': 'anti-symmetric'}

        return label_dict

    def get_labels(self) -> List[str]:
        """
        Retrieves a list of mode labels filtered by the boundary conditions and limited by `n_mode`.

        Returns:
            List[str]: A list of mode labels that match the current boundary conditions.
        """
        filtered_labels = [
            label for label, symmetries in self.label_dictionary.items() if self.boundaries.get_x_parity() in [symmetries['x'], 'zero'] and self.boundaries.get_y_parity() in [symmetries['y'], 'zero']
        ]

        return filtered_labels[:self.n_mode]
