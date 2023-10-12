#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dataclasses import dataclass
from PyFinitDiff.boundaries import Boundaries2D


@dataclass
class ModeLabel:
    boundaries: Boundaries2D
    n_mode: int = 2

    label_list_str = ['01', '11', '21', '02', '31', '12', '41', '22', '03', '51', '32', '13']

    def __post_init__(self):
        self.label_dictionnary = {}

        for mode_number in self.label_list_str:
            azimuthal_number = int(mode_number[0])
            if azimuthal_number == 0:
                self.label_dictionnary[f'LP{mode_number}'] = {'x': 'symmetric', 'y': 'symmetric'}

            elif azimuthal_number % 2 == 1:
                self.label_dictionnary[f'LP{mode_number}_a'] = {'x': 'anti-symmetric', 'y': 'symmetric'}
                self.label_dictionnary[f'LP{mode_number}_b'] = {'x': 'symmetric', 'y': 'anti-symmetric'}

            elif azimuthal_number % 2 == 0:
                self.label_dictionnary[f'LP{mode_number}_a'] = {'x': 'symmetric', 'y': 'symmetric'}
                self.label_dictionnary[f'LP{mode_number}_b'] = {'x': 'anti-symmetric', 'y': 'anti-symmetric'}

    def get_labels(self) -> list:
        label_dictionnary = {
            label: symmetries for label, symmetries in self.label_dictionnary.items() if self.boundaries.x_symmetry in [symmetries['x'], 'zero']

        }

        label_dictionnary = {
            label: symmetries for label, symmetries in label_dictionnary.items() if self.boundaries.y_symmetry in [symmetries['y'], 'zero']

        }

        return list(label_dictionnary.keys())[:self.n_mode]

    def remove_non_degenerate_suffix(self, label_list):
        for label in label_list:
            lp_number = label.split('_')[0]
# -
