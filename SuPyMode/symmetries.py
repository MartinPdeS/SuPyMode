#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyFinitDiff.boundaries import Boundaries2D

mode_dictionary = {
    0: {'label': 'LP01', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},

    1: {'label': 'LP11_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'anti-symmetric'},
    2: {'label': 'LP11_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'symmetric'},

    3: {'label': 'LP21_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},
    4: {'label': 'LP21_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'anti-symmetric'},

    5: {'label': 'LP02', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},

    6: {'label': 'LP31_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'anti-symmetric'},
    7: {'label': 'LP31_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'symmetric'},

    8: {'label': 'LP12_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'anti-symmetric'},
    9: {'label': 'LP12_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'symmetric'},

    10: {'label': 'LP41_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},
    11: {'label': 'LP41_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'anti-symmetric'},

    12: {'label': 'LP22_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},
    13: {'label': 'LP22_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'anti-symmetric'},

    14: {'label': 'LP03', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},

    15: {'label': 'LP51_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'anti-symmetric'},
    16: {'label': 'LP51_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'symmetric'},

    17: {'label': 'LP23_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},
    18: {'label': 'LP23_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'anti-symmetric'},

    19: {'label': 'LP61_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},
    20: {'label': 'LP61_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'anti-symmetric'},

    21: {'label': 'LP13_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'anti-symmetric'},
    22: {'label': 'LP13_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'symmetric'},

    23: {'label': 'LP24_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},
    24: {'label': 'LP24_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'anti-symmetric'},

    25: {'label': 'LP71_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'anti-symmetric'},
    26: {'label': 'LP71_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'symmetric'},

    27: {'label': 'LP32_a', 'x_symmetry': 'symmetric', 'y_symmetry': 'anti-symmetric'},
    28: {'label': 'LP32_b', 'x_symmetry': 'anti-symmetric', 'y_symmetry': 'symmetric'},

    29: {'label': 'LP04', 'x_symmetry': 'symmetric', 'y_symmetry': 'symmetric'},

}


def get_modes_label(boundaries: Boundaries2D, n_mode: int = None) -> list:
    sub_mode = mode_dictionary

    if boundaries.x_symmetry != 'zero':
        sub_mode = {
            k: v for k, v in sub_mode.items() if v['x_symmetry'] == boundaries.x_symmetry
        }

    if boundaries.y_symmetry != 'zero':
        sub_mode = {
            k: v for k, v in sub_mode.items() if v['y_symmetry'] == boundaries.y_symmetry
        }

    return [v['label'] for k, v in sub_mode.items()][:n_mode]


# -
