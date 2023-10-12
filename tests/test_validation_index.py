#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy

# from SuPyMode.tools.fibermodes_validation import get_SuPyMode_smf28_taper, get_fibermodes_indexes


# def test_validation_fibermodes_index():
#     wavelength = 1550e-9

#     superset = get_SuPyMode_smf28_taper(
#         wavelength=wavelength,
#         resolution=80,
#         n_step=100
#     )

#     fibermodes_effective_indexes = get_fibermodes_indexes(
#         mode_number_list=['LP01', 'LP02', 'LP21'],
#         itr_list=superset.itr_list,
#         wavelength=wavelength,
#     )

#     for fibermode_name, supymode_name in zip(['LP01', 'LP02', 'LP21'], ['LP01', 'LP02', 'LP21_a']):
#         fibermode_neff = fibermodes_effective_indexes[fibermode_name]
#         supymode_neff = getattr(superset, supymode_name).index.get_values()

#         discrepency = numpy.isclose(fibermode_neff, supymode_neff, rtol=1e-3)

#         if not numpy.all(discrepency):
#             raise ValueError('Mismatch between fibermodes and SuPyMode for effective index')

# -
