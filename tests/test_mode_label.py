#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest

from SuPyMode.workflow import Boundaries
from SuPyMode.mode_label import ModeLabel

configuration_list = [
    dict(boundaries=Boundaries(left='symmetric', top='symmetric'), mode_number=0, name='LP01'),
    dict(boundaries=Boundaries(left='symmetric', top='symmetric'), mode_number=1, name='LP21'),
    dict(boundaries=Boundaries(left='symmetric'), mode_number=1, name='LP11'),
    dict(boundaries=Boundaries(left='symmetric'), mode_number=2, name='LP21'),
    dict(boundaries=Boundaries(), mode_number=5, name='LP02')
]


@pytest.mark.parametrize("configuration", configuration_list)
def test_configurations(configuration: dict):
    mode_name = configuration.pop('name')

    mode_label = ModeLabel(**configuration)

    assert mode_name == mode_label.raw_label, f"Mismatch between expected mode_label for auto-labeler: {mode_name} vs {mode_label.label}"

# -
