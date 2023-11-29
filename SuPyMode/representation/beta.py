# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import numpy

from SuPyMode.tools import plot_style
from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot


class Beta(InheritFromSuperMode, BaseSingleModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self._data = self.parent_supermode.binded_supermode.get_betas()
        self.plot_style = plot_style.beta

    def get_values(self) -> numpy.ndarray:
        return self._data

# -
