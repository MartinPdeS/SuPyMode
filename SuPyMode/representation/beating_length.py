# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import numpy

from SuPyMode.tools import plot_style
from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot


class BeatingLength(InheritFromSuperMode, BaseMultiModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self.plot_style = plot_style.beating_length

    def get_values(self, other_supermode) -> numpy.ndarray:
        """
        Return the array of the modal coupling for the mode
        """
        return self.parent_supermode.binded_supermode.get_beating_length_with_mode(other_supermode.binded_supermode)


# -
