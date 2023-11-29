# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import numpy

from SuPyMode.tools import plot_style
from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot


class Adiabatic(InheritFromSuperMode, BaseMultiModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self.plot_style = plot_style.adiabatic

    def get_values(self, other_supermode) -> numpy.ndarray:
        """
        Return the array of the modal coupling for the mode
        """
        output = self.parent_supermode.binded_supermode.get_adiabatic_with_mode(other_supermode.binded_supermode)

        if not self.parent_supermode.is_computation_compatible(other_supermode):
            output *= numpy.inf

        return output


# -
