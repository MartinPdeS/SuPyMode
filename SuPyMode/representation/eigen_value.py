# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

from typing import NoReturn
from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot
import matplotlib.pyplot as plt


class EigenValue(InheritFromSuperMode, BaseSingleModePlot):
    """
    Represents the eigenvalues of a mode derived from a supermode in a waveguide or optical fiber simulation.

    This class extends from `InheritFromSuperMode` to access supermode-related data and from `BaseSingleModePlot`
    to provide plotting capabilities tailored to eigenvalue visualization.

    Attributes:
        parent_supermode (InheritFromSuperMode): The parent supermode object from which eigenvalue data is derived.
    """

    def __init__(self, parent_supermode: SuperMode):
        """
        Initializes an EigenValue object with a parent supermode reference.

        Args:
            parent_supermode (InheritFromSuperMode): A reference to the parent supermode object.
        """
        self.parent_supermode = parent_supermode
        self.data = self.parent_supermode.binding.get_eigen_value()

    def _dress_ax(self, ax: plt.Axis) -> NoReturn:
        ax.set_xlabel('Inverse taper ratio')
        ax.set_ylabel('Mode eigen values')


# -
