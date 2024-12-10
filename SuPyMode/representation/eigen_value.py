# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot
import matplotlib.pyplot as plt


class EigenValue(InheritFromSuperMode, BaseSingleModePlot):
    """
    Represents the eigenvalues of a mode derived from a supermode in a waveguide or optical fiber simulation.

    This class extends from `InheritFromSuperMode` to access supermode-related data and from `BaseSingleModePlot`
    to provide plotting capabilities tailored to eigenvalue visualization.

    Attributes
    ----------
    parent_supermode : SuperMode
        The parent supermode object from which eigenvalue data is derived.
    data : numpy.ndarray
        The eigenvalue data retrieved from the parent supermode binding.
    """
    def __init__(self, parent_supermode: SuperMode):
        """
        Initialize an EigenValue object with a parent supermode reference.

        Parameters
        ----------
        parent_supermode : SuperMode
            A reference to the parent supermode object that provides the base eigenvalue data.
        """
        self.parent_supermode = parent_supermode
        self.data = self.parent_supermode.binding.get_eigen_value()

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels for the eigenvalue plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels.
        """
        ax.set_xlabel('Inverse taper ratio')
        ax.set_ylabel('Mode eigenvalues')

# -
