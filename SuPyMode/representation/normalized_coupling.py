# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

import numpy
from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot
import matplotlib.pyplot as plt


class NormalizedCoupling(InheritFromSuperMode, BaseMultiModePlot):
    """
    Represents the normalized mode coupling between modes of different supermodes in optical fiber simulations.

    This class extends from `InheritFromSuperMode` for accessing supermode-related data and `BaseMultiModePlot`
    for plotting functionalities tailored to visualize mode coupling comparisons.

    Attributes
    ----------
    parent_supermode : SuperMode
        The supermode instance to which this NormalizedCoupling object is linked.

    """
    def __init__(self, parent_supermode: SuperMode):
        """
        Initialize a NormalizedCoupling object with a reference to a parent supermode.

        Parameters
        ----------
        parent_supermode : SuperMode
            The parent supermode object that provides the base mode data.
        """
        self.parent_supermode = parent_supermode

    def get_values(self, other_supermode: SuperMode) -> numpy.ndarray:
        """
        Calculate the normalized mode coupling between the parent supermode and another specified supermode.

        Parameters
        ----------
        other_supermode : SuperMode
            The supermode with which to compare the parent supermode.

        Returns
        -------
        numpy.ndarray
            An array of normalized mode coupling values, adjusted for computational compatibility.
        """
        output = self.parent_supermode.binding.get_normalized_coupling_with_mode(other_supermode.binding)

        if not self.parent_supermode.is_computation_compatible(other_supermode):
            output *= 0

        return output

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels for the normalized coupling plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels.
        """
        ax.set_xlabel('Inverse taper ratio')
        ax.set_ylabel('Mode coupling')

# -
