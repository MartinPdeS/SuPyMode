# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

import numpy
from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot
import matplotlib.pyplot as plt


class BeatingLength(InheritFromSuperMode, BaseMultiModePlot):
    """
    Represents the beating lengths between modes of different supermodes in optical fiber simulations.

    This class extends from `InheritFromSuperMode` to utilize supermode-related data and from `BaseMultiModePlot`
    for advanced plotting functionalities tailored to visualize beating length comparisons.

    Attributes
    ----------
    parent_supermode : SuperMode
        The parent supermode object that provides the base mode data.
    """
    def __init__(self, parent_supermode: SuperMode):
        """
        Initialize a BeatingLength object with a reference to a parent supermode.

        Parameters
        ----------
        parent_supermode : SuperMode
            The parent supermode object that provides the base mode data.
        """
        self.parent_supermode = parent_supermode

    def get_values(self, other_supermode: SuperMode) -> numpy.ndarray:
        """
        Calculate the beating length between the parent supermode and another specified supermode.

        Parameters
        ----------
        other_supermode : SuperMode
            The supermode with which to compare the parent supermode.

        Returns
        -------
        numpy.ndarray
            An array of beating lengths calculated between the two supermodes.
        """
        return self.parent_supermode.binding.get_beating_length_with_mode(other_supermode.binding)

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels for the beating length plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels.
        """
        ax.set_xlabel('Inverse taper ratio')
        ax.set_ylabel('Beating length [m]')

# -
