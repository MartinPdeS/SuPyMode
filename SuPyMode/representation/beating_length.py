# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import numpy
import matplotlib.pyplot as plt
from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot
from SuPyMode.binary.interface_supermode import SUPERMODE


class BeatingLength(InheritFromSuperMode, BaseMultiModePlot):
    """
    Represents the beating lengths between modes of different supermodes in optical fiber simulations.

    This class extends from `InheritFromSuperMode` to utilize supermode-related data and from `BaseMultiModePlot`
    for advanced plotting functionalities tailored to visualize beating length comparisons.

    Attributes
    ----------
    supermode : SUPERMODE
        The parent supermode object that provides the base mode data.
    """

    def __init__(self, supermode: SUPERMODE):
        """
        Initialize a BeatingLength object with a reference to a parent supermode.

        Parameters
        ----------
        supermode : SUPERMODE
            The parent supermode object that provides the base mode data.
        """
        self.supermode = supermode

    def get_values(self, other_supermode: SUPERMODE) -> numpy.ndarray:
        """
        Calculate the beating length between the parent supermode and another specified supermode.

        Parameters
        ----------
        other_supermode : SUPERMODE
            The supermode with which to compare the parent supermode.

        Returns
        -------
        numpy.ndarray
            An array of beating lengths calculated between the two supermodes.
        """
        return self.supermode.get_beating_length_with_mode(other_supermode)

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels for the beating length plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels.
        """
        ax.set_xlabel("Inverse taper ratio")
        ax.set_ylabel("Beating length [m]")


# -
