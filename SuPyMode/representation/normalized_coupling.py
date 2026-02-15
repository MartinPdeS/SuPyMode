# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot as plt

from SuPyMode.binary.interface_supermode import SUPERMODE
from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot


class NormalizedCoupling(InheritFromSuperMode, BaseMultiModePlot):
    """
    Represents the normalized mode coupling between modes of different supermodes in optical fiber simulations.

    This class extends from `InheritFromSuperMode` for accessing supermode-related data and `BaseMultiModePlot`
    for plotting functionalities tailored to visualize mode coupling comparisons.

    Attributes
    ----------
    supermode : SUPERMODE
        The supermode instance to which this NormalizedCoupling object is linked.

    """

    def __init__(self, supermode: SUPERMODE):
        """
        Initialize a NormalizedCoupling object with a reference to a parent supermode.

        Parameters
        ----------
        supermode : SUPERMODE
            The parent supermode object that provides the base mode data.
        """
        self.supermode = supermode

    def get_values(self, other_supermode: SUPERMODE) -> numpy.ndarray:
        """
        Calculate the normalized mode coupling between the parent supermode and another specified supermode.

        Parameters
        ----------
        other_supermode : SUPERMODE
            The supermode with which to compare the parent supermode.

        Returns
        -------
        numpy.ndarray
            An array of normalized mode coupling values, adjusted for computational compatibility.
        """
        output = self.supermode.get_normalized_coupling_with_mode(other_supermode)

        if not self.supermode.is_computation_compatible(other_supermode):
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
        ax.set_xlabel("Inverse taper ratio")
        ax.set_ylabel("Mode coupling")


# -
