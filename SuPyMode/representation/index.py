# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

from SuPyMode.binary.interface_supermode import SUPERMODE
from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot


class Index(InheritFromSuperMode, BaseSingleModePlot):
    """
    Represents the effective refractive index of a mode derived from a supermode in optical fiber simulations.

    This class extends from `InheritFromSuperMode` for accessing supermode-related data and `BaseSingleModePlot`
    for plotting functionalities tailored to visualize the effective refractive index.

    Attributes
    ----------
    supermode : SUPERMODE
        The supermode instance to which this Index object is linked.
    data : numpy.ndarray
        The effective refractive index data for the mode derived from the parent supermode.
    """

    def __init__(self, supermode: SUPERMODE):
        """
        Initialize an Index object with a reference to a parent supermode.

        Parameters
        ----------
        supermode : SUPERMODE
            The parent supermode object that provides the base mode data.
        """
        self.supermode = supermode
        self.data = self.supermode.get_index()

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels and limits for the effective refractive index plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels and limits.
        """
        ax.set(
            xlabel="Inverse taper ratio",
            ylabel="Effective refractive index",
            ylim=[1.44, 1.455],
        )


# -
