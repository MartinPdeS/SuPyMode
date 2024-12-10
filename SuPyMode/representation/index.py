# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot
import matplotlib.pyplot as plt


class Index(InheritFromSuperMode, BaseSingleModePlot):
    """
    Represents the effective refractive index of a mode derived from a supermode in optical fiber simulations.

    This class extends from `InheritFromSuperMode` for accessing supermode-related data and `BaseSingleModePlot`
    for plotting functionalities tailored to visualize the effective refractive index.

    Attributes
    ----------
    parent_supermode : SuperMode
        The supermode instance to which this Index object is linked.
    data : numpy.ndarray
        The effective refractive index data for the mode derived from the parent supermode.
    """
    def __init__(self, parent_supermode: SuperMode):
        """
        Initialize an Index object with a reference to a parent supermode.

        Parameters
        ----------
        parent_supermode : SuperMode
            The parent supermode object that provides the base mode data.
        """
        self.parent_supermode = parent_supermode
        self.data = self.parent_supermode.binding.get_index()

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels and limits for the effective refractive index plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels and limits.
        """
        ax.set(
            xlabel='Inverse taper ratio',
            ylabel='Effective refractive index',
            ylim=[1.44, 1.455]
        )

# -
