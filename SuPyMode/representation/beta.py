# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot
import matplotlib.pyplot as plt


class Beta(InheritFromSuperMode, BaseSingleModePlot):
    """
    Represents the propagation constants (beta values) of a mode derived from a supermode in optical simulations.

    This class utilizes inheritance from `InheritFromSuperMode` for accessing supermode-related data and
    `BaseSingleModePlot` for plotting functionalities tailored to propagation constant visualization.

    Attributes
    ----------
    parent_supermode : SuperMode
        A reference to the parent supermode object from which beta data is sourced.
    data : numpy.ndarray
        The propagation constant (beta) data retrieved from the parent supermode binding.
    """
    def __init__(self, parent_supermode: SuperMode):
        """
        Initialize a Beta object with a reference to a parent supermode.

        Parameters
        ----------
        parent_supermode : SuperMode
            The parent supermode object that provides the base beta data.
        """
        self.parent_supermode = parent_supermode
        self.data = self.parent_supermode.binding.get_betas()

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels for the propagation constant plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels.
        """
        ax.set_xlabel('Inverse taper ratio')
        ax.set_ylabel('Propagation constant [rad/M]')

# -
