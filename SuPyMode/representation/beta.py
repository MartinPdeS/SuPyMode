# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from SuPyMode.binary.interface_supermode import SUPERMODE
from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot


class Beta(InheritFromSuperMode, BaseSingleModePlot):
    """
    Represents the propagation constants (beta values) of a mode derived from a supermode in optical simulations.

    This class utilizes inheritance from `InheritFromSuperMode` for accessing supermode-related data and
    `BaseSingleModePlot` for plotting functionalities tailored to propagation constant visualization.

    Attributes
    ----------
    supermode : SUPERMODE
        A reference to the parent supermode object from which beta data is sourced.
    data : numpy.ndarray
        The propagation constant (beta) data retrieved from the parent supermode binding.
    """

    def __init__(self, supermode: SUPERMODE):
        """
        Initialize a Beta object with a reference to a parent supermode.

        Parameters
        ----------
        supermode : SUPERMODE
            The parent supermode object that provides the base beta data.
        """
        self.supermode = supermode
        self.data = self.supermode.get_betas()

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels for the propagation constant plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels.
        """
        ax.set_xlabel("Inverse taper ratio")
        ax.set_ylabel("Propagation constant [rad/M]")


# -
