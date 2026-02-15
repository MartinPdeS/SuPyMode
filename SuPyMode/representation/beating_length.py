# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import numpy
import matplotlib.pyplot as plt
from MPSPlots.helper import post_mpl_plot

from SuPyMode.binary.interface_supermode import SUPERMODE


class BeatingLength:
    """
    Represents the beating lengths between modes of different supermodes in optical fiber simulations.

    This class extends from `InheritFromSuperMode` to utilize supermode-related data and provides plotting functionalities
    tailored to visualize beating length comparisons.

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

    @post_mpl_plot
    def plot(self, other_supermode: SUPERMODE, ax: plt.Axes = None) -> plt.Figure:
        """
        Plot the beating length between the parent supermode and another specified supermode.

        This method generates a plot of the beating length as a function of the inverse taper ratio (ITR),
        formatted according to the predefined plot style.

        Parameters
        ----------
        other_supermode : SUPERMODE
            The supermode to compare against.
        ax : matplotlib.axes.Axes, optional
            The axis on which to plot. If `None`, a new axis is created (default is `None`).

        Returns
        -------
        matplotlib.figure.Figure
            The figure object containing the generated plot.
        """
        if ax is None:
            figure, ax = plt.subplots(nrows=1, ncols=1)
        else:
            figure = ax.figure

        ax.set_xlabel("Inverse taper ratio")
        ax.set_ylabel("Beating length [m]")

        if not self.supermode.is_computation_compatible(other_supermode):
            return

        y = self.get_values(other_supermode=other_supermode)

        label = f"{self.supermode.stylized_label} - {other_supermode.stylized_label}"

        ax.plot(
            self.supermode.model_parameters.itr_list, abs(y), label=label, linewidth=2
        )

        ax.legend()

        return figure


# -
