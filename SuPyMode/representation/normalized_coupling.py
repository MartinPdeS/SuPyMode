# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot as plt
from MPSPlots.helper import post_mpl_plot

from SuPyMode.binary.interface_supermode import SUPERMODE


class NormalizedCoupling:
    """
    Represents the normalized mode coupling between modes of different supermodes in optical fiber simulations.

    This class extends from `InheritFromSuperMode` for accessing supermode-related data and provides plotting functionalities
    tailored to visualize mode coupling comparisons.

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

    @post_mpl_plot
    def plot(self, other_supermode: SUPERMODE, ax: plt.Axes = None) -> plt.Figure:
        """
        Plot the normalized mode coupling between the parent supermode and another specified supermode.

        This method generates a plot of the normalized coupling as a function of the inverse taper ratio (ITR),
        formatted according to the predefined plot style.

        Parameters
        ----------
        other_supermode : SUPERMODE
            The supermode to compare against.
        ax : matplotlib.axes.Axes
            The axis on which to plot.

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
        ax.set_ylabel("Mode coupling")

        if not self.supermode.is_computation_compatible(other_supermode):
            return

        y = self.get_values(other_supermode=other_supermode)

        label = f"{self.supermode.stylized_label} - {other_supermode.stylized_label}"

        ax.plot(self.supermode.model_parameters.itr_list, abs(y), label=label)

        ax.legend()

        return figure


# -
