# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from MPSPlots.helper import post_mpl_plot

from SuPyMode.binary.interface_supermode import SUPERMODE


class Adiabatic:
    """
    Represents the adiabatic criterion between modes of different supermodes in optical fiber simulations.

    This class extends from `InheritFromSuperMode` for accessing supermode-related data and provides plotting functionalities
    tailored to visualize adiabatic transition measurements.

    Attributes
    ----------
    supermode : SUPERMODE
        The parent supermode object that provides the base mode data.

    """

    def __init__(self, supermode: SUPERMODE):
        """
        Initialize an Adiabatic object with a reference to a parent supermode.

        Parameters
        ----------
        supermode : SUPERMODE
            The parent supermode object that provides the base mode data.
        """
        self.supermode = supermode

    def get_values(self, other_supermode: SUPERMODE) -> numpy.ndarray:
        """
        Calculate the adiabatic transition measure between the parent supermode and another specified supermode.

        Parameters
        ----------
        other_supermode : SUPERMODE
            The supermode with which to compare the parent supermode.

        Returns
        -------
        numpy.ndarray
            An array of adiabatic transition measures calculated between the two supermodes, adjusted as needed for
            computational compatibility.
        """
        output = self.supermode.get_adiabatic_with_mode(other_supermode)

        if not self.supermode.is_computation_compatible(other_supermode):
            output *= numpy.inf

        return abs(output)

    @post_mpl_plot
    def plot(self, other_supermode: SUPERMODE, ax: plt.Axes = None) -> plt.Figure:
        """
        Plot the adiabatic criterion between the parent supermode and another specified supermode.
        This method generates a plot of the adiabatic criterion as a function of the inverse taper ratio (ITR),
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

        ax.set(
            yscale="log",
            xlabel="Inverse taper ratio",
            ylabel=r"Adiabatic criterion [$\mu$m$^{-1}$]",
            ylim=[1e-5 * 1e6, 1 * 1e6],
        )

        def log_scale_yaxis(value, tick_position):
            return f"{value / 1e6:.0e}"

        # Apply the custom formatter for log scale
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(log_scale_yaxis))

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
