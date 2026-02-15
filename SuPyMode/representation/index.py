# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from MPSPlots.helper import post_mpl_plot

from SuPyMode.binary.interface_supermode import SUPERMODE


class Index:
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

    @post_mpl_plot
    def plot(self, ax: plt.Axes = None, show: bool = True) -> plt.Figure:
        """
        Plot the propagation constant for a single mode.

        This method generates a plot of the propagation constants as a function of the inverse taper ratio (ITR),
        formatted according to the predefined plot style.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The axis on which to plot. If `None`, a new axis is created (default is `None`).
        show : bool, optional
            Whether to display the plot immediately (default is `True`).

        Returns
        -------
        matplotlib.figure.Figure
            The figure object containing the generated plot.

        """
        if ax is None:
            figure, ax = plt.subplots(1, 1)
        else:
            figure = ax.figure

        ax.plot(
            self.supermode.model_parameters.itr_list,
            self.data,
            label=f"{self.supermode.stylized_label}",
            linewidth=2,
        )

        ax.legend()

        ax.set(
            xlabel="Inverse taper ratio",
            ylabel="Effective refractive index",
            ylim=[1.44, 1.455],
        )

        return figure


# -
