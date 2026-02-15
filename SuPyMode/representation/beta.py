# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from MPSPlots.helper import post_mpl_plot

from SuPyMode.binary.interface_supermode import SUPERMODE


class Beta:
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

        ax.set_xlabel("Inverse taper ratio")
        ax.set_ylabel("Propagation constant [rad/M]")

        return figure


# -
