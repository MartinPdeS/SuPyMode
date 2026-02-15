# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
from MPSPlots.helper import post_mpl_plot

from SuPyMode.binary.interface_supermode import SUPERMODE


class EigenValue:
    """
    Represents the eigenvalues of a mode derived from a supermode in a waveguide or optical fiber simulation.

    This class extends from `InheritFromSuperMode` to access supermode-related data and from `BaseSingleModePlot`
    to provide plotting capabilities tailored to eigenvalue visualization.

    Attributes
    ----------
    supermode : SUPERMODE
        The parent supermode object from which eigenvalue data is derived.
    data : numpy.ndarray
        The eigenvalue data retrieved from the parent supermode binding.
    """

    def __init__(self, supermode: SUPERMODE):
        """
        Initialize an EigenValue object with a parent supermode reference.

        Parameters
        ----------
        supermode : SUPERMODE
            A reference to the parent supermode object that provides the base eigenvalue data.
        """
        self.supermode = supermode
        self.data = self.supermode.get_eigen_value()

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
            ylabel="Eigenvalue",
        )

        return figure


# -
