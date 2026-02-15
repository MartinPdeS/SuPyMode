# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
import matplotlib.pyplot as plt

from SuPyMode.binary.interface_supermode import SUPERMODE
from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot


class EigenValue(InheritFromSuperMode, BaseSingleModePlot):
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

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels for the eigenvalue plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels.
        """
        ax.set_xlabel("Inverse taper ratio")
        ax.set_ylabel("Mode eigenvalues")


# -
