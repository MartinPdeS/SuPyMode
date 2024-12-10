# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

import numpy
from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


class Adiabatic(InheritFromSuperMode, BaseMultiModePlot):
    """
    Represents the adiabatic criterion between modes of different supermodes in optical fiber simulations.

    This class extends from `InheritFromSuperMode` for accessing supermode-related data and from `BaseMultiModePlot`
    for plotting functionalities tailored to visualize adiabatic transition measurements.

    Attributes
    ----------
    parent_supermode : SuperMode
        The parent supermode object that provides the base mode data.

    """
    def __init__(self, parent_supermode: SuperMode):
        """
        Initialize an Adiabatic object with a reference to a parent supermode.

        Parameters
        ----------
        parent_supermode : SuperMode
            The parent supermode object that provides the base mode data.
        """
        self.parent_supermode = parent_supermode

    def get_values(self, other_supermode: SuperMode) -> numpy.ndarray:
        """
        Calculate the adiabatic transition measure between the parent supermode and another specified supermode.

        Parameters
        ----------
        other_supermode : SuperMode
            The supermode with which to compare the parent supermode.

        Returns
        -------
        numpy.ndarray
            An array of adiabatic transition measures calculated between the two supermodes, adjusted as needed for
            computational compatibility.
        """
        output = self.parent_supermode.binding.get_adiabatic_with_mode(other_supermode.binding)

        if not self.parent_supermode.is_computation_compatible(other_supermode):
            output *= numpy.inf

        return abs(output)

    def _dress_ax(self, ax: plt.Axes) -> None:
        """
        Set axis labels for the adiabatic criterion plot.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis object on which to set the labels.
        """
        ax.set(
            yscale='log',
            xlabel='Inverse taper ratio',
            ylabel=r'Adiabatic criterion [$\mu$m$^{-1}$]',
            ylim=[1e-5 * 1e6, 1 * 1e6],
        )

        def log_scale_yaxis(value, tick_position):
            return f'{value / 1e6:.0e}'

        # Apply the custom formatter for log scale
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(log_scale_yaxis))

# -
