# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

import numpy

from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot
from MPSPlots.render2D import Axis, SceneList


class BeatingLength(InheritFromSuperMode, BaseMultiModePlot):
    """
    Represents the beating lengths between modes of different supermodes in optical fiber simulations.

    This class extends from `InheritFromSuperMode` to utilize supermode-related data and from `BaseMultiModePlot`
    for advanced plotting functionalities tailored to visualize beating length comparisons.

    Class Attributes:
        plot_style (dict): Default style settings for plots generated by this class.
    """

    plot_style = dict(
        show_legend=True,
        x_label='Inverse taper ratio',
        y_label='Beating length [m]',
        y_scale="log",
        line_width=2
    )

    def __init__(self, parent_supermode: SuperMode):
        """
        Initializes a BeatingLength object with a reference to a parent supermode.

        Args:
            parent_supermode (SuperMode): The parent supermode object that provides the base mode data.
        """
        self.parent_supermode = parent_supermode

    def get_values(self, other_supermode: SuperMode) -> numpy.ndarray:
        """
        Calculates the beating length between the parent supermode and another specified supermode.

        Args:
            other_supermode (SuperMode): The supermode with which to compare the parent supermode.

        Returns:
            numpy.ndarray: An array of beating lengths calculated between the two supermodes.
        """
        return self.parent_supermode.binding.get_beating_length_with_mode(other_supermode.binding)

    def render_on_ax(self, ax: Axis, other_supermode: SuperMode) -> None:
        """
        Renders beating length data as a line plot on the provided Axis object, comparing the parent supermode
        with another supermode.

        Args:
            ax (Axis): The Axis object on which to plot the beating lengths.
            other_supermode (SuperMode): The other supermode to compare against.

        Note:
            This method utilizes the `plot_style` class attribute to define the appearance of the plot.
        """
        y = self.get_values(other_supermode=other_supermode)

        label = f'{self.parent_supermode.stylized_label} - {other_supermode.stylized_label}'

        ax.add_line(x=self.itr_list, y=numpy.abs(y), label=label)

    def plot(self, other_supermode: SuperMode) -> SceneList:
        """
        Generates a plot of beating lengths between the parent supermode and another specified supermode using a SceneList.

        This method creates a single-axis plot showing the comparative beating lengths as a function of the inverse taper ratio,
        formatted according to the predefined plot style.

        Args:
            other_supermode (SuperMode): The supermode to compare against.

        Returns:
            SceneList: A scene list containing the plot of beating lengths.
        """
        figure = SceneList()

        ax = figure.append_ax(**self.plot_style)

        self.render_on_ax(ax=ax, other_supermode=other_supermode)

        return figure

# -
