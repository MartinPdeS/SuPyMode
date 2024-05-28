# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot
from MPSPlots.render2D import SceneList, Axis


class Beta(InheritFromSuperMode, BaseSingleModePlot):
    """
    Represents the propagation constants (beta values) of a mode derived from a supermode in optical simulations.

    This class utilizes inheritance from `InheritFromSuperMode` for accessing supermode-related data and
    `BaseSingleModePlot` for plotting functionalities tailored to propagation constant visualization.

    Class Attributes:
        plot_style (dict): A dictionary defining the default style settings for plots generated by this class.

    Attributes:
        parent_supermode (InheritFromSuperMode): A reference to the parent supermode object from which beta data is sourced.
    """

    plot_style = dict(
        show_legend=True,
        x_label='Inverse taper ratio',
        y_label='Propagation constant [rad/M]',
        y_scale="linear",
        line_width=2
    )

    def __init__(self, parent_supermode: SuperMode):
        """
        Initializes a Beta object with a reference to a parent supermode.

        Args:
            parent_supermode (InheritFromSuperMode): The parent supermode object.
        """
        self.parent_supermode = parent_supermode

    @property
    def data(self):
        return self.parent_supermode.binded_supermode.get_betas()

    def render_on_ax(self, ax: Axis) -> None:
        """
        Renders the propagation constants as a line plot on the provided Axis object.

        Args:
            ax (Axis): The Axis object on which to plot the propagation constants.

        Note:
            Utilizes the `plot_style` class attribute to define the appearance of the plot.
        """
        ax.add_line(
            x=self.itr_list,
            y=self.data,
            label=f'{self.stylized_label}'
        )

    def plot(self) -> SceneList:
        """
        Generates a plot of the propagation constants using a SceneList to manage multiple plots if necessary.

        This method creates a single-axis plot showing the propagation constants as a function of the inverse taper ratio,
        formatted according to the predefined plot style.

        Returns:
            SceneList: A scene list containing the plot of propagation constants.
        """
        figure = SceneList()

        ax = figure.append_ax()

        ax.set_style(**self.plot_style)

        self.render_on_ax(ax=ax)

        return figure


# -
