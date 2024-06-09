# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot
from MPSPlots.render2D import SceneList, Axis


class EigenValue(InheritFromSuperMode, BaseSingleModePlot):
    """
    Represents the eigenvalues of a mode derived from a supermode in a waveguide or optical fiber simulation.

    This class extends from `InheritFromSuperMode` to access supermode-related data and from `BaseSingleModePlot`
    to provide plotting capabilities tailored to eigenvalue visualization.

    Attributes:
        parent_supermode (InheritFromSuperMode): The parent supermode object from which eigenvalue data is derived.
    """

    plot_style = dict(
        show_legend=True,
        x_label='Inverse taper ratio',
        y_label='Mode eigen values',
        y_scale="linear",
        line_width=2
    )

    def __init__(self, parent_supermode: SuperMode):
        """
        Initializes an EigenValue object with a parent supermode reference.

        Args:
            parent_supermode (InheritFromSuperMode): A reference to the parent supermode object.
        """
        self.parent_supermode = parent_supermode
        self.data = self.parent_supermode.binding.get_eigen_value()

    def render_on_ax(self, ax: Axis) -> None:
        """
        Renders the eigenvalues as a line plot on the provided Axis object.

        Args:
            ax (Axis): The Axis object on which the eigenvalues will be plotted.

        Note:
            This method utilizes the plotting configuration set on the Axis to define the appearance of the plot.
        """
        ax.add_line(
            x=self.itr_list,
            y=self.data,
            label=f'{self.stylized_label}'
        )

    def plot(self) -> SceneList:
        """
        Generates a plot of the eigenvalues using a SceneList to manage multiple plots if necessary.

        This method creates a single-axis plot showing the mode eigenvalues as a function of the inverse taper ratio.

        Returns:
            SceneList: A scene list containing the eigenvalue plot.
        """
        figure = SceneList()

        ax = figure.append_ax()

        ax.set_style(**self.plot_style)

        self.render_on_ax(ax=ax)

        return figure


# -
