# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
import numpy

from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot
from MPSPlots.render2D import Axis, SceneList

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode


ax_style = dict(
    show_legend=True,
    x_label='Inverse taper ratio',
    y_label='Beating length [m]',
    y_scale="log",
    line_width=2
)


class BeatingLength(InheritFromSuperMode, BaseMultiModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode

    def get_values(self, other_supermode) -> numpy.ndarray:
        """
        Return the array of the modal coupling for the mode
        """
        return self.parent_supermode.binded_supermode.get_beating_length_with_mode(other_supermode.binded_supermode)

    def render_on_ax(self, ax: Axis, other_supermode: SuperMode) -> None:
        """
        Renders the normalized-coupling data with another modes on a given ax

        :param      ax:               The axis to which add the plot
        :type       ax:               Axis
        :param      other_supermode:  The other supermode
        :type       other_supermode:  SuperMode

        :returns:   No returns
        :rtype:     None
        """
        y = self.get_values(other_supermode=other_supermode)

        ax.add_line(
            x=self.itr_list,
            y=numpy.abs(y),
            label=f'{self.parent_supermode.stylized_label} - {other_supermode.stylized_label}'
        )

    def plot(self, other_supermode: SuperMode) -> SceneList:
        """
        Plots the normalized coupling of this specific mode with the other one
        given as input.

        :param      other_supermode:  The other supermode
        :type       other_supermode:  SuperMode

        :returns:   The scene list.
        :rtype:     SceneList
        """
        figure = SceneList()

        ax = figure.append_ax(**ax_style)

        self.render_on_ax(ax=ax, other_supermode=other_supermode)

        return figure

# -
