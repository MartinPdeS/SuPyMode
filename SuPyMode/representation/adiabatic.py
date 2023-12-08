# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

import numpy

from SuPyMode.representation.base import InheritFromSuperMode, BaseMultiModePlot
from MPSPlots.render2D import SceneList, Axis




ax_style = dict(
    show_legend=True,
    x_label='Inverse taper ratio',
    y_label=r'Adiabatic criterion [$\mu$m$^{-1}$]',
    y_scale='log',
    y_scale_factor=1e-6,
    y_limits=[1e-5, 1],
    line_width=2
)


class Adiabatic(InheritFromSuperMode, BaseMultiModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode

    def get_values(self, other_supermode: SuperMode) -> numpy.ndarray:
        """
        Return the array of the modal coupling for the mode

        :param      other_supermode:  The other supermode
        :type       other_supermode:  SuperMode

        :returns:   The values
        :rtype:     numpy.ndarray
        """
        output = self.parent_supermode.binded_supermode.get_adiabatic_with_mode(other_supermode.binded_supermode)

        if not self.parent_supermode.is_computation_compatible(other_supermode):
            output *= numpy.inf

        return output

    def render_on_ax(self, ax: Axis, other_supermode: SuperMode) -> None:
        """
        Renders the normalized-coupling data with another modes on a given ax

        :param      ax:               Axis to which add the plot
        :type       ax:               Axis
        :param      other_supermode:  The other supermode
        :type       other_supermode:  SuperMode

        :returns:   No returns
        :rtype:     None
        """
        if not self.parent_supermode.is_computation_compatible(other_supermode):
            return

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
