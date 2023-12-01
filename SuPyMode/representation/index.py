# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import numpy

from SuPyMode.representation.base import InheritFromSuperMode, BaseSingleModePlot
from MPSPlots.render2D import SceneList, Axis


ax_style = dict(
    show_legend=True,
    x_label='Inverse taper ratio',
    y_label='Effective refraction index',
    y_scale="linear",
    y_limits=[1.44, 1.455],
    line_width=2
)


class Index(InheritFromSuperMode, BaseSingleModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self._data = self.parent_supermode.binded_supermode.get_index()

    def get_values(self) -> numpy.ndarray:
        return self._data

    def render_on_ax(self, ax: Axis) -> None:
        """
        Render the instance values to given ax.

        :param      ax:   The axis to which add the values
        :type       ax:   Axis

        :returns:   No returns
        :rtype:     None
        """
        ax.add_line(
            x=self.itr_list,
            y=self._data,
            label=f'{self.stylized_label}'
        )

    def plot(self) -> SceneList:
        """
        Plot the mode effective index.

        :returns:   The figure
        :rtype:     SceneMatrix
        """
        figure = SceneList()

        ax = figure.append_ax()

        ax.set_style(**ax_style)

        self.render_on_ax(ax=ax)

        return figure


# -
