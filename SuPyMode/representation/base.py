# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode


import numpy
from MPSPlots.render2D import SceneList, Axis


class BaseMultiModePlot():
    def _render_on_ax_(self, ax: Axis, other_supermode: SuperMode = None):
        if other_supermode is None:
            other_supermode = self.parent_supermode.parent_set.supermodes
        else:
            other_supermode = numpy.atleast_1d(other_supermode)

        for mode in other_supermode:
            if mode.ID == self.ID or mode.solver_number != self.solver_number:
                continue

            ax.add_line(
                x=self.itr_list,
                y=self.get_values(mode),
                label=f'{self.stylized_label} - {mode.stylized_label}'
            )

            ax.set_style(**self.plot_style)

    def plot(
            self,
            other_supermode: SuperMode = None,
            row: int = 0,
            col: int = 0) -> None:
        """
        Plotting method for the index.

        :param      slice_list:  Value reprenting the slice where the mode field is evaluated.
        :type       slice_list:  list
        :param      itr_list:    Value of itr value to evaluate the mode field.
        :type       itr_list:    list

        :returns:   the figure containing all the plots.
        :rtype:     SceneList
        """
        figure = SceneList(unit_size=(10, 4), tight_layout=True)

        ax = figure.append_ax()

        self._render_on_ax_(ax=ax, other_supermode=other_supermode)

        return figure


class BaseSingleModePlot():
    def __getitem__(self, idx):
        return self._data[idx]

    def _render_on_ax_(self, ax: Axis) -> None:
        self._set_axis_(ax)

        ax.set_style(self.plot_style)

        ax.add_line(x=self.itr_list, y=self._data, label=self.stylized_label)

    def plot(self, row: int = 0, col: int = 0) -> None:
        """
        Plotting method for the index.

        :param      slice_list:  Value reprenting the slice where the mode field is evaluated.
        :type       slice_list:  list
        :param      itr_list:    Value of itr value to evaluate the mode field.
        :type       itr_list:    list

        :returns:   the figure containing all the plots.
        :rtype:     SceneList
        """
        figure = SceneList(unit_size=(10, 4), tight_layout=True)

        ax = figure.append_ax()

        self._render_on_ax_(ax)

        return figure


class InheritFromSuperMode():
    def _set_axis_(self, ax: Axis):
        for element, value in self.plot_style.items():
            setattr(ax, element, value)

    def __getitem__(self, idx):
        return self._data[idx]

    @property
    def mode_number(self) -> int:
        return self.parent_supermode.mode_number

    @property
    def solver_number(self) -> int:
        return self.parent_supermode.solver_number

    @property
    def axes(self):
        return self.parent_supermode.axes

    @property
    def boundaries(self):
        return self.parent_supermode.boundaries

    @property
    def itr_list(self):
        return self.parent_supermode.itr_list

    @property
    def ID(self):
        return self.parent_supermode.ID

    @property
    def label(self):
        return self.parent_supermode.label

    @property
    def stylized_label(self):
        return self.parent_supermode.stylized_label

    def slice_to_itr(self, slice: list = []):
        return self.parent_supermode.parent_set.slice_to_itr(slice)

    def itr_to_slice(self, itr: list = []):
        return self.parent_supermode.parent_set.itr_to_slice(itr)

    def _get_symmetrize_vector(self, *args, **kwargs):
        return self.parent_supermode._get_symmetrize_vector(*args, **kwargs)

    def _get_axis_vector(self, *args, **kwargs):
        return self.parent_supermode._get_axis_vector(*args, **kwargs)

    def get_axis(self, *args, **kwargs):
        return self.parent_supermode.get_axis(*args, **kwargs)

# -
