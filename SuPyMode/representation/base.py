# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

from typing import NoReturn
import matplotlib.pyplot as plt


class BaseMultiModePlot():
    def render_on_ax(self, ax: plt.Axes, other_supermode: SuperMode) -> None:
        """
        Renders normalized mode coupling data as a line plot on the provided Axes object, comparing the parent supermode
        with another supermode.

        Args:
            ax (plt.Axes): The Axes object on which to plot the normalized mode coupling.
            other_supermode (SuperMode): The other supermode to compare against.

        Note:
            This method is conditioned on computational compatibility between the supermodes.
        """
        if not self.parent_supermode.is_computation_compatible(other_supermode):
            return

        y = self.get_values(other_supermode=other_supermode)

        label = f'{self.parent_supermode.stylized_label} - {other_supermode.stylized_label}'

        ax.plot(self.itr_list, abs(y), label=label, linewidth=2)

    def plot(self, other_supermode: SuperMode) -> NoReturn:
        """
        Generates a plot of specific parameter coupling between the parent supermode and another specified supermode using a SceneList.

        This method creates a single-axis plot showing the comparative mode couplings as a function of the inverse taper ratio,
        formatted according to the predefined plot style.

        Args:
            other_supermode (SuperMode): The supermode to compare against.

        Returns:
            SceneList: A scene list containing the plot of normalized mode couplings.
        """
        figure, ax = plt.subplots(1, 1)

        self._dress_ax(ax)

        self.render_on_ax(ax=ax, other_supermode=other_supermode)

        ax.legend()

        figure.tight_layout()
        plt.show()


class BaseSingleModePlot():
    def __getitem__(self, idx: int):
        return self._data[idx]

    def render_on_ax(self, ax: plt.Axis) -> None:
        """
        Renders the eigenvalues as a line plot on the provided Axis object.

        Args:
            ax (Axis): The Axis object on which the eigenvalues will be plotted.

        Note:
            This method utilizes the plotting configuration set on the Axis to define the appearance of the plot.
        """
        ax.plot(self.itr_list, self.data, label=f'{self.stylized_label}', linewidth=2)

    def plot(self) -> NoReturn:
        """
        Generates a plot of the using matplotlib.

        This method creates a single-axis plot showing the propagation constants as a function of the inverse taper ratio,
        formatted according to the predefined plot style.

        """
        figure, ax = plt.subplots(1, 1)

        self._dress_ax(ax)

        self.render_on_ax(ax=ax)

        ax.legend()

        figure.tight_layout()
        plt.show()


class InheritFromSuperMode():
    def _set_axis_(self, ax: plt.Axes):
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
