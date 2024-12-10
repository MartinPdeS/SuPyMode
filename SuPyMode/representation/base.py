from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

import matplotlib.pyplot as plt
from MPSPlots.styles import mps


class BaseMultiModePlot():
    def plot(self, other_supermode: SuperMode, ax: plt.Axes = None, show: bool = True) -> plt.Figure:
        """
        Plot the coupling between the parent supermode and another supermode.

        This method generates a plot of specific parameter couplings between the parent supermode and the specified
        `other_supermode`, using a single-axis matplotlib plot. The plot shows the normalized coupling as a function of
        the inverse taper ratio (ITR), formatted according to the predefined plot style.

        Parameters
        ----------
        other_supermode : SuperMode
            The supermode to compare against.
        ax : matplotlib.axes.Axes, optional
            The axis on which to plot. If `None`, a new axis is created (default is `None`).
        show : bool, optional
            Whether to display the plot immediately (default is `True`).

        Returns
        -------
        matplotlib.figure.Figure
            The figure object containing the generated plot.

        Examples
        --------
        >>> fig, ax = plt.subplots()
        >>> base_multi_mode_plot.plot(other_supermode=mode2, ax=ax, show=True)
        >>> plt.show()
        """
        if ax is None:
            with plt.style.context(mps):
                figure, ax = plt.subplots(1, 1)
        else:
            figure = ax.figure

        self._dress_ax(ax)

        if not self.parent_supermode.is_computation_compatible(other_supermode):
            return

        y = self.get_values(other_supermode=other_supermode)

        label = f'{self.parent_supermode.stylized_label} - {other_supermode.stylized_label}'

        ax.plot(self.itr_list, abs(y), label=label, linewidth=2)

        ax.legend()

        if show:
            plt.show()

        return figure


class BaseSingleModePlot():
    def __getitem__(self, idx: int):
        return self._data[idx]

    def plot(self, ax: plt.Axes = None, show: bool = True) -> plt.Figure:
        """
        Plot the propagation constant for a single mode.

        This method generates a plot of the propagation constants as a function of the inverse taper ratio (ITR),
        formatted according to the predefined plot style.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            The axis on which to plot. If `None`, a new axis is created (default is `None`).
        show : bool, optional
            Whether to display the plot immediately (default is `True`).

        Returns
        -------
        matplotlib.figure.Figure
            The figure object containing the generated plot.

        Examples
        --------
        >>> fig, ax = plt.subplots()
        >>> base_single_mode_plot.plot(ax=ax, show=True)
        >>> plt.show()
        """
        if ax is None:
            with plt.style.context(mps):
                figure, ax = plt.subplots(1, 1)
        else:
            figure = ax.figure

        ax.plot(self.itr_list, self.data, label=f'{self.stylized_label}', linewidth=2)

        ax.legend()

        self._dress_ax(ax)

        if show:
            plt.show()

        return figure


class InheritFromSuperMode():
    def _set_axis_(self, ax: plt.Axes):
        """
        Set the axis properties according to the predefined plot style.

        This method applies various properties (e.g., labels, limits) to the given axis based on the `plot_style` dictionary.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axis to which the plot style will be applied.
        """
        for element, value in self.plot_style.items():
            setattr(ax, element, value)

    def __getitem__(self, idx):
        """
        Get the data value at the specified index.

        Parameters
        ----------
        idx : int
            The index of the data value to retrieve.

        Returns
        -------
        Any
            The data value at the specified index.
        """
        return self.data[idx]

    @property
    def mode_number(self) -> int:
        """
        Get the mode number of the parent supermode.

        Returns
        -------
        int
            The mode number of the parent supermode.
        """
        return self.parent_supermode.mode_number

    @property
    def solver_number(self) -> int:
        """
        Get the solver number of the parent supermode.

        Returns
        -------
        int
            The solver number of the parent supermode.
        """
        return self.parent_supermode.solver_number

    @property
    def axes(self):
        """
        Get the axes of the parent supermode.

        Returns
        -------
        Any
            The axes associated with the parent supermode.
        """
        return self.parent_supermode.axes

    @property
    def boundaries(self):
        """
        Get the boundary conditions of the parent supermode.

        Returns
        -------
        Any
            The boundaries of the parent supermode.
        """
        return self.parent_supermode.boundaries

    @property
    def itr_list(self):
        """
        Get the list of inverse taper ratio (ITR) values.

        Returns
        -------
        list
            The list of ITR values associated with the parent supermode.
        """
        return self.parent_supermode.itr_list

    @property
    def ID(self):
        """
        Get the identifier (ID) of the parent supermode.

        Returns
        -------
        Any
            The identifier of the parent supermode.
        """
        return self.parent_supermode.ID

    @property
    def label(self):
        """
        Get the label of the parent supermode.

        Returns
        -------
        str
            The label of the parent supermode.
        """
        return self.parent_supermode.label

    @property
    def stylized_label(self):
        """
        Get the stylized label of the parent supermode.

        Returns
        -------
        str
            The stylized label of the parent supermode.
        """
        return self.parent_supermode.stylized_label

    def slice_to_itr(self, slice: list = []):
        """
        Convert slice indices to inverse taper ratio (ITR) values.

        Parameters
        ----------
        slice : list of int, optional
            A list of slice indices to convert (default is an empty list).

        Returns
        -------
        list
            A list of ITR values corresponding to the provided slice indices.
        """
        return self.parent_supermode.parent_set.slice_to_itr(slice)

    def itr_to_slice(self, itr: list = []):
        """
        Convert inverse taper ratio (ITR) values to slice indices.

        Parameters
        ----------
        itr : list of float, optional
            A list of ITR values to convert (default is an empty list).

        Returns
        -------
        list
            A list of slice indices corresponding to the provided ITR values.
        """
        return self.parent_supermode.parent_set.itr_to_slice(itr)

    def _get_symmetrize_vector(self, *args, **kwargs):
        """
        Get the symmetrization vector from the parent supermode.

        Returns
        -------
        Any
            The symmetrization vector computed by the parent supermode.
        """
        return self.parent_supermode._get_symmetrize_vector(*args, **kwargs)

    def _get_axis_vector(self, *args, **kwargs):
        """
        Get the axis vector from the parent supermode.

        Returns
        -------
        Any
            The axis vector computed by the parent supermode.
        """
        return self.parent_supermode._get_axis_vector(*args, **kwargs)

    def get_axis(self, *args, **kwargs):
        """
        Get the axis information from the parent supermode.

        Returns
        -------
        Any
            The axis information retrieved from the parent supermode.
        """
        return self.parent_supermode.get_axis(*args, **kwargs)

# -
