import matplotlib.pyplot as plt
from MPSPlots.styles import mps

from SuPyMode.binary.interface_boundaries import BoundaryValue
from SuPyMode.binary.interface_supermode import SUPERMODE
from SuPyMode.utils import get_symmetrized_vector


class BaseMultiModePlot:
    def plot(
        self, other_supermode: SUPERMODE, ax: plt.Axes = None, show: bool = True
    ) -> plt.Figure:
        """
        Plot the coupling between the parent supermode and another supermode.

        This method generates a plot of specific parameter couplings between the parent supermode and the specified
        `other_supermode`, using a single-axis matplotlib plot. The plot shows the normalized coupling as a function of
        the inverse taper ratio (ITR), formatted according to the predefined plot style.

        Parameters
        ----------
        other_supermode : SUPERMODE
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

        if not self.supermode.is_computation_compatible(other_supermode):
            return

        y = self.get_values(other_supermode=other_supermode)

        label = f"{self.supermode.stylized_label} - {other_supermode.stylized_label}"

        ax.plot(
            self.supermode.model_parameters.itr_list, abs(y), label=label, linewidth=2
        )

        ax.legend()

        if show:
            plt.show()

        return figure


class BaseSingleModePlot:
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

        """
        if ax is None:
            with plt.style.context(mps):
                figure, ax = plt.subplots(1, 1)
        else:
            figure = ax.figure

        ax.plot(
            self.supermode.model_parameters.itr_list,
            self.data,
            label=f"{self.supermode.stylized_label}",
            linewidth=2,
        )

        ax.legend()

        self._dress_ax(ax)

        if show:
            plt.show()

        return figure


class InheritFromSuperMode:
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
        return self.supermode.parent_set.slice_to_itr(slice)

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
        return self.supermode.parent_set.itr_to_slice(itr)

    def _get_symmetrize_vector(self, *args, **kwargs):
        """
        Get the symmetrization vector from the parent supermode.

        Returns
        -------
        Any
            The symmetrization vector computed by the parent supermode.
        """
        return self.supermode._get_symmetrize_vector(*args, **kwargs)

    def _get_axis_vector(self, supermode: object, add_symmetries: bool = True) -> tuple:
        """
        Computes the full axis vectors, optionally including symmetries.

        Parameters
        ----------
        add_symmetries : bool, optional, default=True
            Whether to include symmetries when computing the axis vectors.

        Returns
        -------
        tuple
            A tuple containing the full x-axis and y-axis vectors.
        """
        full_x_axis = self.supermode.model_parameters.x_vector
        full_y_axis = self.supermode.model_parameters.y_vector

        if not add_symmetries:
            return full_x_axis, full_y_axis

        if supermode.boundaries.right in [
            BoundaryValue.Symmetric,
            BoundaryValue.AntiSymmetric,
        ]:
            full_x_axis = get_symmetrized_vector(full_x_axis, symmetry_type="last")
            full_x_axis.sort()

        if supermode.boundaries.left in [
            BoundaryValue.Symmetric,
            BoundaryValue.AntiSymmetric,
        ]:
            full_x_axis = get_symmetrized_vector(full_x_axis, symmetry_type="first")
            full_x_axis.sort()

        if supermode.boundaries.top in [
            BoundaryValue.Symmetric,
            BoundaryValue.AntiSymmetric,
        ]:
            full_y_axis = get_symmetrized_vector(full_y_axis, symmetry_type="last")
            full_y_axis.sort()

        if supermode.boundaries.bottom in [
            BoundaryValue.Symmetric,
            BoundaryValue.AntiSymmetric,
        ]:
            full_y_axis = get_symmetrized_vector(full_y_axis, symmetry_type="first")
            full_y_axis.sort()

        return full_x_axis, full_y_axis

    def get_axis(
        self, supermode: object, slice_number: int, add_symmetries: bool = True
    ) -> tuple:
        """
        Computes the scaled axis vectors for a specific slice, optionally including symmetries.

        Parameters
        ----------
        slice_number : int
            The slice index for which to compute the axis vectors.
        add_symmetries : bool, optional, default=True
            Whether to include symmetries in the computed axis vectors.

        Returns
        -------
        tuple
            A tuple containing the scaled x-axis and y-axis vectors for the given slice.
        """
        itr = supermode.model_parameters.itr_list[slice_number]
        x_axis, y_axis = self._get_axis_vector(
            supermode=supermode, add_symmetries=add_symmetries
        )
        return (x_axis * itr, y_axis * itr)


# -
