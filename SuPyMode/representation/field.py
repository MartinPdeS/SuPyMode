# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from SuPyMode.supermode import SuperMode

import numpy
from MPSPlots.render2D import Axis
from MPSPlots import colormaps
from MPSPlots.render2D import SceneMatrix

from SuPyMode.utils import interpret_slice_number_and_itr, slice_to_itr

from SuPyMode.representation.base import InheritFromSuperMode


class Field(InheritFromSuperMode):
    """
    Represents a field derived from a supermode in a modal analysis framework.

    This class extends functionality from a parent supermode class to manage field data operations,
    including retrieving and processing field data for visualization and analysis.

    Attributes:
        parent_supermode (InheritFromSuperMode): Reference to the parent supermode object that provides source data.
    """

    plot_style = dict(
        show_legend=False,
        x_label=r'',
        y_label=r'',
        show_ticks=False,
        x_scale_factor=1e6,
        y_scale_factor=1e6
    )

    def __init__(self, parent_supermode: SuperMode):
        """
        Initialize the Field object with a reference to a parent supermode.

        Args:
            parent_supermode (InheritFromSuperMode): The parent supermode from which this field is derived.
        """
        self.parent_supermode = parent_supermode
        self.data = self.parent_supermode.binding.get_fields()

    def get_norm(self, slice_number: int) -> float:
        """
        Calculate the norm of the field for a specific slice.

        Args:
            slice_number (int): The slice number for which to calculate the norm.

        Returns:
            float: The norm of the field.
        """
        return self.parent_supermode.binding.get_norm(slice_number)

    @property
    def itr_list(self) -> numpy.ndarray:
        """
        Provides a list of iteration indices available for the fields.

        Returns:
            numpy.ndarray: An array of iteration indices.
        """
        return self.parent_supermode.binding.model_parameters.itr_list

    @property
    def parent_superset(self) -> object:
        """
        Access the parent set of the supermode.

        Returns:
            object: The parent set object.
        """
        return self.parent_supermode.parent_set

    def get_field(self, slice_number: int = None, itr: float = None, add_symmetries: bool = True) -> numpy.ndarray:
        """
        Retrieve a specific field adjusted for boundary conditions and optionally add symmetries.

        Args:
            slice_number (int): The slice number to retrieve.
            itr (float): The iteration to use for retrieving the field.
            add_symmetries (bool): Whether to add boundary symmetries to the field.

        Returns:
            numpy.ndarray: The requested field as a numpy array.

        Raises:
            AssertionError: If neither or both of slice_number and itr are defined.
        """
        slice_list, itr_list = interpret_slice_number_and_itr(
            itr_baseline=self.itr_list,
            itr_list=itr,
            slice_list=slice_number
        )

        fields = self.parent_supermode.binding.get_fields()

        fields = numpy.take(fields, slice_list, axis=0)
        if add_symmetries:
            fields = self._get_symmetrized_field(fields)

        return fields

    def normalize_field(self, field: numpy.ndarray, itr: float, norm_type: str = 'L2') -> numpy.ndarray:
        """
        Normalize a field array based on a specified normalization method.

        Currently, this method is deprecated.

        Args:
            field (numpy.ndarray): The field to normalize.
            itr (float): The iteration value for normalization scaling.
            norm_type (str): The type of normalization ('max', 'center', 'L2', 'cmt').

        Returns:
            numpy.ndarray: The normalized field.
        """
        match norm_type.lower():
            case 'max':
                norm = abs(field).max()
            case 'center':
                idx_x_center = numpy.argmin(abs(self.parent_supermode.parent_set.coordinate_system.x_vector))
                idx_y_center = numpy.argmin(abs(self.parent_supermode.parent_set.coordinate_system.y_vector))
                center_value = field[idx_x_center, idx_y_center]
                norm = center_value
            case 'l2':
                dx_scaled = self.parent_supermode.parent_set.coordinate_system.dx * itr
                dy_scaled = self.parent_supermode.parent_set.coordinate_system.dy * itr
                norm = numpy.sqrt(numpy.trapz(numpy.trapz(numpy.square(field), dx=dy_scaled, axis=0), dx=dx_scaled, axis=0))
            case 'cmt':
                dx_scaled = self.parent_supermode.parent_set.coordinate_system.dx * itr
                dy_scaled = self.parent_supermode.parent_set.coordinate_system.dy * itr
                norm = numpy.sqrt(numpy.trapz(numpy.trapz(numpy.square(field), dx=dx_scaled, axis=0), dx=dy_scaled, axis=0))

        return field / norm

    def _get_symmetrized_field_and_axis(self, field: numpy.ndarray) -> tuple:
        """
        Generate a symmetrical version of the input field mesh according to defined boundary conditions.

        Args:
            field (numpy.ndarray): The 2D field mesh to be symmetrized.

        Returns:
            numpy.ndarray: The symmetrized field mesh.

        Raises:
            AssertionError: If the input is not a 2D array.
        """
        x_axis, y_axis = self._get_axis_vector(add_symmetries=True)

        field = self._get_symmetrized_field(field=field)

        return x_axis, y_axis, field

    def _get_symmetrized_field(self, field: numpy.ndarray) -> numpy.ndarray:
        """
        Retrieve the field and axis data adjusted for symmetry.

        This method generates symmetrized versions of the field and its corresponding axis vectors.

        Args:
            field (numpy.ndarray): The field data array to be symmetrized.

        Returns:
            tuple: A tuple containing the x-axis, y-axis, and the symmetrized field data.
        """
        field = field.squeeze()
        assert field.ndim == 2, f"Expected a 2-dimensional array, but got {field.ndim}-dimensional."

        symmetric_field = field[:, -2::-1]
        match self.boundaries.left.lower():
            case 'symmetric':
                field = numpy.c_[symmetric_field, field]

            case 'anti-symmetric':
                field = numpy.c_[-symmetric_field, field]

        match self.boundaries.right.lower():
            case 'symmetric':
                field = numpy.c_[field, symmetric_field]

            case 'anti-symmetric':
                field = numpy.c_[field, -symmetric_field]

        symmetric_field = field[-2::-1, :]
        match self.boundaries.top.lower():
            case 'symmetric':
                field = numpy.r_[field, symmetric_field]

            case 'anti-symmetric':
                field = numpy.r_[field, -symmetric_field]

        match self.boundaries.bottom.lower():
            case 'symmetric':
                field = numpy.r_[symmetric_field, field]

            case 'anti-symmetric':
                field = numpy.r_[-symmetric_field, field]

        return field

    def _render_on_ax_(self, ax: Axis, slice: int) -> None:
        """
        Renders a specified slice of the field data on a provided Axis object.

        This method is typically used to add a field mesh plot to a plotting axis as part of a larger figure or plot matrix.

        Args:
            ax (Axis): The plotting Axis object to render the field on.
            slice (int): The slice index of the field data to render.

        Note:
            This method applies internal colormap settings and adjusts the axis based on the field data.
        """
        field = self.parent_supermode.binding.get_field[slice]

        x, y, field = self._get_symmetrized_field_and_axis(field=field)

        artist = ax.add_mesh(
            x=x,
            y=y,
            scalar=field,
        )

        ax.add_colorbar(
            artist=artist,
            colormap=colormaps.blue_black_red,
            symmetric=True,
            position='right'
        )

    def plot(
            self,
            itr_list: list[float] = [],
            slice_list: list[int] = [0, -1],
            add_symmetries: bool = True,
            show_mode_label: bool = True,
            show_itr: bool = True,
            show_slice: bool = True) -> SceneMatrix:
        """
        Plot the field for specified iterations or slice numbers.

        Args:
            itr_list (list[float]): List of iterations to evaluate the field.
            slice_list (list[int]): List of slices to evaluate the field.
            add_symmetries (bool): Whether to include boundary symmetries in the plot.
            show_mode_label (bool): Whether to show the mode label.
            show_itr (bool): Whether to show the iteration value.
            show_slice (bool): Whether to show the slice number.

        Returns:
            SceneMatrix: A matrix scene containing the plots.
        """
        figure = SceneMatrix(unit_size=(3, 3))

        slice_list, itr_list = interpret_slice_number_and_itr(
            itr_baseline=self.itr_list,
            itr_list=itr_list,
            slice_list=slice_list
        )

        for n, (itr, slice_number) in enumerate(zip(itr_list, slice_list)):
            ax = figure.append_ax(row=n, column=0)

            ax.set_style(**self.plot_style)

            self.render_on_ax(
                ax=ax,
                slice_number=slice_number,
                add_symmetries=add_symmetries
            )

        return figure

    def render_on_ax(
            self,
            ax: object,
            slice_number: int,
            show_mode_label: bool = True,
            show_itr: bool = True,
            show_slice: bool = True,
            add_symmetries: bool = True) -> None:
        """
        Render the mode field at given slice number into input ax.

        :param      ax:               { parameter_description }
        :type       ax:               object
        :param      slice_number:     The slice number
        :type       slice_number:     int
        :param      show_mode_label:  The show mode label
        :type       show_mode_label:  bool
        :param      show_itr:         The show itr
        :type       show_itr:         bool
        :param      show_slice:       The show slice
        :type       show_slice:       bool
        :param      add_symmetries:   Indicates if the symmetries is added
        :type       add_symmetries:   bool

        :returns:   No returns
        :rtype:     None
        """
        ax.title = self.get_plot_mode_field_title(
            slice_number=slice_number,
            show_mode_label=show_mode_label,
            show_itr=show_itr,
            show_slice=show_slice
        )

        field = self.get_field(
            slice_number=slice_number,
            add_symmetries=add_symmetries
        )

        x, y = self.get_axis(
            slice_number=slice_number,
            add_symmetries=add_symmetries
        )

        artist = ax.add_mesh(
            x=x,
            y=y,
            scalar=field,
        )

        ax.add_colorbar(
            artist=artist,
            colormap=colormaps.blue_black_red,
            symmetric=True
        )

    def get_plot_mode_field_title(self, slice_number: int, show_mode_label: bool, show_itr: bool, show_slice: bool) -> str:
        """
        Constructs a title for the field plot based on the mode, iteration, and slice number.

        Args:
            slice_number (int): The slice number for which the field is plotted.
            show_mode_label (bool): Flag to include the mode label in the title.
            show_itr (bool): Flag to include the iteration number in the title.
            show_slice (bool): Flag to include the slice number in the title.

        Returns:
            str: The constructed title for the plot.
        """
        title = ''

        if show_mode_label:
            title += f'{self.stylized_label}'

        if show_itr or show_slice:
            itr = slice_to_itr(itr_list=self.itr_list, slice_number=slice_number)
            title += '\n'

        if show_slice:
            title += f'slice: {slice_number}'

        if show_itr:
            title += f'  itr: {itr:.3f}'

        return title


# -
