# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import numpy
from MPSPlots.render2D import Axis
from MPSPlots import colormaps
from MPSPlots.render2D import SceneMatrix

from SuPyMode.tools.utils import interpret_slice_number_and_itr, slice_to_itr

from SuPyMode.representation.base import InheritFromSuperMode

ax_style = dict(
    show_legend=False,
    x_label=r'X-Direction [$\mu m$]',
    y_label=r'Y-direction [$\mu m$]',
    x_scale_factor=1e6,
    y_scale_factor=1e6,
    equal=True
)


class Field(InheritFromSuperMode):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self._data = self.parent_supermode.binded_supermode.get_fields()

    def get_values(self):
        return self._data

    def get_norm(self, slice_number: int) -> float:
        return self.parent_supermode.binded_supermode.get_norm(slice_number)

    @property
    def itr_list(self) -> numpy.ndarray:
        return self.parent_supermode.binded_supermode.model_parameters.itr_list

    @property
    def parent_superset(self) -> object:
        return self.parent_supermode.parent_set

    def get_field(
            self,
            slice_number: int = [],
            itr: float = [],
            add_symmetries: bool = True) -> numpy.ndarray:
        """
        Returns the field with the predefined boundary conditions for a certain slice number.
        The normalization type must either be square integration (L2), max value set to one (max), center value set to one (center)
        or the normalization provided with coupled-mode theory (cmt).

        :param      slice_number:    The slice number
        :type       slice_number:    int
        :param      add_symmetries:  Add or not the boundary symmetries
        :type       add_symmetries:  bool

        :returns:   The field mesh.
        :rtype:     numpy.ndarray
        """
        slice_number = numpy.array(slice_number)
        itr = numpy.array(itr)

        assert (slice_number.size == 0) ^ (itr.size == 0), 'Exactly one of the two values [slice_number, itr] has to be defined'

        slice_list, itr_list = interpret_slice_number_and_itr(
            itr_baseline=self.itr_list,
            itr_list=itr,
            slice_list=slice_number
        )

        field = numpy.take(self._data, slice_number, axis=0)

        if add_symmetries:
            field = self._get_symmetrized_field(field=field)

        return field

    def normalize_field(self, field: numpy.ndarray, itr: float, norm_type: str = 'L2') -> numpy.ndarray:
        """
        Deprecated at the moment.

        :param      field:      The field
        :type       field:      { type_description }
        :param      itr:        The itr
        :type       itr:        float
        :param      norm_type:  The normalize type
        :type       norm_type:  str

        :returns:   { description_of_the_return_value }
        :rtype:     { return_type_description }
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
        x_axis, y_axis = self._get_axis_vector(add_symmetries=True)

        field = self._get_symmetrized_field(field=field)

        return x_axis, y_axis, field

    def _get_symmetrized_field(self, field: numpy.ndarray) -> numpy.ndarray:
        """
        Take as input a mesh and return a symmetric version of that mesh.
        The symmetry can be setted mirror the top or bottom, left or right.

        :param      vector:          The 2-d mesh
        :type       vector:          numpy.ndarray

        :returns:   The symmetrized mesh
        :rtype:     numpy.ndarray

        :raises     AssertionError:  Verify that input vector is 2-dimensionnal.
        """
        assert field.ndim == 2, f'Mesh should be 2d, instead {field.ndim} dimensional is provided.'

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
        x, y, field = self._get_symmetrized_field_and_axis(field=self._data[slice])

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
        Plot each of the mode field for different itr value or slice number.

        :param      itr_list:    List of itr value to evaluate the mode field
        :type       itr_list:    list
        :param      slice_list:  List of integer reprenting the slice where the mode field is evaluated
        :type       slice_list:  list

        :returns:   The figure
        :rtype:     SceneMatrix
        """
        figure = SceneMatrix(unit_size=(3, 3))

        slice_list, itr_list = interpret_slice_number_and_itr(
            itr_baseline=self.itr_list,
            itr_list=itr_list,
            slice_list=slice_list
        )

        for n, (itr, slice_number) in enumerate(zip(itr_list, slice_list)):
            ax = figure.append_ax(row=n, column=0)

            ax.set_style(**ax_style)

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
        Gets the title for the plot_field outputed subplots.

        :param      supermode:         The supermode corresponding to the specific subplot.
        :type       supermode:         SuperMode
        :param      itr:               The itr value
        :type       itr:               float
        :param      slice_number:      The slice number
        :type       slice_number:      int
        :param      show_mode_label:   If True the mode label will be shown.
        :type       show_mode_label:   bool
        :param      show_itr:          If True the title contains the itr value.
        :type       show_itr:          bool
        :param      show_slice:        If True the title contains the slice number of the evaluated ITR
        :type       show_slice:        bool

        :returns:   The plot mode field title.
        :rtype:     str
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
