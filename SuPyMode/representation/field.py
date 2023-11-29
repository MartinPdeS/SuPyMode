# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import numpy
from MPSPlots.render2D import Axis
from MPSPlots import colormaps
from MPSPlots.render2D import SceneMatrix

from SuPyMode.tools import plot_style
from SuPyMode.representation.base import InheritFromSuperMode


class Field(InheritFromSuperMode):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self._data = self.parent_supermode.binded_supermode.get_fields()

    def get_values(self):
        return self._data

    def get_norm(self, slice_number: int) -> float:
        return self.parent_supermode.binded_supermode.get_norm(slice_number)

    def get_field(self, slice_number: int = None, itr: float = None, add_symmetries: bool = True, normalization: str = 'L2') -> numpy.ndarray:
        """
        Returns the field with the predefined boundary conditions for a certain slice number.
        The normalization type must either be square integration (L2), max value set to one (max), center value set to one (center)
        or the normalization provided with coupled-mode theory (cmt).

        :param      slice_number:    The slice number
        :type       slice_number:    int
        :param      add_symmetries:  Add or not the boundary symmetries
        :type       add_symmetries:  bool
        :param      normalization:   The normalization type ['L2', 'max', 'center', 'cmt']
        :type       normalization:   str

        :returns:   The field mesh.
        :rtype:     numpy.ndarray
        """
        assert (slice_number is None) ^ (itr is None), 'Exactly one of the two values [slice_number, itr] has to be defined'
        assert normalization.lower() in ['l2', 'max', 'center', 'cmt'], "Normalization type must be in ['l2', 'max', 'center', 'cmt']"

        slice_number, itr = self.parent_supermode.parent_set._interpret_itr_slice_list_(
            slice_list=[] if slice_number is None else slice_number,
            itr_list=[] if itr is None else itr
        )

        field = self._data[slice_number[0]]

        if add_symmetries:
            field = self._get_symmetrized_field(field=field)

        return self.normalize_field(field=field, norm_type=normalization, itr=itr)

    def normalize_field(self, field: numpy.ndarray, itr: float, norm_type: str = 'L2') -> numpy.ndarray:
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

        ax.set_style(**plot_style.field)

    def plot_field(self,
            mode_of_interest: list = 'all',
            itr_list: list[float] = [],
            slice_list: list[int] = [],
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

        slice_list, itr_list = self._interpret_itr_slice_list_(slice_list=slice_list, itr_list=itr_list)

        mode_of_interest = self.interpret_mode_of_interest(mode_of_interest=mode_of_interest)

        for m, mode in enumerate(mode_of_interest):
            for n, (itr, slice_number) in enumerate(zip(itr_list, slice_list)):
                title = self.get_plot_mode_field_title(
                    supermode=mode,
                    itr=itr,
                    slice_number=slice_number,
                    show_mode_label=show_mode_label,
                    show_itr=show_itr,
                    show_slice=show_slice
                )

                ax = figure.append_ax(
                    row=n,
                    column=m,
                    title=title
                )

                field = mode.field.get_field(slice_number=slice_number, add_symmetries=True)

                x, y = mode.field.get_axis(slice_number=slice_number)

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

                ax.set_style(**plot_style.field)

        return figure

    # def plot(self, slice_list: list = [], itr_list: list = []) -> SceneList:
    #     """
    #     Plotting method for the fields.

    #     :param      slice_list:  Value reprenting the slice where the mode field is evaluated.
    #     :type       slice_list:  list
    #     :param      itr_list:    Value of itr value to evaluate the mode field.
    #     :type       itr_list:    list

    #     :returns:   the figure containing all the plots.
    #     :rtype:     SceneList
    #     """
    #     figure = SceneList(unit_size=(3, 3), tight_layout=True, ax_orientation='horizontal')

    #     slice_list, itr_list = self._interpret_itr_slice_list_(
    #         slice_list=slice_list,
    #         itr_list=itr_list
    #     )

    #     slice_list = numpy.atleast_1d(slice_list)
    #     itr_list = numpy.atleast_1d(itr_list)

    #     for n, (slice, itr) in enumerate(zip(slice_list, itr_list)):
    #         ax = figure.append_ax(
    #             title=f'{self.parent_supermode.stylized_label}\n[slice: {slice}  ITR: {itr:.4f}]'
    #         )

    #         self._render_on_ax_(ax=ax, slice=slice)

    #     return figure

# -
