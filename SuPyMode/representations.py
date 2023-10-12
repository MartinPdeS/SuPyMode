# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

import numpy
from MPSPlots.render2D import SceneList, Axis
from MPSPlots import colormaps

from SuPyMode.tools import plot_style


class BaseMultiModePlot():
    def _render_on_ax_(self, ax: Axis, other_supermode: 'SuperMode' = None):
        if other_supermode is None:
            other_supermode = self.parent_supermode.parent_set.supermodes

        for mode in other_supermode:
            if mode.ID == self.ID or mode.solver_number != self.solver_number:
                continue

            ax.add_line(
                x=self.itr_list,
                y=self.get_values(mode),
                label=f'{self.stylized_label} - {mode.stylized_label}'
            )

            ax.set_style(self.plot_style)

    def plot(self, other_supermode=None, row: int = 0, col: int = 0) -> None:
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
    def mode_number(self):
        return self.parent_supermode.mode_number

    @property
    def solver_number(self):
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

    def _interpret_itr_slice_list_(self, *args, **kwargs):
        return self.parent_supermode.parent_set._interpret_itr_slice_list_(*args, **kwargs)

    def _get_symmetrize_vector(self, *args, **kwargs):
        return self.parent_supermode._get_symmetrize_vector(*args, **kwargs)

    def _get_axis_vector(self, *args, **kwargs):
        return self.parent_supermode._get_axis_vector(*args, **kwargs)

    def get_axis(self, *args, **kwargs):
        return self.parent_supermode.get_axis(*args, **kwargs)


class NameSpace():
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


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

        field = self._data[slice_number]

        if add_symmetries:
            field = self._get_symmetrized_field(field=field)

        match normalization.lower():
            case 'max':
                norm = abs(field).max()
            case 'center':
                idx_x_center = numpy.argmin(abs(self.parent_supermode.parent_set.coordinate_system.x_vector))
                idx_y_center = numpy.argmin(abs(self.parent_supermode.parent_set.coordinate_system.y_vector))
                center_value = field[idx_x_center, idx_y_center]
                norm = center_value
            case 'l2':
                norm = numpy.sqrt(numpy.trapz(numpy.trapz(numpy.square(field), dx=1, axis=0), dx=1, axis=0))
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
            case 'symmetric': field = numpy.c_[symmetric_field, field]

            case 'anti-symmetric': field = numpy.c_[-symmetric_field, field]

        match self.boundaries.right.lower():
            case 'symmetric': field = numpy.c_[field, symmetric_field]

            case 'anti-symmetric': field = numpy.c_[field, -symmetric_field]

        symmetric_field = field[-2::-1, :]
        match self.boundaries.top.lower():
            case 'symmetric': field = numpy.r_[field, symmetric_field]

            case 'anti-symmetric': field = numpy.r_[field, -symmetric_field]

        match self.boundaries.bottom.lower():
            case 'symmetric': field = numpy.r_[symmetric_field, field]

            case 'anti-symmetric': field = numpy.r_[-symmetric_field, field]

        return field

    def _render_on_ax_(self, ax: Axis, slice: int) -> None:
        x, y, field = self._get_symmetrized_field_and_axis(field=self._data[slice])

        ax.add_mesh(
            x=x,
            y=y,
            scalar=field,
            colormap=colormaps.blue_black_red
        )

        ax.add_colorbar(symmetric=True, position='right')

        ax.set_style(plot_style.field)

    def plot(self, slice_list: list = [], itr_list: list = []) -> SceneList:
        """
        Plotting method for the fields.

        :param      slice_list:  Value reprenting the slice where the mode field is evaluated.
        :type       slice_list:  list
        :param      itr_list:    Value of itr value to evaluate the mode field.
        :type       itr_list:    list

        :returns:   the figure containing all the plots.
        :rtype:     SceneList
        """
        figure = SceneList(unit_size=(3, 3), tight_layout=True)

        slice_list, itr_list = self._interpret_itr_slice_list_(
            slice_list=slice_list,
            itr_list=itr_list
        )

        slice_list = numpy.atleast_1d(slice_list)
        itr_list = numpy.atleast_1d(itr_list)

        for n, (slice, itr) in enumerate(zip(slice_list, itr_list)):
            ax = figure.append_ax(
                title=f'{self.parent_supermode.stylized_label}\n[slice: {slice}  ITR: {itr:.4f}]'
            )

            self._render_on_ax_(ax=ax, slice=slice)

        return figure


class EigenValue(InheritFromSuperMode, BaseSingleModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self._data = self.parent_supermode.binded_supermode.get_eigen_value()
        self.plot_style = plot_style.eigen_value

    def get_values(self):
        return self._data


class Index(InheritFromSuperMode, BaseSingleModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self._data = self.parent_supermode.binded_supermode.get_index()
        self.plot_style = plot_style.index

    def get_values(self):
        return self._data


class Beta(InheritFromSuperMode, BaseSingleModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self._data = self.parent_supermode.binded_supermode.get_betas()
        self.plot_style = plot_style.beta

    def get_values(self) -> numpy.ndarray:
        return self._data


class Overlap(InheritFromSuperMode, BaseMultiModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self.plot_style = plot_style.overlap

    def get_values(self, other_supermode) -> numpy.ndarray:
        """
        Return the array of the modal coupling for the mode
        """
        return self.parent_supermode.binded_supermode.get_overlap_integrals_with_mode(other_supermode.binded_supermode)


class NormalizedCoupling(InheritFromSuperMode, BaseMultiModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self.plot_style = plot_style.normalized_coupling

    def get_values(self, other_supermode) -> numpy.ndarray:
        """
        Return the array of the modal coupling for the mode
        """

        output = self.parent_supermode.binded_supermode.get_normalized_coupling_with_mode(other_supermode.binded_supermode)

        if not self.parent_supermode.is_computation_compatible(other_supermode):
            output *= 0

        return output


class BeatingLength(InheritFromSuperMode, BaseMultiModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self.plot_style = plot_style.beating_length

    def get_values(self, other_supermode) -> numpy.ndarray:
        """
        Return the array of the modal coupling for the mode
        """
        return self.parent_supermode.binded_supermode.get_beating_length_with_mode(other_supermode.binded_supermode)


class Adiabatic(InheritFromSuperMode, BaseMultiModePlot):
    def __init__(self, parent_supermode):
        self.parent_supermode = parent_supermode
        self.plot_style = plot_style.adiabatic

    def get_values(self, other_supermode) -> numpy.ndarray:
        """
        Return the array of the modal coupling for the mode
        """
        output = self.parent_supermode.binded_supermode.get_adiabatic_with_mode(other_supermode.binded_supermode)

        if not self.parent_supermode.is_computation_compatible(other_supermode):
            output *= numpy.inf

        return output
# -
