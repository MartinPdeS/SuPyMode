#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
import pickle
import numpy
import logging
from dataclasses import dataclass
from pathlib import Path
from itertools import combinations, product
from pathvalidate import sanitize_filepath

# Third-party imports
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pyvista

# Local imports
from SuPyMode.supermode import SuperMode
from SuPyMode.slice_structure import SliceStructure
from SuPyMode.tools import plot_style
from SuPyMode.tools.utils import test_valid_input, get_intersection
from SuPyMode.profiles import AlphaProfile
from SuPyMode.tools import directories
from MPSPlots.render2D import SceneMatrix, SceneList, Axis, Multipage
from MPSPlots import colormaps


@dataclass
class SuperSet(object):
    """
    Solver to which is associated the computed SuperSet Modes

    .. note::
        This class is a representation of the fiber optic structures set of supermodes, hence the name.
        This class has not ling to c++ codes, it is pure Python.
        The items of this class are the supermodes generated from within the SuPySolver

    """
    parent_solver: object
    wavelength: float

    def __post_init__(self):
        self.wavenumber = 2 * numpy.pi / self.wavelength
        self._transmission_matrix = None
        self.supermodes = []
        self._itr_to_slice = interp1d(self.itr_list, numpy.arange(self.itr_list.size))
        self._slice_to_itr = interp1d(numpy.arange(self.itr_list.size), self.itr_list)

    def __getitem__(self, idx: int):
        return self.supermodes[idx]

    def __setitem__(self, idx: int, value):
        self.supermodes[idx] = value

    @property
    def geometry(self):
        """
        Return geometry of the coupler structure
        """
        return self.parent_solver.geometry

    @property
    def itr_list(self):
        """
        Return list of itr value that are used to compute the supermodes
        """
        return self.parent_solver.itr_list

    @property
    def coordinate_system(self):
        """
        Return axes object of the geometry
        """
        return self.parent_solver.geometry.coordinate_system

    @property
    def fundamental_supermodes(self) -> list[SuperMode]:
        """
        Returns a list of the fundamental supermodes.
        Those supermodes are defined as the highest beta-value supermodes associated
        with a particular fiber

        :returns:   The fundamental supermodes
        :rtype:     list[SuperMode]
        """
        return self.get_fundamental_supermodes(tolerance=1e-3)

    @property
    def non_fundamental_supermodes(self) -> list[SuperMode]:
        """
        Returns a list of the non-fundamental supermodes.
        The fundamental supermodes are defined as the highest beta-value supermodes associated
        with a particular fiber

        :returns:   The non-fundamental supermodes
        :rtype:     list[SuperMode]
        """
        return self.get_non_fundamental_supermodes(tolerance=1e-3)

    @property
    def transmission_matrix(self) -> numpy.ndarray:
        """
        Return supermode transfert matrix
        """
        if self._transmission_matrix is None:
            self.compute_transmission_matrix()

        return self._transmission_matrix

    def itr_to_slice(self, itr_list: list[float]) -> list[int]:
        """
        Return slice number associated to itr value

        :param      itr_list:      Inverse taper ration value to evaluate the slice.
        :type       itr_list:      list[float]

        :returns:   List of itr values,
        :rtype:     list[int]
        """
        itr_list = numpy.asarray(itr_list)

        return numpy.floor(self._itr_to_slice(itr_list)).astype(int)

    def slice_to_itr(self, slice_list: list[int]) -> list[float]:
        """
        Return slice number associated to itr value

        :param      slice_list:      Value of the slice to which evaluate the itr.
        :type       slice_list:      list[int]

        :returns:   List of itr values,
        :rtype:     list[float]
        """
        slice_list = numpy.asarray(slice_list) % self.itr_list.size

        return self._slice_to_itr(slice_list)

    def get_fundamental_supermodes(self, *, tolerance: float = 0.1) -> list[SuperMode]:
        """
        Returns list of modes that do not spatially overlap and that have the highest
        propagation constant values.

        :param      tolerance:  The tolerance to which consider the spatial overlap
        :type       tolerance:  float

        :returns:   List of the fundamental modes.
        :rtype:     list
        """
        self.sorting_modes_beta()

        fundamental_supermodes = [self.supermodes[0]]

        def coupling(mode_0: SuperMode, mode_1: SuperMode):
            field_0 = numpy.abs(mode_0.field[0])
            field_1 = numpy.abs(mode_1.field[0])

            return numpy.sum(field_0 * field_1)

        for mode_0 in self.supermodes:
            couplings = [
                coupling(mode_0, mode_1) for mode_1 in fundamental_supermodes
            ]

            couplings = numpy.asarray(couplings)

            if numpy.any(couplings > tolerance):
                continue

            fundamental_supermodes.append(mode_0)

        return fundamental_supermodes

    def get_non_fundamental_supermodes(self, *, tolerance: float = 0.1) -> list[SuperMode]:
        """
        Returns list of modes that do not spatially don't overlap with the fundamental modes.
        Those mode are usually related to higher-order or cladding supermodes.

        :param      tolerance:  The tolerance to which consider the spatial overlap
        :type       tolerance:  float

        :returns:   List of the non-fundamental modes.
        :rtype:     list
        """
        non_fundamental_supermodes = self.supermodes

        for supermodes in self.get_fundamental_supermodes(tolerance=tolerance):
            non_fundamental_supermodes.remove(supermodes)

        return non_fundamental_supermodes

    def get_mode_solver_classification(self) -> list[list[SuperMode]]:
        """
        Returns a list containing the modes ordered per solver number.

        :returns:   The mode solver classification.
        :rtype:     list[list[SuperMode]]
        """
        solver_numbers = [mode.solver_number for mode in self]

        number_of_solvers = len(set(solver_numbers))

        mode_solver_array = [
            [] for i in range(number_of_solvers)
        ]

        for mode in self:
            mode_solver_array[mode.solver_number].append(mode)

        return mode_solver_array

    def label_supermodes(self, *label_list) -> None:
        for n, label in enumerate(label_list):
            self[n].label = label

            setattr(self, label, self[n])

    def reset_labels(self) -> None:
        for n, super_mode in self:
            super_mode.label = f'mode_{n}'

    def swap_supermode_order(self, idx0: int, idx1: int) -> "SuperSet":
        """
        Swap two supermodes.
        it doesn't change any of their characteristic, it only changes the
        order on whihc they will appear, notably for the plots.

        :param      idx0:            Index of the first mode to swap
        :type       idx0:            int
        :param      idx1:            Index of the second mode to swap
        :type       idx1:            int
        """
        self.supermodes[idx0], self.supermodes[idx1] = self.supermodes[idx1], self.supermodes[idx0]

        return self

    def get_slice_structure(self, *, itr: int, add_symmetries: bool = True) -> SliceStructure:
        x, y = self.supermodes[0].get_axis_vector()

        output_slice = SliceStructure(
            parent_superset=self,
            itr=itr,
            supermodes=self.supermodes,
            add_symmetries=add_symmetries
        )

        return output_slice

    def compute_transmission_matrix(self) -> None:
        """
        Calculates the transmission matrix with only the propagation constant included.

        :returns:   The transmission matrix.
        :rtype:     numpy.ndarray
        """
        shape = [
            len(self.supermodes),
            len(self.supermodes),
            len(self.itr_list)
        ]

        self._transmission_matrix = numpy.zeros(shape)

        for mode in self.supermodes:
            self._transmission_matrix[mode.mode_number, mode.mode_number, :] = mode.beta._data * 2.0 * numpy.pi

    def add_coupling_to_t_matrix(self, *, t_matrix: numpy.ndarray, adiabatic_factor: numpy.ndarray) -> numpy.ndarray:
        """
        Add the coupling coefficients to the transmission matrix.

        :param      t_matrix:          The t matrix to which add the coupling values
        :type       t_matrix:          numpy.ndarray
        :param      adiabatic_factor:  The adiabatic factor, if None, it is set to one meaning normalized coupling [z-independent]
        :type       adiabatic_factor:  numpy.ndarray

        :returns:   The transmission matrix.
        :rtype:     numpy.ndarray
        """
        size = t_matrix.shape[-1]

        t_matrix = t_matrix.astype(complex)

        for mode_0, mode_1 in combinations(self.supermodes, 2):

            coupling = mode_0.normalized_coupling.get_values(mode_1)[:size]

            coupling *= adiabatic_factor

            t_matrix[mode_0.mode_number, mode_1.mode_number, :] = - coupling
            t_matrix[mode_1.mode_number, mode_0.mode_number, :] = + coupling

        if numpy.isnan(t_matrix).any() or numpy.isinf(t_matrix).any():
            raise ValueError('Nan or inf values detected in transmission matrix')

        return t_matrix

    def compute_coupling_factor(self, *, coupler_length: float) -> numpy.ndarray:
        r"""
        Compute the coupling factor defined as:
        .. math::
            f_c = \frac{1}{\rho} \frac{d \rho}{d z}

        :param      coupler_length:     The length of the coupler
        :type       coupler_length:     float

        :returns:   The amplitudes as a function of the distance in the coupler
        :rtype:     numpy.ndarray
        """

        dx = coupler_length / (self.itr_list.size)

        ditr = numpy.gradient(numpy.log(self.itr_list), axis=0)

        return ditr / dx

    def get_transmision_matrix_from_profile(self, *, profile: AlphaProfile, add_coupling: bool = True) -> tuple:
        """
        Gets the transmision matrix from profile.

        :param      profile:          The z-profile of the coupler
        :type       profile:          object
        :param      add_coupling:     Add coupling to the transmission matrix
        :type       add_coupling:     bool
        """
        final_slice = self.itr_to_slice(itr_list=profile.smallest_itr)

        sub_t_matrix = self.transmission_matrix[..., :final_slice]

        sub_itr_vector = self.itr_list[: final_slice]

        if add_coupling:
            sub_t_matrix = self.add_coupling_to_t_matrix(
                t_matrix=sub_t_matrix,
                adiabatic_factor=profile.evaluate_adiabatic_factor(itr=sub_itr_vector)
            )

        sub_distance = profile.evaluate_distance_vs_itr(sub_itr_vector)

        return sub_distance, sub_itr_vector, sub_t_matrix

    def propagate(self, *,
            profile: AlphaProfile,
            initial_amplitude: list,
            max_step: float = None,
            n_step: int = None,
            add_coupling: bool = True,
            method: str = 'RK45',
            **kwargs) -> numpy.ndarray:
        """
        Returns the amplitudes value of the supermodes in the coupler.

        :param      initial_amplitude:  The initial amplitude
        :type       initial_amplitude:  list
        :param      profile:            The z-profile of the coupler
        :type       profile:            object
        :param      max_step:           The maximum stride to use in the solver
        :type       max_step:           float
        :param      add_coupling:       Add coupling to the transmission matrix
        :type       add_coupling:       bool
        :param      kwargs:             The keywords arguments to be passed to the solver
        :type       kwargs:             dictionary

        :returns:   The amplitudes as a function of the distance in the coupler
        :rtype:     numpy.ndarray
        """
        initial_amplitude = numpy.asarray(initial_amplitude).astype(complex)

        if max_step is None:
            max_step = self.parent_solver.wavelength / 200

        sub_distance, sub_itr_vector, sub_t_matrix = self.get_transmision_matrix_from_profile(
            profile=profile,
            add_coupling=add_coupling,
        )

        z_to_itr = interp1d(
            profile.distance,
            profile.itr_list,
            bounds_error=False,
            fill_value='extrapolate',
            axis=-1
        )

        itr_to_t_matrix = interp1d(
            sub_itr_vector,
            sub_t_matrix,
            bounds_error=False,
            fill_value='extrapolate',
            axis=-1
        )

        def model(z, y):
            itr = z_to_itr(z)
            return 1j * itr_to_t_matrix(itr).dot(y)

        sol = solve_ivp(
            model,
            y0=initial_amplitude,
            t_span=[0, profile.length],
            vectorized=True,
            max_step=max_step,
            method=method,
            **kwargs
        )

        norm = (numpy.abs(sol.y)**2).sum(axis=0)

        if not numpy.all(numpy.isclose(norm, 1.0, atol=1e-1)):
            logging.warning(f'Warning Power conservation is not acheived [tol: 1e-1]. You should consider reducing the max step size [{max_step = }]')

        return sol.t, sol.y, z_to_itr(sol.t)

    def interpret_initial_input(self, initial_amplitude: list) -> numpy.ndarray:
        amplitude_size = len(initial_amplitude)
        number_of_supermodes = len(self.supermodes)
        assert len(initial_amplitude) == len(self.supermodes), f'Amplitudes size: {amplitude_size} do not match with the number of supermodes: {number_of_supermodes}'

        if isinstance(initial_amplitude, SuperMode):
            return initial_amplitude.amplitudes
        else:
            return numpy.asarray(initial_amplitude).astype(complex)

    def plot_propagation(self, *,
            profile: AlphaProfile,
            initial_amplitude,
            max_step: float = None,
            add_coupling: bool = True,
            method: str = 'RK45',
            sub_sampling: int = 5,
            show_energy: bool = True,
            show_amplitudes: bool = True,
            **kwargs) -> tuple:

        initial_amplitude = self.interpret_initial_input(
            initial_amplitude=initial_amplitude
        )

        z, amplitudes, itr_list = self.propagate(
            initial_amplitude=initial_amplitude,
            profile=profile,
            add_coupling=add_coupling,
            max_step=max_step,
            method=method
        )

        figure = SceneList(unit_size=(12, 4))

        ax = figure.append_ax(
            line_width=2,
            show_legend=True,
            x_label='Propagation distance z',
            y_label='Inverse taper ratio [ITR]'
        )

        for idx, mode in enumerate(self.supermodes):
            color = f"C{idx}"
            if show_energy:
                ax.add_line(
                    x=z[::sub_sampling],
                    y=abs(amplitudes[idx, ::sub_sampling])**2,
                    label=mode.stylized_label,
                    line_width=2.0,
                    line_style='-',
                    color=color
                )

            if show_amplitudes:
                ax.add_line(
                    x=z[::sub_sampling],
                    y=amplitudes[idx, ::sub_sampling].real,
                    label=mode.stylized_label,
                    line_width=2.0,
                    line_style='--',
                    color=color
                )

        if show_energy:
            total_energy = abs(amplitudes)**2
            total_energy = total_energy.sum(axis=0)
            total_energy = numpy.sqrt(total_energy)

            ax.add_line(
                x=z[::sub_sampling],
                y=total_energy[::sub_sampling],
                label='Total energy',
                line_width=3.0,
                line_style='--',
                color='black'
            )

        return figure, (z, amplitudes, itr_list)

    def generate_propagation_gif(self, *,
            profile: AlphaProfile,
            initial_amplitude,
            max_step: float = None,
            coupling: str = 'normalized',
            method: str = 'RK45',
            sub_sampling: int = 5,
            mutliplicative_factor: float = 1,
            save_directory: str = 'new_figure.gif',
            delta_azimuth: float = 0,
            **kwargs) -> tuple:
        """
        Generates a gif video of the mode propagation.

        :param      initial_amplitude:  The initial amplitude
        :type       initial_amplitude:  list
        :param      coupler_length:     The length of the coupler
        :type       coupler_length:     float
        :param      max_step:           The maximum stride to use in the solver
        :type       max_step:           float
        :param      sub_sampling:       Propagation undersampling factor for the video production
        :type       sub_sampling:       int
        :param      kwargs:             The keywords arguments
        :type       kwargs:             dictionary
        """

        initial_amplitude = self.interpret_initial_input(
            initial_amplitude=initial_amplitude
        )

        z_list, amplitudes_list, itr_list = self.propagate(
            initial_amplitude=initial_amplitude,
            profile=profile,
            coupling=coupling,
            max_step=max_step,
            method=method
        )

        self.generate_propagation_gif_from_values(
            amplitudes_list=amplitudes_list,
            itr_list=itr_list,
            z_list=z_list,
            mutliplicative_factor=mutliplicative_factor,
            save_directory=save_directory,
            delta_azimuth=delta_azimuth,
            sub_sampling=sub_sampling
        )

        return z_list, amplitudes_list, itr_list

    def generate_propagation_gif_from_values(self, *,
            amplitudes_list: numpy.ndarray,
            itr_list: numpy.ndarray,
            z_list: numpy.ndarray,
            sub_sampling: int = 10000,
            mutliplicative_factor: float = -100,
            delta_azimuth: float = 0,
            save_directory: str = 'new_figure.gif',
            colormap: str = 'bwr',
            **kwargs) -> None:
        """
        Generates a gif video of the mode propagation.

        :param      initial_amplitude:  The initial amplitude
        :type       initial_amplitude:  list
        :param      coupler_length:     The length of the coupler
        :type       coupler_length:     float
        :param      max_step:           The maximum stride to use in the solver
        :type       max_step:           float
        :param      sub_sampling:       Propagation undersampling factor for the video production
        :type       sub_sampling:       int
        :param      kwargs:             The keywords arguments
        :type       kwargs:             dictionary
        """
        amplitudes_list = amplitudes_list[:, ::sub_sampling]
        itr_list = itr_list[::sub_sampling]
        z_list = z_list[::sub_sampling]

        structure = self.get_slice_structure(itr=1.0, add_symmetries=True)
        total_field = structure.get_field_combination(amplitudes_list[:, 0], Linf_normalization=True) * mutliplicative_factor

        x, y = numpy.mgrid[0: total_field.shape[0], 0: total_field.shape[1]]
        grid = pyvista.StructuredGrid(x, y, total_field)

        plotter = pyvista.Plotter(notebook=False, off_screen=True)
        plotter.open_gif(save_directory, fps=20)
        plotter.view_isometric()
        # plotter.set_background('black', top='white')

        plotter.add_mesh(
            grid,
            scalars=total_field,
            style='surface',
            show_edges=True,
            edge_color='k',
            colormap=colormap,
            show_scalar_bar=False,
            clim=[-100, 100]
        )

        pts = grid.points.copy()
        azimuth = 0
        for z, amplitudes, itr in zip(z_list, amplitudes_list.T, itr_list):
            print(f'itr: {itr}')
            plotter.camera.elevation = -20
            plotter.camera.azimuth = azimuth
            azimuth += delta_azimuth

            structure = self.get_slice_structure(itr=itr, add_symmetries=True)
            total_field = structure.get_field_combination(amplitudes, Linf_normalization=True) * mutliplicative_factor

            pts[:, -1] = total_field.T.ravel()
            plotter.update_coordinates(pts, render=True)
            plotter.update_scalars(total_field.T.ravel(), render=False)
            plotter.add_title(f'ITR: {itr: .3f}\t  z: {z: .3e}', font='courier', color='w', font_size=20)

            plotter.write_frame()

        plotter.close()

    def _sorting_modes_(self, *ordering_list) -> None:
        """
        Generic mode sorting method

        :param      ordering_parameters:  The ordering list to sort the supermodes
        :type       ordering_parameters:  list
        """
        order = numpy.lexsort(ordering_list)

        supermodes = [self.supermodes[idx] for idx in order]

        for n, supermode in enumerate(supermodes):
            supermode.mode_number = n

        return supermodes

    def sorting_modes_beta(self) -> None:
        """
        Re-order modes to sort them in descending value of propagation constant.
        """
        return self._sorting_modes_([-mode.beta[-1] for mode in self.supermodes])

    def sorting_modes(self, *, sorting_method: str = "beta", keep_only: int = None) -> None:
        """
        Re-order modes according to a sorting method, either "beta" or "symmetry+beta".
        The final mode selection will also be filter to keep only a certain number of modes
        """
        assert sorting_method.lower() in ["beta", "symmetry+beta"], \
            f"Unrecognized sortingmethod: {sorting_method}, accepted values are ['beta', 'symmetry+beta']"

        match sorting_method.lower():
            case "beta":
                supermodes = self.sorting_modes_beta()
            case "symmetry+beta":
                supermodes = self.sorting_modes_solver_beta()

        self.all_supermodes = supermodes

        self.supermodes = supermodes[:keep_only]

    def sorting_modes_solver_beta(self):
        """
        Re-order modes to sort them in with two parameters:
        ascending cpp_solver number and descending value of propagation constant.
        """
        return self._sorting_modes_(
            [-mode.beta[-1] for mode in self.supermodes],
            [mode.solver_number for mode in self.supermodes],
        )

    def _interpret_itr_slice_list_(self, *, slice_list: list[int] | int = [], itr_list: list[float] | float = [], sorting: bool = False) -> tuple:
        """
        Interpret and returns list of slice_number and associated itr values for a certain
        slice_list and itr_list input.

        :param      slice_list:  The slice list
        :type       slice_list:  list
        :param      itr_list:    The itr list
        :type       itr_list:    list

        :returns:   A tuple containing interpreted slice_number and itr_list
        :rtype:     tuple
        """

        if slice_list is not None:
            slice_list = numpy.atleast_1d(slice_list)

        if itr_list is not None:
            itr_list = numpy.atleast_1d(itr_list)

        itr_total_list = numpy.concatenate((itr_list, self.slice_to_itr(slice_list)))

        if sorting:
            itr_total_list = numpy.sort(itr_total_list)

        slice_total_list = self.itr_to_slice(itr_list=itr_total_list)

        if len(slice_total_list) == 0:
            slice_total_list = numpy.array([0, -1])
            itr_total_list = self.slice_to_itr(slice_total_list)
            return slice_total_list, itr_total_list

        if len(itr_total_list) == 1:
            return slice_total_list[0], itr_total_list[0]

        else:
            return slice_total_list, itr_total_list

    @staticmethod
    def single_plot(plot_function):
        def wrapper(self, *args, **kwargs):
            figure = SceneList(unit_size=(16, 6), ax_orientation='vertical')
            figure.append_ax()
            plot_function(self, ax=figure[0], *args, **kwargs)

            return figure

        return wrapper

    @single_plot
    def plot_index(self,
            ax: Axis,
            show_crossings: bool = False,
            mode_of_interest: list[SuperMode] = 'all',
            **artist_kwargs) -> SceneList:
        """
        Plot effective index for each mode as a function of itr

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**plot_style.index)

        self._render_index_vs_itr_on_ax_(
            ax=ax,
            show_crossings=show_crossings,
            mode_of_interest=mode_of_interest,
            **artist_kwargs
        )

    def _render_index_vs_itr_on_ax_(self,
            ax: Axis,
            show_crossings: bool = False,
            mode_of_interest: list[SuperMode] = 'all',
            **artist_kwargs) -> None:
        """
        Render the propagation constant of the modes vs ITR on a specific user given axis.

        :param      ax:                The ax to which add the artists
        :type       ax:                Axis
        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   No returns
        :rtype:     None
        """
        mode_of_interest = self.interpret_mode_of_interest(mode_of_interest)

        for mode in mode_of_interest:
            y = mode.index.get_values()

            ax.add_line(
                x=self.itr_list,
                y=y,
                label=f'{mode.stylized_label}'
            )

        if show_crossings:
            crossings = self.get_index_crossing(mode_of_interest=mode_of_interest)
            for crossing in crossings.values():
                ax.add_scatter(
                    x=crossing['itr'],
                    y=crossing['index'],
                    marker='o',
                    color='black',
                    marker_size=5,
                    label='mode crossing'
                )

    @single_plot
    def plot_beta(self,
            ax: Axis,
            show_crossings: bool = False,
            mode_of_interest: list[SuperMode] = 'all',
            **artist_kwargs) -> SceneList:
        """
        Plot propagation constant for each mode as a function of itr

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**plot_style.beta)

        self._render_beta_vs_itr_on_ax_(
            ax=ax,
            show_crossings=show_crossings,
            mode_of_interest=mode_of_interest,
            **artist_kwargs
        )

    def _render_beta_vs_itr_on_ax_(self,
            ax: Axis,
            show_crossings: bool = False,
            mode_of_interest: list[SuperMode] = 'all',
            **artist_kwargs) -> None:
        """
        Render the propagation constant of the modes vs ITR on a specific user given axis.

        :param      ax:                The ax to which add the artists
        :type       ax:                Axis
        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   No returns
        :rtype:     None
        """
        mode_of_interest = self.interpret_mode_of_interest(mode_of_interest)

        for mode in mode_of_interest:
            y = mode.beta.get_values()

            ax.add_line(
                x=self.itr_list,
                y=y,
                label=f'{mode.stylized_label}',
                **artist_kwargs
            )

        if show_crossings:
            crossings = self.get_beta_crossing()
            for crossing in crossings.values():
                ax.add_scatter(
                    x=crossing['itr'],
                    y=crossing['beta'],
                    marker='o',
                    color='black',
                    marker_size=5
                )

    @single_plot
    def plot_eigen_value(self,
            ax: Axis,
            mode_of_interest: list[SuperMode] = 'all',
            **artist_kwargs) -> SceneList:
        """
        Plot propagation constant for each mode as a function of itr

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**plot_style.eigen_value)

        self._render_eigen_values_vs_itr_on_ax_(
            ax=ax,
            mode_of_interest=mode_of_interest,
            **artist_kwargs
        )

    def _render_eigen_values_vs_itr_on_ax_(self,
            ax: Axis,
            mode_of_interest: list[SuperMode] = 'all',
            **artist_kwargs) -> None:
        """
        Render the eigen values of the modes vs ITR on a specific user given axis.

        :param      ax:                The ax to which add the artists
        :type       ax:                Axis
        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   No returns
        :rtype:     None
        """
        mode_of_interest = self.interpret_mode_of_interest(mode_of_interest)

        for mode in mode_of_interest:
            y = mode.eigen_value.get_values()

            ax.add_line(
                x=self.itr_list,
                y=y,
                label=f'{mode.stylized_label}',
                **artist_kwargs
            )

    @single_plot
    def plot_normalized_coupling(self,
            ax: Axis,
            mode_of_interest: list = 'all',
            mode_selection='pairs',
            **artist_kwargs) -> SceneList:
        """
        Plot coupling value for each mode as a function of itr

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      mode_selection:    The mode selection
        :type       mode_selection:    str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**plot_style.normalized_coupling)

        self._render_normalized_coupling_vs_itr_on_ax_(
            ax=ax,
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection,
            **artist_kwargs
        )

    def _render_normalized_coupling_vs_itr_on_ax_(self,
            ax: Axis,
            mode_of_interest: list[SuperMode] = 'all',
            mode_selection: str = 'pairs',
            **artist_kwargs) -> None:
        """
        Render the beating length vs ITR on a specific user given axis.

        :param      ax:                The ax to which add the artists
        :type       ax:                Axis
        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      mode_selection:    The mode selection
        :type       mode_selection:    str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   No returns
        :rtype:     None
        """
        combination = self.interpret_combinations(
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection
        )

        for mode_0, mode_1 in combination:
            if mode_0.is_computation_compatible(mode_1):
                y = numpy.abs(mode_0.normalized_coupling.get_values(other_supermode=mode_1))

                ax.add_line(
                    x=self.itr_list,
                    y=y,
                    label=f'{mode_0.stylized_label} - {mode_1.stylized_label}',
                    **artist_kwargs
                )

    @single_plot
    def plot_overlap(self,
            ax: Axis,
            mode_of_interest: list = 'all',
            mode_selection='pairs',
            **artist_kwargs) -> SceneList:
        """
        Plot overlap value for each mode as a function of itr

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      mode_selection:    The mode selection
        :type       mode_selection:    str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**plot_style.overlap)

        self._render_overlap_vs_itr_on_ax_(
            ax=ax,
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection,
            **artist_kwargs
        )

    def _render_overlap_vs_itr_on_ax_(self,
            ax: Axis,
            mode_of_interest: list[SuperMode] = 'all',
            mode_selection: str = 'pairs',
            **artist_kwargs) -> None:
        """
        Render the beating length vs ITR on a specific user given axis.

        :param      ax:                The ax to which add the artists
        :type       ax:                Axis
        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      mode_selection:    The mode selection
        :type       mode_selection:    str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   No returns
        :rtype:     None
        """
        combination = self.interpret_combinations(
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection
        )

        for mode_0, mode_1 in combination:
            if mode_0.is_computation_compatible(mode_1):
                y = mode_0.overlap.get_values(other_supermode=mode_1)

                ax.add_line(
                    x=self.itr_list,
                    y=y,
                    label=f'{mode_0.stylized_label} - {mode_1.stylized_label}',
                    **artist_kwargs
                )

    @single_plot
    def plot_beating_length(self,
            ax: Axis,
            mode_of_interest: list = 'all',
            mode_selection='pairs',
            add_profile: list[AlphaProfile] = [],
            core_radius: float = None,
            **artist_kwargs) -> SceneList:
        """
        Plot coupling value for each mode as a function of itr

        :param      mode_of_interest:  List of the mode that are to be considered in the adiabatic criterion plotting.
        :type       mode_of_interest:  list

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**plot_style.beating_length)

        self._render_beating_length_vs_itr_on_ax_(
            ax=ax,
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection,
            add_profile=add_profile,
            core_radius=core_radius,
            **artist_kwargs
        )

        for profile in numpy.atleast_1d(add_profile):
            profile._render_taper_length_scale_vs_itr_on_ax_(ax=ax, core_radius=core_radius)

    def _render_beating_length_vs_itr_on_ax_(self,
            ax: Axis,
            mode_of_interest: list[SuperMode] = 'all',
            mode_selection: str = 'pairs',
            add_profile: list[AlphaProfile] = [],
            **artist_kwargs) -> None:
        """
        Render the beating length vs ITR on a specific user given axis.

        :param      ax:                The ax to which add the artists
        :type       ax:                Axis
        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      mode_selection:    The mode selection
        :type       mode_selection:    str
        :param      add_profile:       The add profile
        :type       add_profile:       Array
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   No returns
        :rtype:     None
        """
        combination = self.interpret_combinations(
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection
        )

        for mode_0, mode_1 in combination:
            if mode_0.is_computation_compatible(mode_1):
                y = mode_0.beating_length.get_values(other_supermode=mode_1)

                ax.add_line(
                    x=self.itr_list,
                    y=y,
                    label=f'{mode_0.stylized_label} - {mode_1.stylized_label}',
                    **artist_kwargs
                )

    @single_plot
    def plot_adiabatic(self,
            ax: Axis,
            mode_of_interest: list[SuperMode] = 'all',
            mode_selection: str = 'pairs',
            add_profile: list[AlphaProfile] = [],
            **artist_kwargs) -> SceneList:
        """
        Plot adiabatic criterion for each mode as a function of itr

        :param      pair_of_interest:  List of the mode that are to be considered in the adiabatic criterion plotting.
        :type       pair_of_interest:  list
        :param      mode_selection:    The type of combination to be plotted, either 'specific/all/pairs'
        :type       mode_selection:    str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**plot_style.adiabatic)

        self._render_adiabatic_vs_itr_on_ax_(
            ax=ax,
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection,
            add_profile=add_profile,
            **artist_kwargs
        )

        for profile in numpy.atleast_1d(add_profile):
            profile._render_adiabatic_factor_vs_itr_on_ax_(ax, line_style='--', color='black')

    def _render_adiabatic_vs_itr_on_ax_(self,
            ax: Axis,
            mode_of_interest: list[SuperMode] = 'all',
            mode_selection: str = 'pairs',
            add_profile: list[AlphaProfile] = [],
            **artist_kwargs) -> None:
        """
        Render the adiabatic criterion vs ITR on a specific user given axis.

        :param      ax:                The ax to which add the artists
        :type       ax:                Axis
        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      mode_selection:    The mode selection
        :type       mode_selection:    str
        :param      add_profile:       The add profile
        :type       add_profile:       Array
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   No returns
        :rtype:     None
        """
        combination = self.interpret_combinations(
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection
        )

        for mode_0, mode_1 in combination:
            if mode_0.is_computation_compatible(mode_1):
                y = mode_0.adiabatic.get_values(other_supermode=mode_1)

                ax.add_line(
                    x=self.itr_list,
                    y=y,
                    label=f'{mode_0.stylized_label} - {mode_1.stylized_label}',
                    **artist_kwargs
                )

        for profile in numpy.atleast_1d(add_profile):
            profile._render_adiabatic_factor_vs_itr_on_ax_(ax, line_style='--', color='black')

    @single_plot
    def plot_normalized_adiabatic(self,
            ax: Axis,
            mode_of_interest: list = 'all',
            mode_selection: str = 'pairs') -> SceneList:
        """
         Plot adiabatic criterion for each mode as a function of itr

        :param      pair_of_interest:  List of the mode that are to be considered in the adiabatic criterion plotting.
        :type       pair_of_interest:  list
        :param      mode_selection:    The type of combination to be plotted, either 'specific/all/pairs'
        :type       mode_selection:    str

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**plot_style.adiabatic)

        combination = self.interpret_combinations(
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection
        )

        for mode_0, mode_1 in combination:
            if mode_0.is_computation_compatible(mode_1):
                n0 = mode_0.index._data
                n1 = mode_1.index._data
                beating_length = mode_0.wavelength / abs(n0 - n1)
                y = mode_0.adiabatic.get_values(other_supermode=mode_1) * beating_length

                ax.add_line(
                    x=self.itr_list,
                    y=y,
                    label=f'{mode_0.stylized_label} - {mode_1.stylized_label}'
                )

    def interpret_combinations(self, mode_of_interest: list, mode_selection: str):
        mode_of_interest = self.interpret_mode_of_interest(mode_of_interest)

        combination = self.interpret_mode_selection(
            mode_of_interest=mode_of_interest,
            mode_selection=mode_selection
        )

        return combination

    def is_compute_compatible(self, pair_of_mode: tuple) -> bool:
        """
        Determines whether the specified pair of mode is compatible for computation.

        :param      pair_of_mode:  The pair of mode
        :type       pair_of_mode:  tuple

        :returns:   True if the specified pair of mode is compute compatible, False otherwise.
        :rtype:     bool
        """
        mode_0, mode_1 = pair_of_mode
        return mode_0.is_computation_compatible(mode_1)

    def remove_duplicate_combination(self, supermodes_list: list) -> list[SuperMode]:
        """
        Removes a duplicate combination in the mode combination list irrespectively of the order.

        :param      supermodes_list:  The supermodes list
        :type       supermodes_list:  list

        :returns:   The reduced supermode list
        :rtype:     list
        """
        output_list = []

        for mode0, mode1 in supermodes_list:
            if (mode0, mode1) not in output_list and (mode1, mode0) not in output_list:
                output_list.append((mode0, mode1))

        return output_list

    def interpret_mode_selection(self, mode_of_interest: list, mode_selection: str) -> set:
        """
        Interpret user input for mode selection and return the combination of mode to consider.

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  list
        :param      mode_selection:    The mode selection
        :type       mode_selection:    str
        """
        test_valid_input(
            variable_name='mode_selection',
            user_input=mode_selection,
            valid_inputs=['pairs', 'specific']
        )

        match mode_selection:
            case 'pairs':
                mode_combinations = product(mode_of_interest, mode_of_interest)
            case 'specific':
                mode_combinations = product(mode_of_interest, self.supermodes)

        mode_combinations = filter(self.is_compute_compatible, mode_combinations)

        mode_combinations = self.remove_duplicate_combination(mode_combinations)

        return set(mode_combinations)

    def interpret_mode_of_interest(self, mode_of_interest: list) -> list:
        if isinstance(mode_of_interest, SuperMode):
            return [mode_of_interest]

        if isinstance(mode_of_interest, list):
            return mode_of_interest

        test_valid_input(
            variable_name='mode_of_interest',
            user_input=mode_of_interest,
            valid_inputs=['all', 'fundamental', 'non-fundamental']
        )

        match mode_of_interest:
            case 'fundamental':
                mode_of_interest = self.fundamental_supermodes
            case 'non-fundamental':
                mode_of_interest = self.non_fundamental_supermodes
            case 'all':
                mode_of_interest = self.supermodes

        return mode_of_interest

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
                    colormap=colormaps.blue_black_red,
                    show_colorbar=False
                )

                artist.colorbar.symmetric = True
                artist.colorbar.position = 'right'

                ax.set_style(**plot_style.field)

        return figure

    def get_plot_mode_field_title(self, supermode: SuperMode, itr: float, slice_number: int, show_mode_label: bool, show_itr: bool, show_slice: bool) -> str:
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
            title += f'{supermode.stylized_label}'

        if show_itr or show_slice:
            title += '\n'

        if show_slice:
            title += f'slice: {slice_number}'

        if show_itr:
            title += f'  itr: {itr:.3f}'

        return title

    def plot(self, plot_type: str, **kwargs) -> SceneList:
        """
        Generic plot function.

        Args:
            type: Plot type ['index', 'beta', 'adiabatic', 'normalized-adiabatic', 'overlap', 'coupling', 'field', 'beating-length']
        """
        test_valid_input(
            variable_name='plot_type',
            user_input=plot_type,
            valid_inputs=['index', 'beta', 'eigen-value', 'adiabatic', 'overlap', 'normalized-adiabatic', 'normalized-coupling', 'field', 'beating-length']
        )

        match plot_type.lower():
            case 'index':
                return self.plot_index(**kwargs)
            case 'beta':
                return self.plot_beta(**kwargs)
            case 'eigen-value':
                return self.plot_eigen_value(**kwargs)
            case 'normalized-coupling':
                return self.plot_normalized_coupling(**kwargs)
            case 'overlap':
                return self.plot_overlap(**kwargs)
            case 'adiabatic':
                return self.plot_adiabatic(**kwargs)
            case 'field':
                return self.plot_field(**kwargs)
            case 'beating-length':
                return self.plot_beating_length(**kwargs)
            case 'normalized-adiabatic':
                return self.plot_normalized_adiabatic(**kwargs)

    def generate_pdf_report(self,
            filename: str = "report",
            directory: str = '.',
            itr_list: list[float] = [],
            slice_list: list[int] = [],
            dpi: int = 200,
            mode_of_interest: list = 'all',
            mode_selection: str = 'specific') -> None:
        """
        Generate a full report of the coupler properties as a .pdf file

        :param      filename:          Name of the Report file to be outputed.
        :type       filename:          str
        :param      itr_list:          List of itr value to evaluate the mode field.
        :type       itr_list:          Array
        :param      slice_list:        List of slice value to evaluate the mode field.
        :type       slice_list:        Array
        :param      dpi:               Pixel density for the image included in the report.
        :type       dpi:               int
        :param      mode_of_interest:  List of the mode that are to be considered in the adiabatic criterion plotting.
        :type       mode_of_interest:  list

        :returns:   { description_of_the_return_value }
        :rtype:     None
        """
        if directory == 'auto':
            directory = directories.reports_path

        filename = Path(directory).joinpath(filename).with_suffix('.pdf')

        logging.info(f"Saving report pdf into: {filename}")

        figures = []
        figures.append(self.geometry.plot()._render_())

        figures.append(self.plot_field(itr_list=itr_list, slice_list=slice_list)._render_())

        figures.append(self.plot_index()._render_())

        figures.append(self.plot_beta()._render_())

        figures.append(self.plot_normalized_coupling(mode_of_interest=mode_of_interest, mode_selection=mode_selection)._render_())

        figures.append(self.plot_adiabatic(mode_of_interest=mode_of_interest, mode_selection=mode_selection)._render_())

        Multipage(filename, figs=figures, dpi=dpi)

        for figure in figures:
            figure.close()

    def save_instance(self, filename: str, directory: str = '.') -> Path:
        """
        Saves the superset instance as a serialized pickle file.

        :param      filename:  The directory where to save the file, 'auto' options means the superset_instance folder
        :type       filename:  str
        :param      filename:  The filename
        :type       filename:  str

        :returns:   The path directory of the saved instance
        :rtype:     Path
        """
        if directory == 'auto':
            directory = directories.instance_directory

        filename = Path(filename).with_suffix('.pickle')

        filename = sanitize_filepath(filename)

        filename = Path(directory).joinpath(filename)

        logging.info(f"Saving pickled superset into: {filename}")

        with open(filename, 'wb') as output_file:
            pickle.dump(self, output_file, pickle.HIGHEST_PROTOCOL)

        return filename

    def get_beta_crossing(self, mode_of_interest: list) -> dict:
        """
        Returns a dictionnay of the beta mode-crossing, meaning points where the modes propagation
        constant do cross.

        :returns:   The mode crossing dictionnary.
        :rtype:     dict
        """
        output_dictionnary = {}

        combination = self.interpret_combinations(
            mode_of_interest=mode_of_interest,
            mode_selection='pairs'
        )

        n_crossing = 0
        for mode_0, mode_1 in combination:
            itr, beta = get_intersection(
                x=self.itr_list,
                y0=mode_0.beta.get_values(),
                y1=mode_1.beta.get_values(),
                average=True
            )

            if len(beta) != 0:
                output_dictionnary[n_crossing] = {
                    'mode0': mode_0,
                    'mode1': mode_1,
                    'itr': itr,
                    'beta': beta
                }
                n_crossing += 1

        return output_dictionnary

    def get_index_crossing(self, mode_of_interest: list) -> dict:
        """
        Returns a dictionnay of the beta mode-crossing, meaning points where the modes propagation
        constant do cross.

        :returns:   The mode crossing dictionnary.
        :rtype:     dict
        """
        output_dictionnary = {}

        combination = self.interpret_combinations(
            mode_of_interest=mode_of_interest,
            mode_selection='pairs'
        )

        n_crossing = 0
        for mode_0, mode_1 in combination:
            itr, index = get_intersection(
                x=self.itr_list,
                y0=mode_0.index.get_values(),
                y1=mode_1.index.get_values(),
                average=True
            )

            if len(index) != 0:
                output_dictionnary[n_crossing] = {
                    'mode0': mode_0,
                    'mode1': mode_1,
                    'itr': itr,
                    'index': index
                }
                n_crossing += 1

        return output_dictionnary


# -
