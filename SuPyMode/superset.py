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
from SuPyMode import representation
from SuPyMode.tools.utils import test_valid_input, get_intersection, interpret_slice_number_and_itr, interpret_mode_of_interest
from SuPyMode.profiles import AlphaProfile
from SuPyMode.tools import directories
from MPSPlots.render2D import SceneMatrix, SceneList, Axis, Multipage


@dataclass
class SuperSet(object):
    """
    Solver to which is associated the computed SuperSet Modes.
    This class is a representation of the fiber optic structures set of supermodes, hence the name.
    The items of this class are the supermodes generated from within the SuPySolver.
    It doesn't link to any c++ binding, it is pure Python.

    """
    parent_solver: object
    wavelength: float

    def __post_init__(self):
        self.wavenumber = 2 * numpy.pi / self.wavelength
        self._transmission_matrix = None
        self.supermodes = []
        self._itr_to_slice = interp1d(self.itr_list, numpy.arange(self.itr_list.size))

    def __getitem__(self, idx: int) -> SuperMode:
        return self.supermodes[idx]

    def __setitem__(self, idx: int, value: SuperMode) -> None:
        self.supermodes[idx] = value

    @property
    def geometry(self):
        """
        Return geometry of the coupler structure
        """
        return self.parent_solver.geometry

    @property
    def itr_list(self) -> numpy.ndarray:
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
        return self.get_fundamental_supermodes(tolerance=1e-2)

    @property
    def non_fundamental_supermodes(self) -> list[SuperMode]:
        """
        Returns a list of the non-fundamental supermodes.
        The fundamental supermodes are defined as the highest beta-value supermodes associated
        with a particular fiber

        :returns:   The non-fundamental supermodes
        :rtype:     list[SuperMode]
        """
        return self.get_non_fundamental_supermodes(tolerance=1e-2)

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
        profile.initialize()

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

    def propagate(
            self, *,
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
        profile.initialize()

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
            t_span=[0, profile.total_length],
            vectorized=True,
            max_step=max_step,
            method=method,
            **kwargs
        )

        norm = (numpy.abs(sol.y)**2).sum(axis=0)

        if not numpy.all(numpy.isclose(norm, 1.0, atol := 1e-1)):
            logging.warning(f'Power conservation is not acheived [{atol = }]. You should consider reducing the max step size [{max_step = }]')

        return sol.t, sol.y, z_to_itr(sol.t)

    def interpret_initial_input(self, initial_amplitude: list) -> numpy.ndarray:
        amplitude_size = len(initial_amplitude)
        number_of_supermodes = len(self.supermodes)
        assert len(initial_amplitude) == len(self.supermodes), f'Amplitudes size: {amplitude_size} do not match with the number of supermodes: {number_of_supermodes}'

        if isinstance(initial_amplitude, SuperMode):
            return initial_amplitude.amplitudes
        else:
            return numpy.asarray(initial_amplitude).astype(complex)

    def plot_propagation(
            self, *,
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

    def generate_propagation_gif(
            self, *,
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

    def generate_propagation_gif_from_values(
            self, *,
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

    def sorting_modes_solver_beta(self) -> list[SuperMode]:
        """
        Re-order modes to sort them in with two parameters:
        ascending cpp_solver number and descending value of propagation constant.

        :returns:   list of supermode in ordered beta
        :rtype:     list[SuperMode]
        """
        return self._sorting_modes_(
            [-mode.beta[-1] for mode in self.supermodes],
            [mode.solver_number for mode in self.supermodes],
        )

    @staticmethod
    def single_plot(plot_function):
        def wrapper(self, *args, mode_of_interest='all', **kwargs):
            mode_of_interest = interpret_mode_of_interest(
                superset=self,
                mode_of_interest=mode_of_interest
            )

            figure = SceneList(unit_size=(16, 6), ax_orientation='vertical')

            ax = figure.append_ax()

            plot_function(self, ax=ax, *args, mode_of_interest=mode_of_interest, **kwargs)

            return figure

        return wrapper

    @staticmethod
    def combination_plot(plot_function):
        def wrapper(self, *args, mode_of_interest='all', mode_selection: str = 'pairs', **kwargs):
            mode_of_interest = interpret_mode_of_interest(
                superset=self,
                mode_of_interest=mode_of_interest
            )

            combination = self.interpret_mode_selection(
                mode_of_interest=mode_of_interest,
                mode_selection=mode_selection
            )

            figure = SceneList(unit_size=(16, 6), ax_orientation='vertical')

            ax = figure.append_ax()

            plot_function(self, ax=ax, *args, mode_of_interest=mode_of_interest, combination=combination, **kwargs)

            return figure

        return wrapper

    @single_plot
    def plot_index(
            self,
            ax: Axis,
            show_crossings: bool = False,
            mode_of_interest: str | list[SuperMode] = 'all') -> SceneList:
        """
        Plot effective index for each mode as a function of itr

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**representation.index.ax_style)

        for mode in mode_of_interest:
            mode.index.render_on_ax(ax=ax)

        if show_crossings:
            self.add_crossings_to_ax(ax=ax, mode_of_interest=mode_of_interest, data_type='index')

    @single_plot
    def plot_beta(
            self,
            ax: Axis,
            show_crossings: bool = False,
            mode_of_interest: str | list[SuperMode] = 'all') -> SceneList:
        """
        Plot propagation constant for each mode as a function of itr

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**representation.beta.ax_style)

        for mode in mode_of_interest:
            mode.beta.render_on_ax(ax=ax)

        if show_crossings:
            self.add_crossings_to_ax(ax=ax, mode_of_interest=mode_of_interest, data_type='beta')

    @single_plot
    def plot_eigen_value(
            self,
            ax: Axis,
            mode_of_interest: str | list[SuperMode] = 'all',
            show_crossings: bool = False) -> SceneList:
        """
        Plot propagation constant for each mode as a function of itr

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**representation.eigen_value.ax_style)

        for mode in mode_of_interest:
            mode.index.render_on_ax(ax=ax)

        if show_crossings:
            self.add_crossings_to_ax(ax=ax, mode_of_interest=mode_of_interest, data_type='eigen_value')

    @combination_plot
    def plot_normalized_coupling(
            self,
            ax: Axis,
            mode_of_interest: list[SuperMode],
            combination: list) -> SceneList:
        """
        Plot normalized coupling value for each mode as a function of itr.

        :param      mode_of_interest:  The mode of interest
        :type       mode_of_interest:  str
        :param      mode_selection:    The mode selection
        :type       mode_selection:    str
        :param      artist_kwargs:     The keywords arguments
        :type       artist_kwargs:     dictionary

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        ax.set_style(**representation.normalized_coupling.ax_style)

        for mode_0, mode_1 in combination:
            mode_0.normalized_coupling.render_on_ax(ax=ax, other_supermode=mode_1)

    @combination_plot
    def plot_beating_length(
            self,
            ax: Axis,
            mode_of_interest: list[SuperMode],
            combination: list) -> SceneList:
        """
        Plot coupling value for each mode as a function of itr

        :param      mode_of_interest:  List of the mode that are to be considered in the adiabatic criterion plotting.
        :type       mode_of_interest:  list

        :returns:   figure instance, to plot the show() method.
        :rtype:     SceneList
        """
        for mode_0, mode_1 in combination:
            ax.set_style(**mode_0.beating_length.ax_style)
            mode_0.beating_length.render_on_ax(ax=ax, other_supermode=mode_1)

    @combination_plot
    def plot_adiabatic(
            self,
            ax: Axis,
            mode_of_interest: list[SuperMode],
            combination: list,
            add_profile: list[AlphaProfile] = []) -> SceneList:
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
        ax.set_style(**representation.adiabatic.ax_style)
        for mode_0, mode_1 in combination:
            mode_0.adiabatic.render_on_ax(ax=ax, other_supermode=mode_1)

        for profile in numpy.atleast_1d(add_profile):
            profile.render_adiabatic_factor_vs_itr_on_ax(ax=ax, line_style='--', line_color='black')

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

    def plot_field(
            self,
            mode_of_interest: list = 'all',
            itr_list: list[float] = [],
            slice_list: list[int] = [0, -1],
            show_mode_label: bool = True,
            show_itr: bool = True,
            show_slice: bool = True) -> SceneList:
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

        mode_of_interest = interpret_mode_of_interest(
            superset=self,
            mode_of_interest=mode_of_interest
        )

        for m, mode in enumerate(mode_of_interest):
            for n, slice_number in enumerate(slice_list):
                ax = figure.append_ax(
                    row=n,
                    column=m,
                )

                ax.set_style(**representation.field.ax_style)

                mode.field.render_on_ax(
                    ax=ax,
                    slice_number=slice_number,
                    show_mode_label=show_mode_label,
                    show_itr=show_itr,
                    show_slice=show_slice
                )

        return figure

    def plot(self, plot_type: str, **kwargs) -> SceneList:
        """
        Generic plot function.

        Args:
            type: Plot type ['index', 'beta', 'adiabatic', 'normalized-adiabatic', 'coupling', 'field', 'beating-length']
        """
        test_valid_input(
            variable_name='plot_type',
            user_input=plot_type,
            valid_inputs=['index', 'beta', 'eigen-value', 'adiabatic', 'normalized-adiabatic', 'normalized-coupling', 'field', 'beating-length']
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

    def generate_pdf_report(
            self,
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

        :returns:   No return
        :rtype:     None
        """
        if directory == 'auto':
            directory = directories.reports_path

        filename = Path(directory).joinpath(filename).with_suffix('.pdf')

        logging.info(f"Saving report pdf into: {filename}")

        figure_list = []

        figure_list.append(self.geometry.plot()._render_())

        figure_list.append(self.plot_field(itr_list=itr_list, slice_list=slice_list)._render_())

        figure_list.append(self.plot_index()._render_())

        figure_list.append(self.plot_beta()._render_())

        figure_list.append(self.plot_normalized_coupling(mode_of_interest=mode_of_interest, mode_selection=mode_selection)._render_())

        figure_list.append(self.plot_adiabatic(mode_of_interest=mode_of_interest, mode_selection=mode_selection)._render_())

        Multipage(filename, figs=figure_list, dpi=dpi)

        for figure in figure_list:
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

    def add_crossings_to_ax(self, ax: Axis, mode_of_interest: list, data_type: str) -> None:
        combination = self.interpret_mode_selection(
            mode_of_interest=mode_of_interest,
            mode_selection='pairs'
        )

        for mode_0, mode_1 in combination:
            x, y = get_intersection(
                x=self.itr_list,
                y0=getattr(mode_0, data_type)._data,
                y1=getattr(mode_1, data_type)._data,
                average=True
            )

            if x is not None:
                ax.add_scatter(
                    x=x,
                    y=y,
                    marker='o',
                    color='black',
                    marker_size=20,
                    label='mode crossing'
                )


# -
