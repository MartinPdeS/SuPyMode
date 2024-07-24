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
from typing import Optional, List, Callable, NoReturn
from FiberFusing.geometry import Geometry

# Third-party imports
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import pyvista

# Local imports
from SuPyMode.binary.ModelParameters import ModelParameters
from SuPyMode.supermode import SuperMode
from SuPyMode.utils import test_valid_input, get_intersection, interpret_slice_number_and_itr, interpret_mode_of_interest
from SuPyMode.profiles import AlphaProfile
from SuPyMode import directories
from MPSPlots.render2D import Multipage
import matplotlib.pyplot as plt


class Propagation:
    def __init__(self, superset: object, distance: numpy.ndarray, profile: AlphaProfile, amplitudes: numpy.ndarray, itr: numpy.ndarray):
        """
        Args:
            profile (AlphaProfile): The profile to propagate.
            amplitude: The modes amplitudes.
        """
        self.superset = superset
        self.supermodes = superset.supermodes
        self.itr = itr
        self.distance = distance
        self.profile = profile
        self.amplitudes = amplitudes

    def plot(self, sub_sampling: int = 5, show_energy: bool = True, show_amplitudes: bool = True, **kwargs: dict) -> NoReturn:
        """
        Plots the propagation of amplitudes over a given profile, showing energy and amplitude plots.

        Args:
            sub_sampling (int): The factor for sub-sampling data for plotting.
            show_energy (bool): Whether to plot the energy of the modes.
            show_amplitudes (bool): Whether to plot the real part of the amplitudes.
            **kwargs (dict): Additional keyword arguments for solver.

        Returns:
            Tuple[plt.Figure, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
            A tuple containing the matplotlib figure object and a tuple with propagation distances, amplitudes, and inverse taper ratios.
        """
        figure, ax = plt.subplots(1, 1)
        ax.set_xlabel('Propagation distance z')
        ax.set_ylabel('Inverse taper ratio [ITR]')

        for idx, mode in enumerate(self.supermodes):
            color = f"C{idx}"
            x_values = self.distance[::sub_sampling]
            y_energy = numpy.abs(self.amplitudes[idx, ::sub_sampling])**2
            y_amplitude = self.amplitudes[idx, ::sub_sampling].real

            if show_energy:
                ax.plot(x_values, y_energy, label=mode.stylized_label, linewidth=2.0, linestyle='-', color=color)
            if show_amplitudes:
                ax.plot(x_values, y_amplitude, label=mode.stylized_label + ' Amplitude', linewidth=2.0, linestyle='--', color=color)

        if show_energy:
            total_energy = numpy.sqrt(numpy.sum(numpy.abs(self.amplitudes)**2, axis=0))[::sub_sampling]
            ax.plot(x_values, total_energy, label='Total energy', linewidth=3.0, linestyle='--', color='black')

        ax.legend()
        plt.show()

    def generate_gif(
            self, *,
            sub_sampling: int = 5,
            mutliplicative_factor: float = -100,
            delta_azimuth: float = 0,
            save_directory: str = 'new_figure.gif',
            colormap: str = 'bwr',
            **kwargs) -> NoReturn:
        """
        Generates a gif video of the mode propagation.

        Args:
            sub_sampling (int): Propagation undersampling factor for the video production.
            mutliplicative_factor (float): Multiplicative factor for scaling.
            save_directory (str): Directory to save the generated GIF.
            delta_azimuth (float): Azimuthal change per frame.
            **kwargs (dict): Additional keyword arguments for solver.

        Returns:
            Tuple: Propagation distances, amplitudes, and ITR list.
        """
        amplitudes_list = self.amplitudes[:, ::sub_sampling]
        itr_list = self.itr[::sub_sampling]
        z_list = self.distance[::sub_sampling]

        structure = self.superset.get_slice_structure(itr=1.0, add_symmetries=True)
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

            structure = self.superset.get_slice_structure(itr=itr, add_symmetries=True)
            total_field = structure.get_field_combination(amplitudes, Linf_normalization=True) * mutliplicative_factor

            pts[:, -1] = total_field.T.ravel()
            plotter.update_coordinates(pts, render=True)
            plotter.update_scalars(total_field.T.ravel(), render=False)
            plotter.add_title(f'ITR: {itr: .3f}\t  z: {z: .3e}', font='courier', color='w', font_size=20)

            plotter.write_frame()

        plotter.close()


@dataclass
class SuperSet(object):
    """
    A class representing a set of supermodes calculated for a specific optical fiber configuration.
    It facilitates operations on supermodes like sorting, plotting, and computations related to fiber optics simulations.

    Attributes:
        model_parameters (ModelParameters):
        wavelength (float): The wavelength used in the solver, in meters.
        geometry (object):
    """
    model_parameters: ModelParameters
    wavelength: float
    geometry: Geometry

    def __post_init__(self):
        self._transmission_matrix = None
        self.supermodes = []
        self._itr_to_slice = interp1d(self.model_parameters.itr_list, numpy.arange(self.model_parameters.n_slice))

    def __getitem__(self, idx: int) -> SuperMode:
        return self.supermodes[idx]

    def __setitem__(self, idx: int, value: SuperMode) -> None:
        self.supermodes[idx] = value

    @property
    def coordinate_system(self):
        """
        Return axes object of the geometry
        """
        return self.geometry.coordinate_system

    @property
    def fundamental_supermodes(self) -> list[SuperMode]:
        """
        Identifies and returns fundamental supermodes based on the highest beta values and minimal spatial overlap.

        Args:
            tolerance (float): The spatial overlap tolerance for mode distinction.

        Returns:
            list[SuperMode]: A list of fundamental supermodes.
        """
        return self.get_fundamental_supermodes(tolerance=1e-2)

    @property
    def non_fundamental_supermodes(self) -> list[SuperMode]:
        """
        Identifies and returns non-fundamental supermodes based on the specified spatial overlap tolerance.

        Args:
            tolerance (float): The spatial overlap tolerance for distinguishing between fundamental and other modes.

        Returns:
            list[SuperMode]: A list of non-fundamental supermodes.
        """
        return self.get_non_fundamental_supermodes(tolerance=1e-2)

    @property
    def transmission_matrix(self) -> numpy.ndarray:
        """
        Return the supermode transmission matrix.
        """
        if self._transmission_matrix is None:
            self.compute_transmission_matrix()

        return self._transmission_matrix

    def itr_to_slice(self, itr_list: list[float]) -> list[int]:
        """
        Convert ITR values to corresponding slice numbers.

        Args:
            itr_list (list[float]): Inverse taper ratio values.

        Returns:
            list[int]: List of slice numbers corresponding to the ITR values.
        """
        itr_list = numpy.asarray(itr_list)

        return numpy.floor(self._itr_to_slice(itr_list)).astype(int)

    def get_fundamental_supermodes(self, *, tolerance: float = 0.1) -> list[SuperMode]:
        """
        Returns a list of fundamental supermodes with the highest propagation constant values and minimal spatial overlap.

        Args:
            tolerance (float): Tolerance for spatial overlap.

        Returns:
            list[SuperMode]: List of fundamental supermodes.
        """
        self.sort_modes_by_beta()

        fundamental_supermodes = [self.supermodes[0]]

        def absolute_overlap(mode_0: SuperMode, mode_1: SuperMode) -> float:
            field_0 = numpy.abs(mode_0.field.data[0])
            norm_0 = field_0.sum()
            field_0 /= numpy.sqrt(norm_0)

            field_1 = numpy.abs(mode_1.field.data[0])
            norm_1 = field_1.sum()
            field_1 /= numpy.sqrt(norm_1)

            overlap = numpy.sum(field_0 * field_1)

            return overlap

        for mode_0 in self.supermodes:
            abs_overlap = [
                absolute_overlap(mode_0, mode_1) for mode_1 in fundamental_supermodes
            ]

            abs_overlaps = numpy.asarray(abs_overlap)

            if numpy.any(abs_overlaps > tolerance):
                continue

            fundamental_supermodes.append(mode_0)

        return fundamental_supermodes

    def get_non_fundamental_supermodes(self, *, tolerance: float = 0.1) -> list[SuperMode]:
        """
        Returns a list of non-fundamental supermodes that do not overlap with the fundamental modes.

        Args:
            tolerance (float): Tolerance for spatial overlap.

        Returns:
            list[SuperMode]: List of non-fundamental supermodes.
        """
        non_fundamental_supermodes = self.supermodes

        for supermodes in self.get_fundamental_supermodes(tolerance=tolerance):
            non_fundamental_supermodes.remove(supermodes)

        return non_fundamental_supermodes

    def get_mode_solver_classification(self) -> list[list[SuperMode]]:
        """
        Returns a list of modes classified by solver number.

        Returns:
            list[list[SuperMode]]: List of lists containing modes classified by solver number.
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
        """
        Assigns labels to the supermodes.

        Args:
            label_list (tuple): Labels to assign to the supermodes.
        """
        for n, label in enumerate(label_list):
            self[n].label = label

            setattr(self, label, self[n])

    def reset_labels(self) -> None:
        """
        Resets labels for all supermodes to default values.
        """
        for n, super_mode in self:
            super_mode.label = f'mode_{n}'

    def compute_transmission_matrix(self) -> None:
        """
        Calculates the transmission matrix with only the propagation constant included.
        """
        shape = [
            len(self.supermodes),
            len(self.supermodes),
            len(self.model_parameters.itr_list)
        ]

        self._transmission_matrix = numpy.zeros(shape)

        for mode in self.supermodes:
            self._transmission_matrix[mode.mode_number, mode.mode_number, :] = mode.beta.data * 2.0 * numpy.pi

    def add_coupling_to_t_matrix(self, *, t_matrix: numpy.ndarray, adiabatic_factor: numpy.ndarray) -> numpy.ndarray:
        """
        Add the coupling coefficients to the transmission matrix.

        Args:
            t_matrix (np.ndarray): Transmission matrix to which coupling values are added.
            adiabatic_factor (np.ndarray): Adiabatic factor, set to one if None (normalized coupling).

        Returns:
            np.ndarray: Updated transmission matrix with coupling values.
        """
        size = t_matrix.shape[-1]

        t_matrix = t_matrix.astype(complex)

        for mode_0, mode_1 in combinations(self.supermodes, 2):

            coupling = mode_0.normalized_coupling.get_values(mode_1)[:size]

            coupling *= adiabatic_factor

            t_matrix[mode_0.mode_number, mode_1.mode_number, :] = - coupling
            t_matrix[mode_1.mode_number, mode_0.mode_number, :] = + coupling

        if numpy.isnan(t_matrix).any():
            raise ValueError('Nan values detected in transmission matrix.')
        if numpy.isinf(t_matrix).any():
            raise ValueError('Inf values detected in transmission matrix, verify that there is no hybrid mode in the computation.')

        return t_matrix

    def compute_coupling_factor(self, *, coupler_length: float) -> numpy.ndarray:
        """
        Compute the coupling factor defined as:

        .. math::
            f_c = \frac{1}{\rho} \frac{d \rho}{d z}

        Args:
            coupler_length (float): Length of the coupler.

        Returns:
            np.ndarray: Coupling factor as a function of distance in the coupler.
        """

        dx = coupler_length / (self.model_parameters.n_slice)

        ditr = numpy.gradient(numpy.log(self.model_parameters.itr_list), axis=0)

        return ditr / dx

    def get_transmision_matrix_from_profile(self, *, profile: AlphaProfile, add_coupling: bool = True) -> tuple:
        """
        Get the transmission matrix from the profile.

        Args:
            profile (AlphaProfile): Z-profile of the coupler.
            add_coupling (bool): Add coupling to the transmission matrix. Defaults to True.

        Returns:
            tuple: Distance, ITR vector, and transmission matrix.
        """
        profile.initialize()

        final_slice = self.itr_to_slice(itr_list=profile.smallest_itr)

        sub_t_matrix = self.transmission_matrix[..., :final_slice]

        sub_itr_vector = self.model_parameters.itr_list[: final_slice]

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
            max_step: Optional[float] = None,
            n_step: Optional[int] = None,
            add_coupling: bool = True,
            method: str = 'RK45',
            **kwargs: dict) -> Propagation:
        """
        Propagates the amplitudes of the supermodes in a coupler based on a given profile.

        Args:
            profile (AlphaProfile): The z-profile of the coupler.
            initial_amplitude (list): The initial amplitude as a list.
            max_step (float, optional): The maximum step size used by the solver. Defaults to None.
            n_step (int, optional): Number of steps used by the solver (not currently used in this method).
            add_coupling (bool): Flag to add coupling to the transmission matrix. Defaults to True.
            method (str): Integration method to be used by the solver. Defaults to 'RK45'.
            **kwargs (Dict[str, Any]): Additional keyword arguments to be passed to the solver.

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]: A tuple containing the times of the solution,
                                                       the solution array of amplitudes, and the interpolated
                                                       index of refraction at those times.
        """
        initial_amplitude = numpy.asarray(initial_amplitude, dtype=complex)

        if max_step is None:
            max_step = self.model_parameters.wavelength / 200

        sub_distance, sub_itr_vector, sub_t_matrix = self.get_transmision_matrix_from_profile(
            profile=profile,
            add_coupling=add_coupling
        )

        z_to_itr = interp1d(profile.distance, profile.itr_list, bounds_error=False, fill_value='extrapolate')
        itr_to_t_matrix = interp1d(sub_itr_vector, sub_t_matrix, bounds_error=False, fill_value='extrapolate')

        def model(z, y):
            itr = z_to_itr(z)
            return 1j * itr_to_t_matrix(itr) @ y

        sol = solve_ivp(
            fun=model,
            y0=initial_amplitude,
            t_span=[0, profile.total_length],
            method=method,
            vectorized=True,
            max_step=max_step,
            **kwargs
        )

        # Check power conservation across the propagation
        norm = numpy.sum(numpy.abs(sol.y)**2, axis=0)
        if not numpy.allclose(norm, 1.0, atol=1e-1):
            logging.warning(f'Power conservation not achieved [{max_step = }, atol = 1e-1].')

        propagation = Propagation(superset=self, amplitudes=sol.y, distance=sol.t, profile=profile, itr=z_to_itr(sol.t))
        return propagation
        # return sol.t, sol.y, z_to_itr(sol.t)

    def interpret_initial_input(self, initial_amplitude: list | SuperMode) -> numpy.ndarray:
        """
        Interprets the initial amplitude input, ensuring compatibility with the expected number of supermodes.

        Args:
            initial_amplitude (list | SuperMode): The initial amplitude as either a list of complex numbers or a SuperMode object.

        Returns:
            np.ndarray: The initial amplitudes as a NumPy array of complex numbers.

        Raises:
            ValueError: If the length of the initial amplitude list does not match the number of supermodes.
        """
        if isinstance(initial_amplitude, SuperMode):
            amplitudes = initial_amplitude.amplitudes
        else:
            amplitudes = initial_amplitude

        amplitude_size = len(amplitudes)
        number_of_supermodes = len(self.supermodes)

        if amplitude_size != number_of_supermodes:
            raise ValueError(f'Amplitudes size: {amplitude_size} does not match with the number of supermodes: {number_of_supermodes}')

        return numpy.asarray(amplitudes, dtype=complex)

    def _sort_modes(self, *ordering_keys) -> List[SuperMode]:
        """
        Sorts supermodes using specified keys provided as tuples in ordering_keys.

        Args:
            ordering_keys (tuple): Tuple containing keys to sort by.

        Returns:
            List[SuperMode]: Sorted list of supermodes.
        """
        order = numpy.lexsort(ordering_keys)
        sorted_supermodes = [self.supermodes[idx] for idx in order]
        for i, supermode in enumerate(sorted_supermodes):
            supermode.mode_number = i
        return sorted_supermodes

    def sort_modes_by_beta(self) -> None:
        """
        Sorts supermodes in descending order of their propagation constants (beta).
        """
        self.all_supermodes = self._sort_modes([-mode.beta.data[-1] for mode in self.supermodes])

    def sort_modes(self, sorting_method: str = "beta", keep_only: Optional[int] = None) -> None:
        """
        Sorts supermodes according to the specified method, optionally limiting the number of modes retained.

        Args:
            sorting_method (str): Sorting method to use, either "beta" or "symmetry+beta".
            keep_only (int, optional): Number of supermodes to retain after sorting.

        Raises:
            ValueError: If an unrecognized sorting method is provided.
        """
        match sorting_method.lower():
            case 'beta':
                self.sort_modes_by_beta()
            case 'symmetry+beta':
                self.sort_modes_by_solver_and_beta()
            case _:
                raise ValueError(f"Unrecognized sorting method: {sorting_method}, accepted values are ['beta', 'symmetry+beta']")

        self.supermodes = self.all_supermodes[:keep_only] if keep_only is not None else self.all_supermodes

    def sort_modes_by_solver_and_beta(self) -> None:
        """
        Sorts supermodes primarily by solver number and secondarily by descending propagation constant (beta).
        """
        self.all_supermodes = self._sort_modes(
            ([mode.solver_number for mode in self.supermodes],
             [-mode.beta[-1] for mode in self.supermodes])
        )

    @staticmethod
    def combination_plot(plot_function) -> Callable:
        def wrapper(self, *args, mode_of_interest='all', mode_selection: str = 'pairs', **kwargs):
            mode_of_interest = interpret_mode_of_interest(
                superset=self,
                mode_of_interest=mode_of_interest
            )

            combination = self.interpret_mode_selection(
                mode_of_interest=mode_of_interest,
                mode_selection=mode_selection
            )

            figure, ax = plt.subplots(1, 1)

            plot_function(self, ax=ax, *args, mode_of_interest=mode_of_interest, combination=combination, **kwargs)

            plt.show()

        return wrapper

    @staticmethod
    def single_plot(plot_function) -> Callable:
        def wrapper(self, *args, mode_of_interest='all', **kwargs):
            mode_of_interest = interpret_mode_of_interest(
                superset=self,
                mode_of_interest=mode_of_interest
            )

            figure, ax = plt.subplots(1, 1)

            plot_function(self, ax=ax, *args, mode_of_interest=mode_of_interest, **kwargs)

            ax.legend()
            plt.show()

        return wrapper

    @single_plot
    def plot_index(
            self,
            ax: plt.Axes,
            show_crossings: bool = False,
            mode_of_interest: str | list[SuperMode] = 'all') -> NoReturn:
        """
        Plot effective index for each mode as a function of itr.

        Args:
            ax (plt.Axes): The matplotlib axis to render the plot on.
            show_crossings (bool): Whether to show crossings in the plot.
            mode_of_interest (str | list[SuperMode]): The mode of interest.

        Returns:
            None
        """
        for mode in mode_of_interest:
            mode.index.render_on_ax(ax=ax)
            mode.index._dress_ax(ax=ax)

        if show_crossings:
            self.add_crossings_to_ax(ax=ax, mode_of_interest=mode_of_interest, data_type='index')

    @single_plot
    def plot_beta(
            self,
            ax: plt.Axes,
            show_crossings: bool = False,
            mode_of_interest: str | list[SuperMode] = 'all') -> NoReturn:
        """
        Plot propagation constant for each mode as a function of itr.

        Args:
            ax (plt.Axes): The matplotlib axis to render the plot on.
            show_crossings (bool): Whether to show crossings in the plot.
            mode_of_interest (str | list[SuperMode]): The mode of interest.

        Returns:
            None
        """
        for mode in mode_of_interest:
            mode.beta.render_on_ax(ax=ax)
            mode.beta._dress_ax(ax=ax)

        if show_crossings:
            self.add_crossings_to_ax(ax=ax, mode_of_interest=mode_of_interest, data_type='beta')

    @single_plot
    def plot_eigen_value(
            self,
            ax: plt.Axes,
            mode_of_interest: str | list[SuperMode] = 'all',
            show_crossings: bool = False) -> NoReturn:
        """
        Plot propagation constant for each mode as a function of itr.

        Args:
            ax (plt.Axes): The matplotlib axis to render the plot on.
            mode_of_interest (str | list[SuperMode]): The mode of interest.
            show_crossings (bool): Whether to show crossings in the plot.

        Returns:
            None
        """
        for mode in mode_of_interest:
            mode.eigen_value.render_on_ax(ax=ax)
            mode.eigen_value._dress_ax(ax=ax)

        if show_crossings:
            self.add_crossings_to_ax(ax=ax, mode_of_interest=mode_of_interest, data_type='eigen_value')

    @combination_plot
    def plot_normalized_coupling(
            self,
            ax: plt.Axes,
            mode_of_interest: list[SuperMode],
            combination: list) -> NoReturn:
        """
        Plot normalized coupling value for each mode as a function of itr.

        Args:
            ax (plt.Axes): The matplotlib axis to render the plot on.
            mode_of_interest (list[SuperMode]): The mode of interest.
            combination (list): The mode combinations.

        Returns:
            None
        """
        for mode_0, mode_1 in combination:
            mode_0.normalized_coupling.render_on_ax(ax=ax, other_supermode=mode_1)
            mode_0.normalized_coupling._dress_ax(ax=ax)

    @combination_plot
    def plot_beating_length(
            self,
            ax: plt.Axes,
            mode_of_interest: list[SuperMode],
            combination: list) -> NoReturn:
        """
        Plot coupling value for each mode as a function of itr.

        Args:
            ax (plt.Axes): The matplotlib axis to render the plot on.
            mode_of_interest (list[SuperMode]): The mode of interest.
            combination (list): The mode combinations.

        Returns:
            None
        """
        for mode_0, mode_1 in combination:
            mode_0.beating_length.render_on_ax(ax=ax, other_supermode=mode_1)
            mode_0.beating_length._dress_ax(ax=ax)

    @combination_plot
    def plot_adiabatic(
            self,
            ax: plt.Axes,
            mode_of_interest: list[SuperMode],
            combination: list,
            add_profile: list[AlphaProfile] = []) -> NoReturn:
        """
        Plot adiabatic criterion for each mode as a function of itr.

        Args:
            ax (plt.Axes): The matplotlib axis to render the plot on.
            mode_of_interest (list[SuperMode]): The mode of interest.
            combination (list): The mode combinations.
            add_profile (list[AlphaProfile]): List of profiles to add to the plot.

        Returns:
            None
        """
        for mode_0, mode_1 in combination:
            mode_0.adiabatic.render_on_ax(ax=ax, other_supermode=mode_1)
            mode_0.adiabatic._dress_ax(ax=ax)

        for profile in numpy.atleast_1d(add_profile):
            profile.render_adiabatic_factor_vs_itr_on_ax(ax=ax, line_style='--')

    def is_compute_compatible(self, pair_of_mode: tuple) -> bool:
        """
        Determines whether the specified pair of mode is compatible for computation.

        Args:
            pair_of_mode (tuple): The pair of modes.

        Returns:
            bool: True if the pair of modes is compute compatible, False otherwise.
        """
        mode_0, mode_1 = pair_of_mode
        return mode_0.is_computation_compatible(mode_1)

    def remove_duplicate_combination(self, supermodes_list: list) -> list[SuperMode]:
        """
        Removes duplicate combinations in the mode combination list irrespective of the order.

        Args:
            supermodes_list (list): List of mode combinations.

        Returns:
            list: Reduced list of unique supermode combinations.
        """
        output_list = []

        for mode0, mode1 in supermodes_list:
            if (mode0, mode1) not in output_list and (mode1, mode0) not in output_list:
                output_list.append((mode0, mode1))

        return output_list

    def interpret_mode_selection(self, mode_of_interest: list, mode_selection: str) -> set:
        """
        Interpret user input for mode selection and return the combination of modes to consider.

        Args:
            mode_of_interest (list): List of modes of interest.
            mode_selection (str): Mode selection method.

        Returns:
            set: Set of mode combinations.
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
            itr_list: list[float] = None,
            slice_list: list[int] = None,
            show_mode_label: bool = True,
            show_itr: bool = True,
            show_slice: bool = True) -> plt.figure:
        """
        Plot each mode field for different ITR values or slice numbers.

        Args:
            mode_of_interest (list): List of modes of interest.
            itr_list (list): List of ITR values to evaluate the mode field.
            slice_list (list): List of slice numbers to evaluate the mode field.
            show_mode_label (bool): Whether to show mode label.
            show_itr (bool): Whether to show ITR.
            show_slice (bool): Whether to show slice.

        Returns:
            plt.Figure: The figure containing the plots.
        """
        slice_list, itr_list = interpret_slice_number_and_itr(
            itr_baseline=self.model_parameters.itr_list,
            itr_list=itr_list,
            slice_list=slice_list
        )

        mode_of_interest = interpret_mode_of_interest(
            superset=self,
            mode_of_interest=mode_of_interest
        )

        unit_size = numpy.array([len(slice_list), len(mode_of_interest)])
        figure, axes = plt.subplots(*unit_size, figsize=3 * numpy.flip(unit_size))

        for m, mode in enumerate(mode_of_interest):
            for n, slice_number in enumerate(slice_list):
                mode.field.render_on_ax(
                    ax=axes[n, m],
                    slice_number=slice_number,
                    show_mode_label=show_mode_label,
                    show_itr=show_itr,
                    show_slice=show_slice
                )

        figure.tight_layout()
        plt.show()

        return figure

    def plot(self, plot_type: str, **kwargs) -> NoReturn:
        """
        General plotting function to handle different types of supermode plots.

        Args:
            plot_type (str): The type of plot to generate. Options include 'index', 'beta', 'eigen-value', etc.
            **kwargs: Additional keyword arguments for specific plot configurations.

        Raises:
            ValueError: If an unrecognized plot type is specified.
        """
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
            case _:
                raise ValueError(f'Invalid plot type: {plot_type}. Options are: index, beta, eigen-value, adiabatic, normalized-adiabatic, normalized-coupling, field, beating-length')

    def generate_pdf_report(
            self,
            filename: str = "report",
            directory: str = '.',
            itr_list: list[float] | None = None,
            slice_list: list[int] | None = None,
            dpi: int = 200,
            mode_of_interest: list = 'all',
            mode_selection: str = 'specific') -> None:
        """
        Generate a full report of the coupler properties as a .pdf file.

        Args:
            filename (str): Name of the report file to be output.
            directory (str): Directory to save the report.
            itr_list (List[float]): List of ITR values to evaluate the mode field.
            slice_list (List[int]): List of slice values to evaluate the mode field.
            dpi (int): Pixel density for the images included in the report.
            mode_of_interest (List): List of modes to consider in the adiabatic criterion plotting.
            mode_selection (str): Method for selecting mode combinations.

        Returns:
            None
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

    def save_instance(self, filename: str, directory: str = 'auto') -> Path:
        """
        Saves the SuperSet instance as a serialized pickle file.

        Args:
            filename (str): Filename for the serialized instance.
            directory (str): Directory to save the file, 'auto' means the instance_directory.

        Returns:
            Path: The path to the saved instance file.
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

    def add_crossings_to_ax(self, ax: plt.Axes, mode_of_interest: list, data_type: str) -> None:
        combination = self.interpret_mode_selection(
            mode_of_interest=mode_of_interest,
            mode_selection='pairs'
        )

        for mode_0, mode_1 in combination:
            x, y = get_intersection(
                x=self.model_parameters.itr_list,
                y0=getattr(mode_0, data_type).data,
                y1=getattr(mode_1, data_type).data,
                average=True
            )

            if x is not None:
                ax.scatter(x=x, y=y, marker='o', color='black', s=20, label='mode crossing')


# -
