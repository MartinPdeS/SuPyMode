#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
import pickle
import numpy
import logging
from dataclasses import dataclass
from pathlib import Path
from itertools import combinations, product
from typing import Optional, List
from FiberFusing.geometry import Geometry


# Third-party imports
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# Local imports
from SuPyMode.binary.ModelParameters import ModelParameters
from SuPyMode.supermode import SuperMode
from SuPyMode.profiles import AlphaProfile
from SuPyMode.propagation import Propagation
from SuPyMode.superset_plots import SuperSetPlots
from SuPyMode.utils import test_valid_input, parse_mode_of_interest, parse_combination, parse_filename

@dataclass
class SuperSet(SuperSetPlots):
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
        for n, super_mode in enumerate(self.supermodes):
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

        return Propagation(
            superset=self,
            amplitudes=sol.y,
            distance=sol.t,
            profile=profile,
            z_to_itr=z_to_itr
        )

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

    def _sort_modes(self, ordering_keys) -> List[SuperMode]:
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
        lexort_index = ([-mode.beta.data[-1] for mode in self.supermodes], )

        self.all_supermodes = self._sort_modes(lexort_index)

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
        lexort_index = (
            [mode.solver_number for mode in self.supermodes],
            [-mode.beta.data[-1] for mode in self.supermodes]
        )

        self.all_supermodes = self._sort_modes(lexort_index)

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

    def interpret_combination(self, mode_of_interest: list, combination: str) -> set:
        """
        Interpret user input for mode selection and return the combination of modes to consider.

        Args:
            mode_of_interest (list): List of modes of interest.
            mode_selection (str): Mode selection method.

        Returns:
            set: Set of mode combinations.
        """
        test_valid_input(
            variable_name='combination',
            user_input=combination,
            valid_inputs=['pairs', 'specific']
        )

        match combination:
            case 'pairs':
                mode_combinations = product(mode_of_interest, mode_of_interest)
            case 'specific':
                mode_combinations = product(mode_of_interest, self.supermodes)

        mode_combinations = filter(self.is_compute_compatible, mode_combinations)

        mode_combinations = self.remove_duplicate_combination(mode_combinations)

        return set(mode_combinations)

    @parse_filename
    def save_instance(self, filename: str) -> Path:
        """
        Saves the SuperSet instance as a serialized pickle file.

        Args:
            filename (str): Filename for the serialized instance.
            directory (str): Directory to save the file, 'auto' means the instance_directory.

        Returns:
            Path: The path to the saved instance file.
        """
        with open(filename.with_suffix('.pickle'), 'wb') as output_file:
            pickle.dump(self, output_file, pickle.HIGHEST_PROTOCOL)

        return filename

    @parse_mode_of_interest
    @parse_combination
    @parse_filename
    def export_data(self,
                    filename: str,
                    mode_of_interest: list = 'all',
                    combination: list = None,
                    export_index: bool = True,
                    export_beta: bool = True,
                    export_eigen_value: bool = False,
                    export_adiabatic: bool = True,
                    export_beating_length: bool = True,
                    export_normalized_coupling: bool = True) -> Path:
        """
        Export the SuperSet data as CSV files, saving specific attributes of the modes or combinations of modes.

        Args:
            filename (str): The directory where the files will be saved.
            mode_of_interest (list): List of modes to be exported. Defaults to 'all'.
            combination (list): List of mode combinations to be exported. Defaults to None.
            export_index (bool): Whether to export the 'index' attribute. Defaults to True.
            export_beta (bool): Whether to export the 'beta' attribute. Defaults to True.
            export_eigen_value (bool): Whether to export the 'eigen_value' attribute. Defaults to False.
            export_adiabatic (bool): Whether to export the 'adiabatic' attribute for combinations. Defaults to True.
            export_beating_length (bool): Whether to export the 'beating_length' attribute for combinations. Defaults to True.
            export_normalized_coupling (bool): Whether to export the 'normalized_coupling' attribute for combinations. Defaults to True.

        Returns:
            Path: The path to the directory where the files were saved.
        """
        from pathlib import Path
        import numpy as np

        # Create the directory if it doesn't exist
        output_dir = Path(filename)
        output_dir.mkdir(parents=True, exist_ok=True)

        def _export_single_data(attribute_name: str):
            """Helper function to export data for a single attribute for all modes of interest."""
            for mode in mode_of_interest:
                data = getattr(mode, attribute_name).data
                sub_filename = output_dir / f"{attribute_name}_{mode.label}".replace("}", "").replace("{", "")
                data = np.vstack([self.model_parameters.itr_list, data])
                np.savetxt(fname=sub_filename.with_suffix('.csv'), X=data, delimiter=',', header=f'ITR, {attribute_name}')

        def _export_combination_data(attribute_name: str):
            """Helper function to export data for a single attribute for all combinations of modes."""
            for mode_0, mode_1 in combination:
                data = getattr(mode_0, attribute_name).get_values(other_supermode=mode_1)
                sub_filename = output_dir / f"{attribute_name}_{mode_0.label}_{mode_1.label}".replace("}", "").replace("{", "")
                data = np.vstack([self.model_parameters.itr_list, data])
                np.savetxt(fname=sub_filename.with_suffix('.csv'), X=data, delimiter=',', header=f'ITR, {attribute_name}')

        # Export single-mode attributes
        if export_index:
            _export_single_data('index')
        if export_beta:
            _export_single_data('beta')
        if export_eigen_value:
            _export_single_data('eigen_value')

        # Export combination-mode attributes
        if export_adiabatic:
            _export_combination_data('adiabatic')
        if export_beating_length:
            _export_combination_data('beating_length')
        if export_normalized_coupling:
            _export_combination_data('normalized_coupling')



# -
