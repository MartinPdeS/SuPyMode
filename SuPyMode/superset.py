#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
import logging
import pickle
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
from typing import List, Optional

import numpy

# Third-party imports
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

# Local imports
from SuPyMode.binary.interface_model_parameters import ModelParameters
from SuPyMode.binary.interface_supermode import SUPERMODE

from SuPyMode.binary.interface_taper import AlphaProfile
from SuPyMode.propagation import Propagation
from SuPyMode.superset_plots import SuperSetPlots
from SuPyMode.utils import parse_filename


@dataclass
class SuperSet(SuperSetPlots):
    """
    A set of supermodes calculated for a specific optical fiber configuration.

    This class manages operations on supermodes, including sorting, plotting, and calculations
    related to fiber optics simulations.

    Parameters
    ----------
    model_parameters : ModelParameters
        Parameters defining the model for simulation.
    """

    model_parameters: ModelParameters

    def __post_init__(self):
        self._transmission_matrix = None
        self.supermodes = []
        itr_list = self.model_parameters.itr_list
        self._itr_to_slice = interp1d(
            itr_list,
            numpy.arange(itr_list.size),
        )

    def __getitem__(self, idx: int) -> SUPERMODE:
        return self.supermodes[idx]

    @property
    def fundamental_supermodes(self) -> list[SUPERMODE]:
        """
        Returns the fundamental supermodes based on the highest beta values and minimal spatial overlap.

        Returns
        -------
        list of SUPERMODE
            A list of fundamental supermodes.
        """
        return self.get_fundamental_supermodes(tolerance=1e-2)

    @property
    def non_fundamental_supermodes(self) -> list[SUPERMODE]:
        """
        Returns the non-fundamental supermodes based on the specified spatial overlap tolerance.

        Returns
        -------
        list of SUPERMODE
            A list of non-fundamental supermodes.
        """
        return self.get_non_fundamental_supermodes(tolerance=1e-2)

    @property
    def transmission_matrix(self) -> numpy.ndarray:
        """
        Returns the supermode transmission matrix.

        Returns
        -------
        numpy.ndarray
            The transmission matrix for the supermodes.
        """
        if self._transmission_matrix is None:
            self.compute_transmission_matrix()

        return self._transmission_matrix

    def itr_to_slice(self, itr_list: list[float]) -> list[int]:
        """
        Convert ITR values to corresponding slice numbers.

        Parameters
        ----------
        itr_list : list of float
            Inverse taper ratio values.

        Returns
        -------
        list of int
            List of slice numbers corresponding to the ITR values.
        """
        itr_list = numpy.asarray(itr_list)

        return numpy.floor(self._itr_to_slice(itr_list)).astype(int)

    def get_fundamental_supermodes(self, *, tolerance: float = 0.1) -> list[SUPERMODE]:
        """
        Returns a list of fundamental supermodes with the highest propagation constant values and minimal spatial overlap.

        Parameters
        ----------
        tolerance : float, optional
            Tolerance for spatial overlap (default is 0.1).

        Returns
        -------
        list of SUPERMODE
            List of fundamental supermodes.
        """
        self.sort_modes_by_beta()

        fundamental_supermodes = [self.supermodes[0]]

        def absolute_overlap(mode_0: SUPERMODE, mode_1: SUPERMODE) -> float:
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

    def get_non_fundamental_supermodes(
        self, *, tolerance: float = 0.1
    ) -> list[SUPERMODE]:
        """
        Returns a list of non-fundamental supermodes that do not overlap with the fundamental modes.

        Parameters
        ----------
        tolerance : float, optional
            Tolerance for spatial overlap (default is 0.1).

        Returns
        -------
        list of SUPERMODE
            List of non-fundamental supermodes.
        """
        non_fundamental_supermodes = self.supermodes

        for supermodes in self.get_fundamental_supermodes(tolerance=tolerance):
            non_fundamental_supermodes.remove(supermodes)

        return non_fundamental_supermodes

    def get_mode_solver_classification(self) -> list[list[SUPERMODE]]:
        """
        Returns a list of modes classified by solver number.

        Returns
        -------
        list of list of SUPERMODE
            List of lists containing modes classified by solver number.
        """
        solver_numbers = [mode.solver_number for mode in self]

        number_of_solvers = len(set(solver_numbers))

        mode_solver_array = [[] for i in range(number_of_solvers)]

        for mode in self:
            mode_solver_array[mode.solver_number].append(mode)

        return mode_solver_array

    def label_supermodes(self, *label_list) -> None:
        """
        Assign labels to the supermodes.

        Parameters
        ----------
        label_list : tuple
            Labels to assign to the supermodes.
        """
        for n, label in enumerate(label_list):
            self[n].label = label

            setattr(self, label, self[n])

    def reset_labels(self) -> None:
        """
        Reset labels for all supermodes to default values.
        """
        for n, super_mode in enumerate(self.supermodes):
            super_mode.label = f"mode_{n}"

    def compute_transmission_matrix(self) -> None:
        """
        Calculate the transmission matrix including only the propagation constant.
        """
        shape = [
            len(self.supermodes),
            len(self.supermodes),
            len(self.model_parameters.itr_list),
        ]

        self._transmission_matrix = numpy.zeros(shape)

        for mode in self.supermodes:
            self._transmission_matrix[mode.mode_number, mode.mode_number, :] = (
                mode.beta.data * 2.0 * numpy.pi
            )

    def add_coupling_to_t_matrix(
        self, *, t_matrix: numpy.ndarray, adiabatic_factor: numpy.ndarray
    ) -> numpy.ndarray:
        """
        Add coupling coefficients to the transmission matrix.

        Parameters
        ----------
        t_matrix : numpy.ndarray
            Transmission matrix to which coupling values are added.
        adiabatic_factor : numpy.ndarray
            Adiabatic factor, set to one if None (normalized coupling).

        Returns
        -------
        numpy.ndarray
            Updated transmission matrix with coupling values.
        """
        size = t_matrix.shape[-1]

        t_matrix = t_matrix.astype(complex)

        for mode_0, mode_1 in combinations(self.supermodes, 2):
            coupling = mode_0.normalized_coupling.get_values(mode_1)[:size]

            coupling *= adiabatic_factor

            t_matrix[mode_0.mode_number, mode_1.mode_number, :] = -coupling
            t_matrix[mode_1.mode_number, mode_0.mode_number, :] = +coupling

        if numpy.isnan(t_matrix).any():
            raise ValueError("Nan values detected in transmission matrix.")
        if numpy.isinf(t_matrix).any():
            raise ValueError(
                "Inf values detected in transmission matrix, verify that there is no hybrid mode in the computation."
            )

        return t_matrix

    def compute_coupling_factor(self, *, coupler_length: float) -> numpy.ndarray:
        """
        Compute the coupling factor defined by the derivative of the inverse taper ratio.

        Parameters
        ----------
        coupler_length : float
            Length of the coupler.

        Returns
        -------
        numpy.ndarray
            Coupling factor as a function of distance in the coupler.
        """
        dx = coupler_length / (self.model_parameters.n_slice)

        ditr = numpy.gradient(numpy.log(self.model_parameters.itr_list), axis=0)

        return ditr / dx

    def get_transmision_matrix_from_profile(
        self, *, profile: AlphaProfile, add_coupling: bool = True
    ) -> tuple:
        """
        Get the transmission matrix from the given profile.

        Parameters
        ----------
        profile : AlphaProfile
            Z-profile of the coupler.
        add_coupling : bool, optional
            Whether to add coupling to the transmission matrix (default is True).

        Returns
        -------
        tuple
            Tuple containing distance, ITR vector, and transmission matrix.
        """
        profile.initialize()

        final_slice = self.itr_to_slice(itr_list=profile.smallest_itr)

        sub_t_matrix = self.transmission_matrix[..., :final_slice]

        sub_itr_vector = self.model_parameters.itr_list[:final_slice]

        if add_coupling:
            sub_t_matrix = self.add_coupling_to_t_matrix(
                t_matrix=sub_t_matrix,
                adiabatic_factor=profile.evaluate_adiabatic_factor(itr=sub_itr_vector),
            )

        sub_distance = profile.evaluate_distance_vs_itr(sub_itr_vector)

        return sub_distance, sub_itr_vector, sub_t_matrix

    def propagate(
        self,
        *,
        profile: AlphaProfile,
        initial_amplitude: list,
        max_step: Optional[float] = None,
        n_step: Optional[int] = None,
        add_coupling: bool = True,
        method: str = "RK45",
        **kwargs: dict,
    ) -> Propagation:
        """
        Propagate the amplitudes of the supermodes in a coupler based on the given profile.

        Parameters
        ----------
        profile : AlphaProfile
            The z-profile of the coupler.
        initial_amplitude : list
            The initial amplitude as a list of complex numbers.
        max_step : float, optional
            The maximum step size used by the solver (default is None).
        n_step : int, optional
            Number of steps used by the solver (default is None).
        add_coupling : bool, optional
            Whether to add coupling to the transmission matrix (default is True).
        method : str, optional
            Integration method to be used by the solver (default is 'RK45').
        **kwargs : dict
            Additional keyword arguments to be passed to the solver.

        Returns
        -------
        Propagation
            A Propagation object containing the results of the propagation.
        """
        initial_amplitude = numpy.asarray(initial_amplitude, dtype=complex)

        if max_step is None:
            max_step = self.model_parameters.wavelength / 200

        sub_distance, sub_itr_vector, sub_t_matrix = (
            self.get_transmision_matrix_from_profile(
                profile=profile, add_coupling=add_coupling
            )
        )

        z_to_itr = interp1d(
            profile.distance,
            profile.itr_list,
            bounds_error=False,
            fill_value="extrapolate",
        )
        itr_to_t_matrix = interp1d(
            sub_itr_vector, sub_t_matrix, bounds_error=False, fill_value="extrapolate"
        )

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
            **kwargs,
        )

        # Check power conservation across the propagation
        norm = numpy.sum(numpy.abs(sol.y) ** 2, axis=0)
        if not numpy.allclose(norm, 1.0, atol=1e-1):
            logging.warning(
                f"Power conservation not achieved [{max_step = }, atol = 1e-1]."
            )

        return Propagation(
            superset=self,
            amplitudes=sol.y,
            distance=sol.t,
            profile=profile,
            z_to_itr=z_to_itr,
        )

    def interpret_initial_input(
        self, initial_amplitude: list | SUPERMODE
    ) -> numpy.ndarray:
        """
        Interpret the initial amplitude input to ensure compatibility with the expected number of supermodes.

        Parameters
        ----------
        initial_amplitude : list or SUPERMODE
            The initial amplitude as either a list of complex numbers or a SUPERMODE object.

        Returns
        -------
        numpy.ndarray
            The initial amplitudes as a NumPy array of complex numbers.

        Raises
        ------
        ValueError
            If the length of the initial amplitude list does not match the number of supermodes.
        """
        if isinstance(initial_amplitude, SUPERMODE):
            amplitudes = initial_amplitude.amplitudes
        else:
            amplitudes = initial_amplitude

        amplitude_size = len(amplitudes)
        number_of_supermodes = len(self.supermodes)

        if amplitude_size != number_of_supermodes:
            raise ValueError(
                f"Amplitudes size: {amplitude_size} does not match with the number of supermodes: {number_of_supermodes}"
            )

        return numpy.asarray(amplitudes, dtype=complex)

    def _sort_modes(self, ordering_keys) -> List[SUPERMODE]:
        """
        Sort supermodes using specified keys provided as tuples in `ordering_keys`.

        Parameters
        ----------
        ordering_keys : tuple
            Tuple containing keys to sort by.

        Returns
        -------
        list of SUPERMODE
            Sorted list of supermodes.
        """
        order = numpy.lexsort(ordering_keys)
        sorted_supermodes = [self.supermodes[idx] for idx in order]
        for i, supermode in enumerate(sorted_supermodes):
            supermode.mode_number = i
        return sorted_supermodes

    def sort_modes_by_beta(self) -> None:
        """
        Sort supermodes in descending order of their propagation constants (beta).
        """
        lexort_index = ([-mode.get_betas()[-1] for mode in self.supermodes],)

        print(lexort_index)

        self.all_supermodes = self._sort_modes(lexort_index)

    def sort_modes(
        self, sorting_method: str = "beta", keep_only: Optional[int] = None
    ) -> None:
        """
        Sort supermodes according to the specified method, optionally limiting the number of modes retained.

        Parameters
        ----------
        sorting_method : str, optional
            Sorting method to use, either "beta" or "symmetry+beta" (default is "beta").
        keep_only : int, optional
            Number of supermodes to retain after sorting (default is None).

        Raises
        ------
        ValueError
            If an unrecognized sorting method is provided.
        """
        match sorting_method.lower():
            case "beta":
                self.sort_modes_by_beta()
            case "symmetry+beta":
                self.sort_modes_by_solver_and_beta()
            case _:
                raise ValueError(
                    f"Unrecognized sorting method: {sorting_method}, accepted values are ['beta', 'symmetry+beta']"
                )

        self.supermodes = (
            self.all_supermodes[:keep_only]
            if keep_only is not None
            else self.all_supermodes
        )

    def sort_modes_by_solver_and_beta(self) -> None:
        """
        Sort supermodes primarily by solver number and secondarily by descending propagation constant (beta).
        """
        lexort_index = (
            [mode.solver_number for mode in self.supermodes],
            [-mode.beta.data[-1] for mode in self.supermodes],
        )

        self.all_supermodes = self._sort_modes(lexort_index)

    def is_compute_compatible(self, pair_of_mode: tuple) -> bool:
        """
        Determine whether the specified pair of modes is compatible for computation.

        Parameters
        ----------
        pair_of_mode : tuple
            The pair of modes to be checked.

        Returns
        -------
        bool
            True if the pair of modes is compatible for computation, False otherwise.
        """
        mode_0, mode_1 = pair_of_mode
        return mode_0.is_computation_compatible(mode_1)

    def remove_duplicate_combination(self, supermodes_list: list) -> list[SUPERMODE]:
        """
        Remove duplicate combinations from the mode combination list irrespective of the order.

        Parameters
        ----------
        supermodes_list : list of tuple
            List of mode combinations.

        Returns
        -------
        list of tuple
            Reduced list of unique supermode combinations.
        """
        output_list = []

        for mode0, mode1 in supermodes_list:
            if (mode0, mode1) not in output_list and (mode1, mode0) not in output_list:
                output_list.append((mode0, mode1))

        return output_list

    def interpret_combination(self, mode_of_interest: list, combination: str) -> set:
        """
        Interpret user input for mode selection and return the combination of modes to consider.

        Parameters
        ----------
        mode_of_interest : list
            List of modes of interest.
        combination : str
            Mode selection method ('pairs' or 'specific').

        Returns
        -------
        set of tuple
            Set of mode combinations.
        """
        self.test_valid_input(
            variable_name="combination",
            user_input=combination,
            valid_inputs=["pairs", "specific"],
        )

        match combination:
            case "pairs":
                mode_combinations = product(mode_of_interest, mode_of_interest)
            case "specific":
                mode_combinations = product(mode_of_interest, self.supermodes)

        mode_combinations = filter(self.is_compute_compatible, mode_combinations)

        mode_combinations = self.remove_duplicate_combination(mode_combinations)

        return set(mode_combinations)

    @parse_filename
    def save_instance(self, filename: str) -> Path:
        """
        Save the SuperSet instance as a serialized pickle file.

        Parameters
        ----------
        filename : str
            Filename for the serialized instance.

        Returns
        -------
        Path
            The path to the saved instance file.
        """
        with open(filename.with_suffix(".pickle"), "wb") as output_file:
            pickle.dump(self, output_file, pickle.HIGHEST_PROTOCOL)

        return filename

    @parse_filename
    def export_data(
        self,
        filename: str,
        mode_of_interest: list[SUPERMODE] | str = "all",
        combination: list | str = "pairs",
    ) -> Path:
        """
        Export the SuperSet data as CSV files, saving specific attributes of the modes or combinations of modes.

        Parameters
        ----------
        filename : str
            The directory where the files will be saved.
        mode_of_interest : list, optional
            List of modes to be exported (default is 'all').
        combination : list, optional
            List of mode combinations to be exported (default is None).

        Returns
        -------
        Path
            The path to the directory where the files were saved.
        """
        mode_of_interest = self.interpret_mode_of_interest(
            superset=self, mode_of_interest=mode_of_interest
        )

        combination = self.interpret_combination(
            mode_of_interest=mode_of_interest, combination=combination
        )

        # Create the directory if it doesn't exist
        output_dir = Path(filename)
        output_dir.mkdir(parents=True, exist_ok=True)

        def _export_single_data(attribute_name: str):
            """Helper function to export data for a single attribute for all modes of interest."""
            for mode in mode_of_interest:
                data = getattr(mode, attribute_name).data
                sub_filename = output_dir / f"{attribute_name}_{mode.label}".replace(
                    "}", ""
                ).replace("{", "")
                data = numpy.vstack([self.model_parameters.itr_list, data])
                numpy.savetxt(
                    fname=sub_filename.with_suffix(".csv"),
                    X=data,
                    delimiter=",",
                    header=f"ITR, {attribute_name}",
                )

        def _export_combination_data(attribute_name: str):
            """Helper function to export data for a single attribute for all combinations of modes."""
            for mode_0, mode_1 in combination:
                data = getattr(mode_0, attribute_name).get_values(
                    other_supermode=mode_1
                )
                sub_filename = (
                    output_dir
                    / f"{attribute_name}_{mode_0.label}_{mode_1.label}".replace(
                        "}", ""
                    ).replace("{", "")
                )
                data = numpy.vstack([self.model_parameters.itr_list, data])
                numpy.savetxt(
                    fname=sub_filename.with_suffix(".csv"),
                    X=data,
                    delimiter=",",
                    header=f"ITR, {attribute_name}",
                )

        # Export single-mode attributes
        _export_single_data("index")
        _export_single_data("beta")
        _export_single_data("eigenvalue")

        # Export combination-mode attributes
        _export_combination_data("adiabatic")
        _export_combination_data("beating_length")
        _export_combination_data("normalized_coupling")
