#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Third-party imports
import numpy
from dataclasses import dataclass, field
from PyFinitDiff.finite_difference_2D import FiniteDifference
from PyFinitDiff.finite_difference_2D import Boundaries
from FiberFusing.geometry import Geometry
from FiberFusing.coordinate_system import CoordinateSystem

# Local imports
from SuPyMode.superset import SuperSet
from SuPyMode.supermode import SuperMode
from SuPyMode.binary.CppSolver import CppSolver
from SuPyMode.binary.ModelParameters import ModelParameters
from SuPyMode.binary.SuperMode import SuperMode as BindingSuperMode  # noqa: F401 It has to be imported in order for pybind11 to know the type
from SuPyMode.mode_label import ModeLabel


@dataclass()
class SuPySolver(object):
    """
    Solver class integrating a C++ eigensolver to compute eigenvalues for optical fiber geometries.
    This class manages the eigenvalue problems and returns collections of computed SuperModes.

    Attributes:
        geometry (Geometry | np.ndarray): The refractive index geometry of the optical structure.
        tolerance (float): Absolute tolerance for the propagation constant computation.
        max_iter (int): Maximum iterations for the C++ eigensolver.
        accuracy (int): Accuracy level of the finite difference method.
        extrapolation_order (int): Order of Taylor series used to extrapolate eigenvalues.
        debug_mode (int): Debug output level from the C++ binding (0, 1, 2).
        coordinate_system (Optional[CoordinateSystem]): The coordinate system linked with the geometry.
    """
    geometry: Geometry | numpy.ndarray = field(repr=False)
    tolerance: float = 1e-8
    max_iter: int = 10_000
    accuracy: int = 2
    extrapolation_order: int = 2
    debug_mode: int = 1
    coordinate_system: CoordinateSystem | None = None

    def __post_init__(self):
        if isinstance(self.geometry, numpy.ndarray):
            assert self.coordinate_system is not None, "Geometry provided without its coordinate system"
            self.mesh = self.geometry
        else:
            self.geometry.generate_coordinate_system()
            self.mesh = self.geometry.generate_mesh()
            self.coordinate_system = self.geometry.coordinate_system

        self.mode_number = 0
        self.solver_number = 0

    def initialize_binding(self, n_sorted_mode: int, boundaries: Boundaries, n_added_mode: int) -> CppSolver:
        """
        Initializes and configures the C++ solver binding for eigenvalue computations.

        Args:
            n_sorted_mode (int): Number of modes to sort and retrieve from the solver.
            boundaries (Boundaries): Boundary conditions for the finite difference system.
            n_added_mode (int): Number of extra modes calculated for accuracy and reliability.

        Returns:
            CppSolver: Configured C++ solver instance.
        """
        self.FD = FiniteDifference(
            n_x=self.coordinate_system.nx,
            n_y=self.coordinate_system.ny,
            dx=self.coordinate_system.dx,
            dy=self.coordinate_system.dy,
            derivative=2,
            accuracy=self.accuracy,
            boundaries=boundaries
        )

        self.FD.construct_triplet()

        new_array = numpy.c_[
            self.FD._triplet.array[:, 1],
            self.FD._triplet.array[:, 0],
            self.FD._triplet.array[:, 2]
        ]

        solver = CppSolver(
            model_parameters=self.model_parameters,
            finit_matrix=new_array.T,
            n_computed_mode=n_sorted_mode + n_added_mode,
            n_sorted_mode=n_sorted_mode,
            max_iter=self.max_iter,
            tolerance=self.tolerance,
            left_boundary=boundaries.left,
            right_boundary=boundaries.right,
            top_boundary=boundaries.top,
            bottom_boundary=boundaries.bottom,
        )

        solver.compute_laplacian()

        return solver

    def init_superset(self, wavelength: float, n_step: int = 300, itr_initial: float = 1.0, itr_final: float = 0.1) -> None:
        """
        Initializes a SuperSet instance containing computed supermodes over a range of inverse taper ratios (ITR).

        Args:
            wavelength (float): Wavelength for the mode computation.
            n_step (int): Number of steps for the ITR interpolation.
            itr_initial (float): Initial ITR value.
            itr_final (float): Final ITR value.
        """
        self.wavelength = wavelength
        self.wavenumber = 2 * numpy.pi / wavelength
        self.itr_list = numpy.linspace(itr_initial, itr_final, n_step)

        self.model_parameters = ModelParameters(
            dx=self.coordinate_system.dx,
            dy=self.coordinate_system.dy,
            wavelength=wavelength,
            itr_list=self.itr_list,
            mesh=self.mesh,
            x_vector=self.coordinate_system.x_vector,
            y_vector=self.coordinate_system.y_vector,
            debug_mode=self.debug_mode
        )

        self.superset = SuperSet(geometry=self.geometry, wavelength=wavelength, model_parameters=self.model_parameters)

    def index_to_eigen_value(self, index: float) -> float:
        """
        Converts a refractive index to the corresponding eigenvalue for the solver.

        Args:
            index (float): Refractive index to convert.

        Returns:
            float: Calculated eigenvalue based on the given index and the wavenumber.
        """
        return -(index * self.wavenumber)**2

    def eigen_value_to_index(self, eigen_value: float) -> float:
        """
        Converts an eigenvalue from the solver to the corresponding refractive index.

        Args:
            eigen_value (float): Eigenvalue to convert.

        Returns:
            float: Equivalent refractive index calculated from the eigenvalue and the wavenumber.
        """
        return numpy.sqrt(eigen_value) / self.wavenumber

    def get_supermode_labels(self, n_modes: int, boundaries: Boundaries, auto_label: bool) -> list:
        """
        Generates labels for supermodes based on boundary conditions and whether auto-labeling is enabled.

        Args:
            n_modes (int): Number of modes for which labels are needed.
            boundaries (Boundaries): Boundary conditions that affect mode symmetries.
            auto_label (bool): If True, automatically generates labels based on mode symmetries; otherwise, generates generic labels.

        Returns:
            list: List of labels for the supermodes.
        """
        if auto_label:
            return [ModeLabel(boundaries=boundaries, mode_number=n).label for n in range(n_modes)]
        else:
            return ["mode_" + "{" + str(n) + "}" for n in range(n_modes)]

    def add_modes(self, n_sorted_mode: int, boundaries: Boundaries, n_added_mode: int = 4, index_guess: float = 0., auto_label: bool = True) -> None:
        """
        Computes and adds a specified number of supermodes to the solver's collection, using given boundary conditions and mode sorting criteria.

        Args:
            n_sorted_mode (int): Number of modes to output and sort from the solver.
            boundaries (Boundaries): Boundary conditions for the finite difference calculations.
            n_added_mode (int): Additional modes computed to ensure mode matching accuracy.
            index_guess (float): Starting guess for the refractive index used in calculations (if 0, auto evaluated).
            auto_label (bool): If True, enables automatic labeling of modes based on symmetry.

        Returns:
            None: This method updates the solver's internal state but does not return any value.
        """
        alpha = self.index_to_eigen_value(index_guess)

        cpp_solver = self.initialize_binding(
            boundaries=boundaries,
            n_added_mode=n_added_mode,
            n_sorted_mode=n_sorted_mode
        )

        self.superset.model_parameters = self.model_parameters

        cpp_solver.loop_over_itr(
            extrapolation_order=self.extrapolation_order,
            alpha=alpha
        )

        mode_labels = self.get_supermode_labels(
            n_modes=n_sorted_mode,
            boundaries=boundaries,
            auto_label=auto_label
        )

        for binding_number, label in enumerate(mode_labels):
            supermode = SuperMode(
                parent_set=self.superset,
                binding=cpp_solver.get_mode(binding_number),
                mode_number=self.mode_number,
                solver_number=self.solver_number,
                boundaries=boundaries,
                label=label
            )

            self.superset.supermodes.append(supermode)

            setattr(self.superset, label, supermode)

            self.mode_number += 1

        self.solver_number += 1

# ---
