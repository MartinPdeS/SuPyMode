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
from SuPyMode.binary.interface_eigensolver import EIGENSOLVER as EigenSolver  # type: ignore
from SuPyMode.binary.interface_model_parameters import MODELPARAMETERS as ModelParameters  # type: ignore
from SuPyMode.binary.interface_supermode import SUPERMODE  # type: ignore
from SuPyMode.mode_label import ModeLabel


@dataclass()
class SuPySolver(object):
    """
    A solver for computing eigenvalues and supermodes of optical fiber geometries using a C++ eigensolver.

    This class manages the eigenvalue problem for optical structures and returns computed supermodes. The solver utilizes
    a C++ backend for efficient eigenvalue computation and integrates with finite difference methods to solve for various
    boundary conditions.

    Parameters
    ----------
    geometry : Geometry or numpy.ndarray
        The refractive index geometry of the optical structure.
    tolerance : float, optional
        Absolute tolerance for the propagation constant computation (default is 1e-8).
    max_iter : int, optional
        Maximum iterations for the C++ eigensolver (default is 10,000).
    accuracy : int, optional
        Accuracy level of the finite difference method (default is 2).
    extrapolation_order : int, optional
        Order of Taylor series used to extrapolate eigenvalues (default is 2).
    debug_mode : int, optional
        Debug output level from the C++ binding, where 0 is no output and higher values provide more detail (default is 1).
    coordinate_system : CoordinateSystem, optional
        The coordinate system linked with the geometry. Must be provided if geometry is given as an array.
    """
    geometry: Geometry | numpy.ndarray = field(repr=False)
    tolerance: float = 1e-8
    max_iter: int = 10_000
    accuracy: int = 2
    extrapolation_order: int = 2
    debug_mode: int = 1
    coordinate_system: CoordinateSystem | None = None

    def __post_init__(self):
        self.mesh = self.geometry.mesh
        self.coordinate_system = self.geometry.coordinate_system

        self.mode_number = 0
        self.solver_number = 0

    def initialize_binding(self, n_sorted_mode: int, boundaries: Boundaries, n_added_mode: int) -> EigenSolver:
        """
        Initialize and configure the C++ solver binding for eigenvalue computations.

        Parameters
        ----------
        n_sorted_mode : int
            Number of modes to sort and retrieve from the solver.
        boundaries : Boundaries
            Boundary conditions for the finite difference system.
        n_added_mode : int
            Number of extra modes calculated for accuracy and reliability.

        Returns
        -------
        EigenSolver
            Configured C++ solver instance.
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

        solver = EigenSolver(
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
        Initialize a SuperSet instance containing computed supermodes over a range of inverse taper ratios (ITR).

        Parameters
        ----------
        wavelength : float
            Wavelength for the mode computation.
        n_step : int, optional
            Number of steps for the ITR interpolation (default is 300).
        itr_initial : float, optional
            Initial ITR value (default is 1.0).
        itr_final : float, optional
            Final ITR value (default is 0.1).
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
        Convert a refractive index to the corresponding eigenvalue for the solver.

        Parameters
        ----------
        index : float
            Refractive index to convert.

        Returns
        -------
        float
            Calculated eigenvalue based on the given index and the wavenumber.
        """
        return -(index * self.wavenumber)**2

    def eigen_value_to_index(self, eigen_value: float) -> float:
        """
        Convert an eigenvalue from the solver to the corresponding refractive index.

        Parameters
        ----------
        eigen_value : float
            Eigenvalue to convert.

        Returns
        -------
        float
            Equivalent refractive index calculated from the eigenvalue and the wavenumber.
        """
        return numpy.sqrt(eigen_value) / self.wavenumber

    def get_supermode_labels(self, n_modes: int, boundaries: Boundaries, auto_label: bool) -> list:
        """
        Generate labels for supermodes based on boundary conditions and whether auto-labeling is enabled.

        Parameters
        ----------
        n_modes : int
            Number of modes for which labels are needed.
        boundaries : Boundaries
            Boundary conditions that affect mode symmetries.
        auto_label : bool
            If True, automatically generates labels based on mode symmetries; otherwise, generates generic labels.

        Returns
        -------
        list
            List of labels for the supermodes.
        """
        if auto_label:
            return [ModeLabel(boundaries=boundaries, mode_number=n).label for n in range(n_modes)]
        else:
            return [f"mode_{{{n}}}" for n in range(n_modes)]

    def add_modes(self, n_sorted_mode: int, boundaries: Boundaries, n_added_mode: int = 4, index_guess: float = 0., auto_label: bool = True) -> None:
        """
        Compute and add a specified number of supermodes to the solver's collection.

        Parameters
        ----------
        n_sorted_mode : int
            Number of modes to output and sort from the solver.
        boundaries : Boundaries
            Boundary conditions for the finite difference calculations.
        n_added_mode : int, optional
            Additional modes computed to ensure mode matching accuracy (default is 4).
        index_guess : float, optional
            Starting guess for the refractive index used in calculations (default is 0., auto-evaluated if set to 0).
        auto_label : bool, optional
            If True, enables automatic labeling of modes based on symmetry (default is True).

        Returns
        -------
        None
            This method updates the solver's internal state but does not return any value.
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
