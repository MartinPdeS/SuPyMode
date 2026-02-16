#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Third-party imports
import numpy
from PyFinitDiff.finite_difference_2D import FiniteDifference
from PyFinitDiff.finite_difference_2D import Boundaries

# Local imports
from SuPyMode.superset import SuperSet
from SuPyMode.binary.interface_supermode import SUPERMODE  # type: ignore
from SuPyMode.binary.interface_eigensolver import EIGENSOLVER  # type: ignore
from SuPyMode.binary.interface_model_parameters import ModelParameters  # type: ignore


class SuPySolver(EIGENSOLVER):
    """
    A solver for computing eigenvalues and supermodes of optical fiber geometries using a C++ eigensolver.

    This class manages the eigenvalue problem for optical structures and returns computed supermodes. The solver utilizes
    a C++ backend for efficient eigenvalue computation and integrates with finite difference methods to solve for various
    boundary conditions.

    Parameters
    ----------
    mesh : numpy.ndarray
        The refractive index geometry of the optical structure.
    x : numpy.ndarray
        The x-coordinate vector of the mesh.
    y : numpy.ndarray
        The y-coordinate vector of the mesh.
    tolerance : float, optional
        Absolute tolerance for the propagation constant computation (default is 1e-8).
    max_iteration : int, optional
        Maximum iterations for the C++ eigensolver (default is 10,000).
    accuracy : int, optional
        Accuracy level of the finite difference method (default is 2).
    extrapolation_order : int, optional
        Order of Taylor series used to extrapolate eigenvalues (default is 2).
    debug_mode : int, optional
        Debug output level from the C++ binding, where 0 is no output and higher values provide more detail (default is 1).
    """

    def __init__(
        self,
        mesh: numpy.typing.NDArray,
        x: numpy.typing.NDArray,
        y: numpy.typing.NDArray,
        tolerance: float = 1e-8,
        max_iteration: int = 10_000,
        accuracy: int = 2,
        extrapolation_order: int = 2,
        debug_mode: int = 1,
    ) -> None:
        """
        Initialize the solver with the given parameters.
        """
        self.mesh = mesh

        dx = abs(x[0] - x[1])
        dy = abs(y[0] - y[1])

        model_parameters = ModelParameters(
            dx=dx,
            dy=dy,
            mesh=mesh,
            x_vector=x,
            y_vector=y,
            debug_mode=debug_mode,
        )

        super().__init__(
            model_parameters=model_parameters,
            max_iteration=max_iteration,
            tolerance=tolerance,
            accuracy=accuracy,
            extrapolation_order=extrapolation_order,
        )

    def get_finit_difference_array(self, boundaries: Boundaries) -> None:
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
        finit_diff_matrix = FiniteDifference(
            n_x=self.model_parameters.nx,
            n_y=self.model_parameters.ny,
            dx=self.model_parameters.dx,
            dy=self.model_parameters.dy,
            derivative=2,
            accuracy=self.accuracy,
            boundaries=boundaries,
        )

        finit_diff_matrix.construct_triplet()

        return numpy.c_[
            finit_diff_matrix._triplet.array[:, 1],
            finit_diff_matrix._triplet.array[:, 0],
            finit_diff_matrix._triplet.array[:, 2],
        ]

    def init_superset(
        self,
        wavelength: float,
        n_step: int = 300,
        itr_initial: float = 1.0,
        itr_final: float = 0.1,
    ) -> None:
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
        self.model_parameters.initialize(
            wavelength=wavelength,
            itr_initial=itr_initial,
            itr_final=itr_final,
            n_step=n_step,
        )

        self.superset = SuperSet(
            model_parameters=self.model_parameters,
        )

    def add_modes(
        self,
        n_sorted_mode: int,
        boundaries: Boundaries,
        n_added_mode: int = 4,
        index_guess: float = 0.0,
        auto_label: bool = True,
    ) -> None:
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
        finit_difference_array = self.get_finit_difference_array(
            boundaries=boundaries,
        )

        self.reset_solver()

        self.set_boundaries(
            left=boundaries.left.value.lower(),
            right=boundaries.right.value.lower(),
            top=boundaries.top.value.lower(),
            bottom=boundaries.bottom.value.lower(),
        )

        self.initialize(
            model_parameters=self.model_parameters,
            finit_matrix=finit_difference_array.T,
            n_computed_mode=n_sorted_mode + n_added_mode,
            n_sorted_mode=n_sorted_mode,
        )

        self.loop_over_itr(guess_index=index_guess, auto_label=auto_label)

        for binding_number in range(n_sorted_mode):
            supermode = self.get_sorted_mode(binding_number)

            self.superset.supermodes.append(supermode)

            setattr(self.superset, supermode.label, supermode)
