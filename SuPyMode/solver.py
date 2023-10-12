#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Third-party imports
import numpy
from dataclasses import dataclass, field
from PyFinitDiff.sparse2D import FiniteDifference2D

# Local imports
from SuPyMode.superset import SuperSet
from SuPyMode.supermode import SuperMode
from PyFinitDiff.boundaries import Boundaries2D
from SuPyMode.binary.CppSolver import CppSolver
from SuPyMode.binary.SuperMode import SuperMode as BindingSuperMode  # It has to be imported in order for pybind11 to know the type
from SuPyMode.tools.mode_label import ModeLabel


@dataclass()
class SuPySolver(object):
    """
    .. note::
        This solver class directly links to a c++ Eigensolver.
        It solves the eigenvalues problems for a given geometry and return a collection of SuperModes.

    """
    geometry: None = field(repr=False)
    """ geometry of the coupler structure """
    tolerance: float = 1e-8
    """ Absolute tolerance on the propagation constant computation """
    max_iter: int = 10000
    """ Maximum iteration for the c++ Eigensolver """
    accuracy: int = 2
    """ Accuracy of the finit difference methode """
    show_iteration: bool = True
    """ Print option. """
    show_eigenvalues: bool = False
    """ Print option. """
    extrapolation_order: int = 1
    """ Order of the taylor serie to extrapolate next eigenvalues . """

    def __post_init__(self):
        self.mode_number = 0
        self.solver_number = 0

    @property
    def axes(self) -> object:
        return self.geometry.coordinate_system

    def get_set(self):
        return self.superset

    def initialize_binding(self,
            wavelength: float,
            n_sorted_mode: int,
            boundaries: Boundaries2D,
            n_added_mode: int) -> CppSolver:
        """
        Initializes the c++ binding of the class.

        :param      wavelength:       Wavelenght for the mode computation
        :type       wavelength:       float
        :param      n_sorted_mode:    Number of mode that are outputed by the c++ solver.
        :type       n_sorted_mode:    int
        :param      n_added_mode:     Number of additional modes that are computed for that will be sorted out, the higher the value the less likely mode mismatch will occur.
        :type       n_added_mode:     int
        :param      boundaries:       Symmetries of the finit-difference system.
        :type       boundaries:       Boundaries2D

        :returns:   The cpp solver.
        :rtype:     CppSolver
        """
        self.geometry.generate_coordinate_mesh_gradient()

        self.FD = FiniteDifference2D(
            n_x=self.geometry.coordinate_system.nx,
            n_y=self.geometry.coordinate_system.ny,
            dx=self.geometry.coordinate_system.dx,
            dy=self.geometry.coordinate_system.dy,
            derivative=2,
            accuracy=self.accuracy,
            boundaries=boundaries
        )

        new_array = numpy.c_[
            self.FD.triplet._array[:, 1],
            self.FD.triplet._array[:, 0],
            self.FD.triplet._array[:, 2]
        ]

        Solver = CppSolver(
            mesh=self.geometry.mesh,
            gradient=self.geometry.gradient * self.geometry.coordinate_system.rho_mesh,
            itr_list=self.itr_list,
            finit_matrix=new_array.T,
            n_computed_mode=n_sorted_mode + n_added_mode,
            n_sorted_mode=n_sorted_mode,
            max_iter=self.max_iter,
            tolerance=self.tolerance,
            wavelength=wavelength,
            show_iteration=self.show_iteration,
            show_eigenvalues=self.show_eigenvalues,
            dx=self.geometry.coordinate_system.dx,
            dy=self.geometry.coordinate_system.dy
        )

        Solver.compute_laplacian()

        return Solver

    def init_superset(self,
            wavelength: float,
            n_step: int = 300,
            itr_initial: float = 1.0,
            itr_final: float = 0.1) -> None:
        """
        Initialize superset instance which contains the computed superodes.

        :param      wavelength:  Wavelenght for the mode computation
        :type       wavelength:  float
        :param      n_step:      Number of stop to iterate through the ITR (inverse taper ration) section.
        :type       n_step:      int
        :param      itr_initial: Initial value of ITR.
        :type       itr_initial: float
        :param      itr_final:   Final value of ITR.
        :type       itr_final:   float
        """
        self.wavelength = wavelength
        self.wavenumber = 2 * numpy.pi / wavelength
        self.itr_list = numpy.linspace(itr_initial, itr_final, n_step)
        self.superset = SuperSet(parent_solver=self, wavelength=wavelength)

    def index_to_eigen_value(self, index):
        return -(index * self.wavenumber)**2

    def eigen_value_to_index(self, eigen_value):
        return numpy.sqrt(eigen_value) / self.wavenumber

    def add_modes(self,
            n_sorted_mode: int,
            boundaries: dict,
            n_added_mode: int = 4,
            auto_labeling: bool = False,
            index_guess: float = 0.,
            auto_label: bool = True) -> None:
        """
        Adds modes to the superset instance. SuperSet is accessible through .get_set().

        :param      n_sorted_mode:    Number of mode that are outputed by the c++ solver.
        :type       n_sorted_mode:    int
        :param      boundaries:       Boundaries of the finit-difference system.
        :type       boundaries:       Boundaries2D
        :param      n_added_mode:     Number of additional modes that are computed for that will be sorted out, the higher the value the less likely mode mismatch will occur.
        :type       n_added_mode:     int
        :param      index_guess:      Initial effective index guess (if 0, auto evaluated).
        :type       index_guess:      float
        """
        alpha = self.index_to_eigen_value(index_guess)

        cpp_solver = self.initialize_binding(
            boundaries=boundaries,
            wavelength=self.wavelength,
            n_added_mode=n_added_mode,
            n_sorted_mode=n_sorted_mode
        )

        cpp_solver.loop_over_itr(
            extrapolation_order=self.extrapolation_order,
            alpha=alpha
        )

        if auto_label:
            mode_labels = ModeLabel(boundaries=boundaries, n_mode=n_sorted_mode).get_labels()
        else:
            mode_labels = ["mode_" + "{" + str(n) + "}" for n in range(n_sorted_mode)]

        for binding_number, label in enumerate(mode_labels):
            supermode = SuperMode(
                parent_set=self.superset,
                binded_supermode=cpp_solver.get_mode(binding_number),
                mode_number=self.mode_number,
                solver_number=self.solver_number,
                wavelength=self.wavelength,
                boundaries=boundaries,
                itr_list=self.itr_list,
                label=label
            )

            self.superset.supermodes.append(supermode)

            setattr(self.superset, label, supermode)

            self.mode_number += 1

        self.solver_number += 1

# ---
