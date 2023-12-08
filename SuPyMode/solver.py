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
    extrapolation_order: int = 2
    """ Order of the taylor serie to extrapolate next eigenvalues . """
    debug_mode: int = 1
    """ Level of debug outprint from the c++ binding [0, 1, 2] """

    def __post_init__(self):
        self.mode_number = 0
        self.solver_number = 0

    @property
    def coordinate_system(self) -> object:
        return self.geometry.coordinate_system

    def initialize_binding(
            self,
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
            self.FD.triplet.array[:, 1],
            self.FD.triplet.array[:, 0],
            self.FD.triplet.array[:, 2]
        ]

        Solver = CppSolver(
            mesh=self.geometry.mesh,
            gradient=self.geometry.n2_gradient * self.geometry.coordinate_system.rho_mesh,
            itr_list=self.itr_list,
            finit_matrix=new_array.T,
            n_computed_mode=n_sorted_mode + n_added_mode,
            n_sorted_mode=n_sorted_mode,
            max_iter=self.max_iter,
            tolerance=self.tolerance,
            wavelength=self.wavelength,
            debug_mode=self.debug_mode,
            dx=self.geometry.coordinate_system.dx,
            dy=self.geometry.coordinate_system.dy
        )

        Solver.compute_laplacian()

        return Solver

    def init_superset(
            self,
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

        :returns:   No returns
        :rtype:     None
        """
        self.wavelength = wavelength
        self.wavenumber = 2 * numpy.pi / wavelength
        self.itr_list = numpy.linspace(itr_initial, itr_final, n_step)
        self.superset = SuperSet(parent_solver=self, wavelength=wavelength)

    def index_to_eigen_value(self, index: float) -> float:
        """
        Converts an effective index to the equivalent eigen value of the
        linear system to be solved.

        :param      eigen_value:  The eigen value
        :type       eigen_value:  float

        :returns:   The equivalent eigen value
        :rtype:     float
        """
        return -(index * self.wavenumber)**2

    def eigen_value_to_index(self, eigen_value: float) -> float:
        """
        Converts an eigen value of the linear equation to solve to
        its equivalent effective index.

        :param      eigen_value:  The eigen value
        :type       eigen_value:  float

        :returns:   The equivalent eigen value
        :rtype:     float
        """
        return numpy.sqrt(eigen_value) / self.wavenumber

    def get_supermode_labels(self, n_modes: int, boundaries: Boundaries2D, auto_label: bool) -> list:
        """
        Generate and returns the supermode label depending if auto_label
        is activated or not.

        :param      n_modes:     The n modes
        :type       n_modes:     int
        :param      boundaries:  The boundaries
        :type       boundaries:  Boundaries2D
        :param      auto_label:  The automatic label option
        :type       auto_label:  bool

        :returns:   The supermode labels.
        :rtype:     list
        """
        if auto_label:
            supermode_labels = ModeLabel(boundaries=boundaries, n_mode=n_modes).get_labels()
        else:
            supermode_labels = ["mode_" + "{" + str(n) + "}" for n in range(n_modes)]

        return supermode_labels

    def add_modes(
            self,
            n_sorted_mode: int,
            boundaries: dict,
            n_added_mode: int = 4,
            index_guess: float = 0.,
            auto_label: bool = True) -> None:
        """
        This methodes compute new set of n_added_mode modes for a given boundaries condition.
        It appends those modes to the one already computed.
        The auto_labeling options works only for almost cylindrical symmetric structure with low itr.
        If wrong label is settle it can be modified with the label_supermode method. Index guess is the effective index
        guess given to the inverse shift power method solver to retrieve the modes with close effective index.
        Adds modes to the superset instance. SuperSet is accessible through .get_set().

        :param      n_sorted_mode:    Number of mode that are outputed by the c++ solver.
        :type       n_sorted_mode:    int
        :param      boundaries:       Boundaries of the finit-difference system.
        :type       boundaries:       Boundaries2D
        :param      n_added_mode:     Number of additional modes that are computed for that will be sorted out, the higher the value the less likely mode mismatch will occur.
        :type       n_added_mode:     int
        :param      index_guess:      Initial effective index guess (if 0, auto evaluated).
        :type       index_guess:      float

        :returns:   No returns
        :rtype:     None
        """
        alpha = self.index_to_eigen_value(index_guess)

        cpp_solver = self.initialize_binding(
            boundaries=boundaries,
            n_added_mode=n_added_mode,
            n_sorted_mode=n_sorted_mode
        )

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
                binded_supermode=cpp_solver.get_mode(binding_number),
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
