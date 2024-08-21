#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Standard imports
import numpy
from dataclasses import dataclass

# Scipy imports
from scipy.sparse import linalg
from scipy.sparse._csr import csr_matrix

# Other imports
from SuPyMode.python_debuging.mode_solver import ModeSolver
from FiberFusing import Geometry
from MPSPlots.render2D import Scene2D, Axis, Line


@dataclass
class EigenSolver:
    laplacian: numpy.ndarray
    """ Laplacian sparse matrix that is used to compute the eigen solutions """
    geometry: Geometry
    """ Geometry object that contain the mesh used connjointly with the laplacian """
    itr_list: numpy.ndarray
    """ Inverse taper ratio list to which evaluate the solutions """
    wavelength: float
    """ Wavelength used for the simulation """
    max_iteration: int = 1000
    """ Maximum number of time the bicgstab is called per mode per slice """
    min_iteration: int = 5
    """ Minimum number of time the bicgstab is called per mode per slice """
    randomize_factor: float = 1e-15
    """ Factor for randomization of eiven vector guess from one slice to the next """
    n_solution: int = 1
    """ Number of mode to compute """
    tolerance: float = 1e-3,
    """ Tolerance for the eigen solution to compute """
    debug: bool = False
    """ Debug parameter, if enabled a lot of verbose will be show """
    debug_guess: bool = False
    """ Debug parameter, if enabled shows the eigen_value guess and computed value """
    debug_bicgstab: bool = False
    """ Debug parameter, if enabled shows the bicgstab debuging """
    fitting_max_degree: int = 5
    """ Fitting maximim degree for eigen value extrapolation """

    def __post_init__(self) -> None:
        self.wavenumber = 2 * numpy.pi / self.wavelength

        self.mode_solvers = []
        for n in range(self.n_solution):
            self.mode_solvers.append(
                ModeSolver(
                    itr_list=self.itr_list,
                    shape=self.geometry.mesh.shape,
                    wavelength=self.wavelength,
                    fitting_max_degree=self.fitting_max_degree,
                    label=f"mode_{n}"
                )
            )

    def __getitem__(self, idx):
        return self.mode_solvers[idx]

    def index_to_eigen_value(self, index: float, itr: float) -> float:
        """
        Return the eigen value associated to the given effective index

        :param      index:  The mode effective index
        :type       index:  float

        :returns:   The associated eigen value
        :rtype:     float
        """
        return (self.wavenumber * index * itr)**2

    def eigen_value_to_index(self, eigen_value: float, itr: float) -> float:
        """
        Return the effective index associated toe the given eigen value.

        :param      eigen_value:  The eigen value
        :type       eigen_value:  float
        :param      itr:          The itr
        :type       itr:          float

        :returns:   The associated effective index
        :rtype:     float
        """
        return numpy.sqrt(eigen_value) / self.wavenumber / itr

    def get_vector_norm(self, vector: numpy.ndarray) -> float:
        """
        Returns the norm of the norm2 of the given vector.

        :param      vector:  The vector
        :type       vector:  float
        """
        return numpy.sqrt((vector**2).sum())

    def normalize_vector(self, vector: numpy.ndarray) -> numpy.ndarray:
        """
        Returns the normalized vector given. The norm is norm2

        :param      vector:  The vector to normalized
        :type       vector:  numpy.ndarray

        :returns:   The normalized vector
        :rtype:
        """
        return vector / self.get_vector_norm(vector)

    def generate_eigen_matrix(self, itr: float = 1) -> csr_matrix:
        """
        Returns a eigen matrix constitued of the laplacian with the diagonal added of the
        mesh value times the wavenumbered squareds.

        :param      itr:  The itr to which evaluate the eigen matrix
        :type       itr:  float
        """
        wavenumber = self.wavenumber * itr
        eigen_matrix = self.laplacian.copy()

        eigen_matrix.setdiag(
            eigen_matrix.diagonal() + (self.geometry.mesh.ravel() * wavenumber)**2
        )

        return eigen_matrix

    def get_eigen_solution_from_scipy(self, eigen_matrix: csr_matrix, index_guess: float, itr: float) -> tuple:
        """
        Return eigen solution for the given eigen matrix input and index guess.

        :param      eigen_matrix:  The eigen matrix
        :type       eigen_matrix:  csr_matrix
        :param      index_guess:   The index guess
        :type       index_guess:   float

        :returns:   The eigen solution from scipy.
        :rtype:     tuple
        """
        eigen_value_guess = self.index_to_eigen_value(index=index_guess, itr=itr)

        print(eigen_value_guess)
        eigen_values, eigen_vectors = linalg.eigs(
            eigen_matrix,
            k=self.n_solution + 5,
            which='LM',
            sigma=eigen_value_guess
        )

        return eigen_values[:self.n_solution], eigen_vectors[:, :self.n_solution]

    def compute_eigen_solution_with_scipy(self, number_of_iteration: int = 1, index_guess: float = None):
        self.pre_computed_slice = 0

        if index_guess is None:
            index_guess = self.geometry.mesh.max()

        for itr in self.itr_list[:number_of_iteration]:
            eigen_matrix = self.generate_eigen_matrix(itr)

            values, vectors = self.get_eigen_solution_from_scipy(
                eigen_matrix=eigen_matrix,
                index_guess=index_guess,
                itr=itr
            )

            for value, vector, mode_solver in zip(values, vectors.T, self.mode_solvers):
                mode_solver.append_next_slice(
                    eigen_vector=vector,
                    eigen_value=value,
                    itr=itr
                )

            self.pre_computed_slice += 1

    def shift_sparse_matrix(self, sparse_matrix: csr_matrix, shift: float) -> csr_matrix:
        """
        Retunrs the eigen matrix given to the same but with the diagonal shifted.

        :param      sparse_matrix:  The sparse matrix
        :type       sparse_matrix:  csr_matrix
        :param      shift:          The shift
        :type       shift:          float
        """
        shifted_sparse_matrix = sparse_matrix.copy()

        shifted_sparse_matrix.setdiag(
            shifted_sparse_matrix.diagonal() - shift
        )

        return shifted_sparse_matrix

    def compute_next_slice(self):
        itr = self.itr_list[self.pre_computed_slice]
        eigen_matrix = self.generate_eigen_matrix(itr=itr)

        for mode_solver in self.mode_solvers:
            if self.debug:
                print(f"Mode processed: {mode_solver.label}")

            eigen_value_guess = mode_solver.extrapolate_eigen_value(itr=itr)
            eigen_vector_guess = mode_solver.extrapolate_eigen_vector(itr=itr)

            eigen_value, eigen_vector = self.power_shift_method(
                eigen_matrix=eigen_matrix,
                eigen_value_guess=eigen_value_guess,
                eigen_vector_guess=eigen_vector_guess
            )

            mode_solver.append_eigen_vector(eigen_vector=eigen_vector)
            mode_solver.append_eigen_value(eigen_value=eigen_value, itr=itr)

        self.pre_computed_slice += 1

    def power_shift_method(self, eigen_matrix: csr_matrix, eigen_value_guess: float = 0, eigen_vector_guess: numpy.ndarray = None) -> tuple:
        """
        Use the inverse power shift method to retrieve the eigen solution
        for a given eigen matrix.
        The method accept an initial eigen value and vector guess

        :param      eigen_matrix:        The eigen matrix
        :type       eigen_matrix:        csr_matrix
        :param      eigen_value_guess:   The eigen value guess
        :type       eigen_value_guess:   float
        :param      eigen_vector_guess:  The eigen vector guess
        :type       eigen_vector_guess:  numpy.ndarray

        :returns:   Eigen value and eigen vector in a tuple
        :rtype:     tuple
        """
        if eigen_vector_guess is None:
            eigen_vector_guess = numpy.random.rand(eigen_matrix.shape[0])

        eigen_vector_guess += numpy.random.rand(*eigen_vector_guess.shape) * self.randomize_factor

        shifted_A = self.shift_sparse_matrix(eigen_matrix, eigen_value_guess)

        solution_k = self.normalize_vector(eigen_vector_guess)

        for iteration in range(self.max_iteration):

            v = self.normalize_vector(solution_k)

            solution_k, _ = linalg.bicgstab(shifted_A, v.real)

            projection = v.dot(solution_k)

            residual = self.get_vector_norm(solution_k - projection * v) / abs(projection)

            eigen_value = eigen_value_guess + 1 / projection
            eigen_vector = solution_k / projection

            if iteration >= self.min_iteration and residual < self.tolerance:
                if self.debug_bicgstab:
                    print(f"Process suceed {residual = :.4e} in {iteration = } iteration")
                return eigen_value, eigen_vector

        print(f" ===> Process failed {residual = :.4e} in {iteration = } iteration")
        return eigen_value, eigen_vector

    def process_mode(self, mode_number: int) -> None:
        """
        Compute the full itr_list iteration for a given mode number.

        :param      mode_number:  The mode number
        :type       mode_number:  int

        :returns:   No return
        :rtype:     None
        """
        mode = self.mode_solvers[mode_number]

        if self.debug:
            print(f"Mode processed: {mode.label}")

        for itr in self.itr_list[self.pre_computed_slice:]:
            if self.debug:
                print(f"itr: {itr}")

            eigen_matrix = self.generate_eigen_matrix(itr=itr)

            eigen_value_guess = mode.extrapolate_eigen_value(itr=itr)

            eigen_vector_guess = mode.last_eigen_vector

            eigen_value, eigen_vector = self.power_shift_method(
                eigen_matrix=eigen_matrix,
                eigen_value_guess=eigen_value_guess,
                eigen_vector_guess=eigen_vector_guess
            )

            if self.debug_guess:
                print(f'eigen value guess: {eigen_value_guess:.8e}\tcomputed eigen value: {eigen_value:.8e}\n')

            mode.append_eigen_vector(eigen_vector=eigen_vector)
            mode.append_eigen_value(eigen_value=eigen_value, itr=itr)

    def render_index_on_ax(self, ax: Axis):
        for mode in self.mode_solvers:
            artist = Line(
                x=mode.itr_list,
                y=mode.effective_index,
            )

            ax.add_artist(artist)

            ax.x_label = 'Inverse taper ratio'
            ax.y_label = 'Effective index'

    def plot(self):
        figure = Scene2D(unit_size=(15, 6))

        ax = Axis(row=0, col=0)

        figure.add_axes(ax)

        self.render_index_on_ax(ax=ax)

        return figure

# -
