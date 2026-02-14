#include <pybind11/pybind11.h>
#include <memory>
#include "eigensolver.h"

/**
 * @brief Python bindings for the EigenSolver C++ module
 *
 * This module provides Python bindings for SuPyMode's eigenmode solver,
 * enabling efficient computation of optical fiber modes using finite
 * difference methods and eigenvalue decomposition.
 */
PYBIND11_MODULE(interface_eigensolver, module)
{
    module.doc() = R"pbdoc(
        SuPyMode EigenSolver Module

        This module provides high-performance eigenmode computation for optical
        waveguides using finite difference methods. It solves the eigenvalue
        problem for the wave equation in cylindrical coordinates to find the
        propagation constants and field distributions of guided modes.

        The solver supports various boundary conditions and can handle both
        step-index and graded-index fiber profiles.
    )pbdoc";

    pybind11::class_<EigenSolver, std::shared_ptr<EigenSolver>>(module, "EIGENSOLVER",
        R"pbdoc(
            Eigenmode solver for optical waveguides.

            This class implements a finite difference eigenvalue solver for computing
            the modes of optical fibers and waveguides. It uses advanced numerical
            methods to find eigenvalues (propagation constants) and eigenvectors
            (mode field distributions) of the wave equation.

            The solver supports:
            - Step-index and graded-index profiles
            - Various boundary conditions (PML, PEC, PMC)
            - Efficient sparse matrix operations
            - Mode sorting and selection
        )pbdoc"
    )
    .def(pybind11::init<const size_t, const double, const size_t, const size_t>(),
        pybind11::arg("max_iteration") = 10'000,
        pybind11::arg("tolerance") = 1e-8,
        pybind11::arg("accuracy") = 2,
        pybind11::arg("extrapolation_order") = 2,
        R"pbdoc(
            Default constructor for EigenSolver.

            Initializes an empty solver instance. Use the initialize method to set
            parameters before solving eigenvalue problems.

            Parameters
            ----------
            max_iteration : int
                Maximum number of iterations for the eigenvalue solver.
            tolerance : float
                Convergence tolerance for the eigenvalue solver.
            accuracy : int
                Accuracy order of the finite difference scheme (e.g., 1, 2, 4)
            extrapolation_order : int
                Order of Richardson extrapolation for iterative eigenvalue computation (e.g., 1, 2, 3)
        )pbdoc"
    )
    .def_readonly(
        "accuracy",
        &EigenSolver::accuracy,
        R"pbdoc(
            Accuracy order of the finite difference scheme.

            This parameter controls the order of accuracy for the finite difference
            approximation used in the eigenvalue solver. Higher values provide better
            accuracy at the cost of increased computational complexity.

            Common values:
            - 1: First-order accuracy (simplest, least accurate)
            - 2: Second-order accuracy (standard choice for many problems)
            - 4: Fourth-order accuracy (more accurate, especially for smooth solutions)

            The optimal choice depends on the specific problem and desired balance
            between accuracy and computational resources.
        )pbdoc"
    )
    .def_readwrite(
        "extrapolation_order",
        &EigenSolver::extrapolation_order,
        R"pbdoc(
            Extrapolation order for iterative eigenvalue computation.
            This parameter specifies the order of Richardson extrapolation used in
            the iterative eigenvalue solver. Extrapolation can accelerate convergence
            and improve accuracy, especially for challenging problems.

            Recommended values:
            - 1: No extrapolation (standard iterative solver)
            - 2: Second-order extrapolation (common choice for improved convergence)
            - 3 or higher: Higher-order extrapolation (can provide better accuracy but may introduce instability)

        )pbdoc"
    )
    .def_readwrite(
        "wavelength",
        &EigenSolver::wavelength,
        R"pbdoc(
            Operating wavelength for the eigenmode solver.

            This parameter specifies the wavelength of light for which the modes are
            being computed. It is used in the calculation of the propagation constants
            and mode field distributions.

            Units: meters (m)

            The wavelength should be chosen based on the application and the physical
            properties of the waveguide being modeled.
        )pbdoc"
    )
    .def_readwrite(
        "mode_number",
        &EigenSolver::mode_number,
        R"pbdoc(
            Current mode number being processed.
            This parameter keeps track of the mode index currently being computed or
            sorted. It is used internally to manage the mode computation process and
            can be accessed for informational purposes.
        )pbdoc"
    )
    .def_readwrite(
        "solver_number",
        &EigenSolver::solver_number,
        R"pbdoc(
            Current solver iteration number.
            This parameter tracks the current iteration of the eigenvalue solver, which
            can be useful for monitoring convergence and debugging. It is updated
            during the iterative eigenvalue computation process.
        )pbdoc"
    )
    .def(
        "index_to_eigenvalue",
        &EigenSolver::index_to_eigenvalue,
        pybind11::arg("mode_index"),
        R"pbdoc(
            Convert mode index to eigenvalue (propagation constant).
            This method takes a mode index (e.g., 0 for the fundamental mode, 1 for
            the first higher-order mode, etc.) and returns the corresponding eigenvalue
            (propagation constant) for that mode.

            Parameters
            ----------
            mode_index : double
                The refractive index of the mode for which to compute the eigenvalue.

            Returns
            -------
            complex
                The eigenvalue (propagation constant) corresponding to the given mode index.
        )pbdoc"
    )
    .def("eigenvalue_to_index",
        &EigenSolver::eigenvalue_to_index,
        pybind11::arg("eigenvalue"),
        R"pbdoc(
            Convert eigenvalue (propagation constant) to mode index.

            This method takes an eigenvalue (propagation constant) and returns the
            corresponding mode index, which can be used to identify the mode order
            (e.g., 0 for fundamental mode, 1 for first higher-order mode, etc.).

            Parameters
            ----------
            eigenvalue : complex
                The eigenvalue (propagation constant) for which to compute the mode index.

            Returns
            -------
            double
                The mode index corresponding to the given eigenvalue.
        )pbdoc"
    )
    .def("_cpp_set_boundaries",
        &EigenSolver::setup_boundaries,
        pybind11::arg("left"),
        pybind11::arg("right"),
        pybind11::arg("top"),
        pybind11::arg("bottom"),
        R"pbdoc(
            Set the boundary conditions for the eigenvalue solver.

            This method configures the boundary conditions for the waveguide
            model. It should be called after initializing the solver and before
            running the eigenvalue computation.

            Parameters
            ----------
            left : str
                Boundary condition for the left edge.
            right : str
                Boundary condition for the right edge.
            top : str
                Boundary condition for the top edge.
            bottom : str
                Boundary condition for the bottom edge.

            Raises
            ------
            ValueError
                If boundary conditions are not recognized.
        )pbdoc"
    )
    .def("_cpp_initialize",
        &EigenSolver::initialize,
        pybind11::arg("model_parameters"),
        pybind11::arg("finit_matrix"),
        pybind11::arg("n_computed_mode"),
        pybind11::arg("n_sorted_mode"),
        R"pbdoc(
            Initialize the eigenvalue solver with the provided parameters.

            This method sets up the solver with the necessary parameters and
            prepares it for mode computation. It should be called before any
            eigenvalue solving methods.

            Parameters
            ----------
            model_parameters : ModelParameters
                Physical and geometric parameters of the waveguide model.
            finit_matrix : numpy.ndarray
                Finite difference matrix representing the discretized wave equation.
            n_computed_mode : int
                Number of eigenmodes to compute.
            n_sorted_mode : int

            Raises
            ------
            ValueError
                If boundary conditions are not recognized or parameters are invalid.
        )pbdoc"
    )
    .def(
        "loop_over_itr",
        &EigenSolver::loop_over_itr,
        pybind11::arg("alpha"),
        R"pbdoc(
            Perform iterative eigenvalue computation with extrapolation.

            This method implements an iterative algorithm to solve the eigenvalue
            problem with Richardson extrapolation for improved convergence and
            accuracy. It's particularly useful for challenging problems where
            direct methods might struggle.

            Parameters
            ----------
            alpha : float
                Extrapolation parameter controlling the step size in the
                iterative process. Optimal values depend on the problem
                characteristics. Typical range: 0.1-2.0.

            Returns
            -------
            bool
                True if convergence was achieved within max_iter iterations,
                False otherwise.

            Notes
            -----
            This method modifies the internal state of the solver and should
            be called after initialization but before retrieving modes.
            The extrapolation helps accelerate convergence especially for
            problems with slowly converging eigenvalues.

            Examples
            --------
            >>> converged = solver.loop_over_itr(extrapolation_order=2, alpha=1.0)
            >>> if not converged:
            ...     print("Warning: Solver did not converge")
        )pbdoc"
    )
    .def(
        "_cpp_compute_laplacian",
        &EigenSolver::compute_laplacian,
        R"pbdoc(
            Compute the Laplacian operator matrix for the finite difference grid.

            This method constructs the discrete Laplacian operator used in the
            wave equation eigenvalue problem. The Laplacian is computed using
            finite difference approximations on the computational grid.

            The method handles:
            - Grid spacing variations
            - Boundary condition implementation
            - Coordinate system transformations (cylindrical)

            Returns
            -------
            None
                The computed Laplacian is stored internally and used by the
                eigenvalue solver.

            Notes
            -----
            This method should be called before solving for eigenvalues.
            The Laplacian matrix is typically sparse and stored efficiently
            using appropriate data structures.

            The finite difference stencil used depends on the accuracy order
            specified in the model parameters, with higher orders providing
            better accuracy at the cost of increased computational complexity.
        )pbdoc"
    )
    .def("_cpp_reset_solver",
        &EigenSolver::reset_solver,
        R"pbdoc(
            Reset the eigenvalue solver to its initial state.

            This method clears all internal data structures and resets the
            solver's state, allowing it to be reused for a new computation
            without needing to reinitialize.

            Returns
            -------
            None
                The solver is reset and ready for a new computation.

            Notes
            -----
            This method is useful when you want to run multiple computations
            with different parameters without creating a new solver instance.
        )pbdoc"
    )
    .def(
        "_cpp_get_mode",
        &EigenSolver::get_sorted_mode,
        R"pbdoc(
            Retrieve computed and sorted eigenmodes.

            Returns the computed eigenmodes sorted by their effective indices
            (propagation constants). Only the first n_sorted_mode modes are
            returned, as specified during initialization.

            Returns
            -------
            tuple
                A tuple containing:
                - eigenvalues (numpy.ndarray): Complex propagation constants
                  of the modes, sorted by decreasing real part (effective index)
                - eigenvectors (numpy.ndarray): Mode field distributions as
                  complex arrays with shape (n_grid_points, n_sorted_mode)

            Notes
            -----
            This method should be called after the eigenvalue computation
            has been performed (either through loop_over_itr or direct solving).

            The eigenvectors represent the electric or magnetic field components
            depending on the formulation used. For vector modes, both
            components may be included.

            The sorting is based on the real part of the eigenvalue, which
            corresponds to the effective refractive index of the mode.

            Examples
            --------
            >>> eigenvals, eigenvecs = solver.get_mode()
            >>> n_eff = numpy.sqrt(eigenvals.real)  # Effective indices
            >>> fundamental_mode = eigenvecs[:, 0]   # Fundamental mode field
        )pbdoc"
    )
    .def_readwrite(
        "n_sorted_mode",
        &EigenSolver::n_sorted_mode,
        R"pbdoc(
            Number of modes to keep after sorting by effective index.

            This attribute controls how many of the computed eigenmodes are
            retained and returned by get_mode(). The modes are sorted by
            their effective refractive index (real part of eigenvalue) in
            descending order.

            Type
            ----
            int
                Must be <= n_computed_mode. Typical values range from 1
                (fundamental mode only) to 10-20 for multimode analysis.

            Notes
            -----
            Changing this value after initialization will affect subsequent
            calls to get_mode() but does not trigger recomputation of
            eigenvalues. The sorting and selection is performed on the
            already computed modes.
        )pbdoc"
    )
    .def_readwrite(
        "n_computed_mode",
        &EigenSolver::n_computed_mode,
        R"pbdoc(
            Number of eigenmodes to compute during eigenvalue solving.

            This attribute determines how many eigenvalues and eigenvectors
            the solver will compute. A larger number provides more modes
            but increases computational cost and memory usage.

            Type
            ----
            int
                Typical range: 10-100 depending on the application.
                For single-mode analysis, 5-10 is often sufficient.
                For multimode analysis, 20-100 may be needed.

            Notes
            -----
            Changing this value after initialization requires rerunning
            the eigenvalue computation to take effect. The computational
            cost scales roughly linearly with this parameter.

            It's recommended to compute more modes than needed (n_computed_mode
            > n_sorted_mode) to ensure that all physically relevant modes
            are captured, especially in the presence of numerical noise.
        )pbdoc"
    )
    .def_readwrite(
        "alpha_vector",
        &EigenSolver::alpha_vector,
        R"pbdoc(
            Vector of alpha parameters used in iterative extrapolation.

            This attribute stores the sequence of alpha values used during
            Richardson extrapolation in the iterative solver. It provides
            insight into the convergence behavior and can be used for
            diagnostic purposes.

            Type
            ----
            numpy.ndarray
                Array of float values representing the alpha parameters
                used in each iteration step.

            Notes
            -----
            This vector is populated during the execution of loop_over_itr()
            and can be examined to understand the solver's convergence
            characteristics. Large variations in alpha values may indicate
            convergence difficulties.

            The optimal alpha value often lies between 0.5 and 2.0, but
            this can vary significantly depending on the problem characteristics
            and the condition number of the eigenvalue problem.
        )pbdoc"
    )

    .def_readonly(
        "boundaries",
        &EigenSolver::boundaries,
        R"pbdoc(
            Boundary condition applied to the left edge of the computational domain.

            This read-only attribute specifies the type of boundary condition
            used at the left boundary of the finite difference grid. The
            boundary condition affects how the wave equation is discretized
            near the domain edges.

            Type
            ----
            str
                One of "PML", "PEC", or "PMC"

            Options
            -------
            - "PML": Perfectly Matched Layer - absorbing boundary that
              minimizes reflections, recommended for most applications
            - "PEC": Perfect Electric Conductor - electric field tangential
              component is zero at the boundary
            - "PMC": Perfect Magnetic Conductor - magnetic field tangential
              component is zero at the boundary
        )pbdoc"
    )
    ;
}
