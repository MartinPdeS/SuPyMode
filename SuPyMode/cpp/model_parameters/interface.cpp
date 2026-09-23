#include <pybind11/pybind11.h>
#include "model_parameters.h"

/**
 * @brief Python bindings for the ModelParameters C++ module
 *
 * This module provides Python bindings for SuPyMode's model parameter
 * management, encapsulating all physical and computational parameters
 * needed for waveguide mode analysis.
 */
PYBIND11_MODULE(interface_model_parameters, module) {
    module.doc() = R"pbdoc(
        SuPyMode ModelParameters Module

        This module provides a comprehensive parameter management system for
        optical waveguide and fiber mode calculations. The ModelParameters
        class encapsulates all physical, geometric, and computational parameters
        required for accurate mode analysis.

        Key functionalities:
        - Physical parameter management (wavelength, refractive index)
        - Geometric parameter handling (grid spacing, coordinates)
        - Computational grid setup and optimization
        - Refractive index mesh management and gradient computation
        - Parameter validation and consistency checking

        The ModelParameters object serves as the foundation for all SuPyMode
        calculations, ensuring consistent parameter usage across different
        computational modules and providing efficient access to derived
        quantities like wave numbers and normalized parameters.
    )pbdoc";

    pybind11::class_<ModelParameters>(module, "ModelParameters",
        R"pbdoc(
            Central parameter management class for waveguide mode analysis.

            This class manages all parameters required for optical waveguide
            mode calculations, including physical properties (wavelength,
            refractive index), geometric specifications (grid dimensions,
            spacing), and computational settings (iteration lists, debug modes).

            The ModelParameters object provides:
            - Centralized parameter storage and validation
            - Automatic computation of derived quantities
            - Efficient mesh and gradient management
            - Consistent coordinate system handling
            - Parameter serialization for data persistence

            All SuPyMode computational modules rely on ModelParameters to
            ensure consistent and accurate calculations across different
            analysis types (single mode, supermode, coupling analysis).

            Key Features:
            - Physical parameter management and validation
            - Computational grid setup and optimization
            - Refractive index mesh handling with gradient computation
            - Wavelength and frequency parameter conversion
            - Cross-module parameter consistency

        )pbdoc")
        .def(
            pybind11::init<
                const pybind11::array_t<double>&,
                const pybind11::array_t<double>&,
                const pybind11::array_t<double>&,
                const double,
                const double,
                const int
            >(),
            pybind11::arg("mesh"),
            pybind11::arg("x_vector"),
            pybind11::arg("y_vector"),
            pybind11::arg("dx"),
            pybind11::arg("dy"),
            pybind11::arg("debug_mode") = 0,
            R"pbdoc(
                Initialize ModelParameters with physical and computational settings.

                Creates a comprehensive parameter object containing all information
                needed for waveguide mode analysis, including the refractive index
                distribution, computational grid, and physical properties.

                Parameters
                ----------
                mesh : numpy.ndarray
                    2D refractive index distribution with shape (ny, nx).
                    Values should be real refractive indices (e.g., 1.0-1.5
                    for typical optical materials). The mesh defines the
                    waveguide structure and material properties.
                x_vector : numpy.ndarray
                    1D array of x-coordinates in meters, length nx.
                    Defines the horizontal spatial grid. Should be uniformly
                    spaced and centered around the waveguide structure.
                y_vector : numpy.ndarray
                    1D array of y-coordinates in meters, length ny.
                    Defines the vertical spatial grid. Should be uniformly
                    spaced and centered around the waveguide structure.
                itr_list : numpy.ndarray
                    1D array of iteration parameters for adaptive computations.
                    Used by iterative solvers for convergence control and
                    adaptive mesh refinement. Typical values: [1.0, 1.5, 2.0].
                wavelength : float
                    Operating wavelength in meters (e.g., 1550e-9 for 1550 nm).
                    This determines the wave number and normalized frequency
                    of the waveguide system.
                dx : float
                    Grid spacing in x-direction in meters. Should match the
                    spacing of x_vector: dx = x_vector[1] - x_vector[0].
                dy : float
                    Grid spacing in y-direction in meters. Should match the
                    spacing of y_vector: dy = y_vector[1] - y_vector[0].
                debug_mode : int, optional
                    Debug level for computational diagnostics (default: 0).
                    - 0: No debug output
                    - 1: Basic convergence information
                    - 2: Detailed computational diagnostics
                    - 3: Full debug output with intermediate results

                Raises
                ------
                ValueError
                    If array dimensions are incompatible or grid spacing
                    doesn't match coordinate vectors.
                RuntimeError
                    If memory allocation fails for large meshes.

                Notes
                -----
                Grid Requirements:
                - The mesh should extend sufficiently beyond the waveguide
                  core to capture evanescent fields (typically 5-10 core radii)
                - Grid spacing should satisfy the Nyquist criterion for the
                  highest spatial frequencies in the mode fields
                - Uniform spacing is assumed for finite difference accuracy

                Physical Considerations:
                - Wavelength should be in vacuum units
                - Refractive indices should be real for lossless analysis
                - The coordinate system assumes z as the propagation direction

            )pbdoc"
        )
        .def(
            "initialize",
            &ModelParameters::initialize,
            pybind11::arg("wavelength"),
            pybind11::arg("itr_initial"),
            pybind11::arg("itr_final"),
            pybind11::arg("n_step"),
            R"pbdoc(
                Initialize derived parameters based on physical settings.

                This method computes all derived parameters from the initial
                physical settings, including the wave number, iteration list,
                and scaled parameters for numerical stability.

                Parameters
                ----------
                wavelength : float
                    Operating wavelength in meters (e.g., 1550e-9 for 1550 nm).
                itr_initial : float
                    Initial iteration parameter for adaptive computations.
                itr_final : float
                    Final iteration parameter for adaptive computations.
                n_step : int
                    Number of steps in the iteration list (>= 2).

            )pbdoc"
        )
        .def_readonly(
            "n_slice",
            &ModelParameters::n_slice,
            R"pbdoc(
                Number of slices in the computational domain.

                This read-only attribute indicates the number of discrete slices
                or layers used in the computational analysis. For 2D problems,
                this is typically 1, while for 3D or multi-layer structures,
                it represents the number of z-direction discretization points.

                Type
                ----
                int
                    Number of computational slices (>= 1)

                Notes
                -----
                The slice count affects memory usage and computational complexity.
                For single-layer waveguides, n_slice = 1 is typical. For
                multi-layer structures or 3D analysis, larger values are used.
            )pbdoc"
        )
        .def_readonly(
            "nx",
            &ModelParameters::nx,
            R"pbdoc(
                Number of grid points in the x-direction.

                This read-only attribute specifies the horizontal resolution
                of the computational grid. It corresponds to the length of
                the x_vector array and the second dimension of the mesh array.

                Type
                ----
                int
                    Number of x-direction grid points

                Notes
                -----
                Higher nx values provide better spatial resolution but increase
                computational cost. Typical values range from 101 to 1001
                depending on the required accuracy and available resources.

                The grid should extend sufficiently beyond the waveguide core
                to capture evanescent field decay.
            )pbdoc"
        )
        .def_readonly(
            "ny",
            &ModelParameters::ny,
            R"pbdoc(
                Number of grid points in the y-direction.

                This read-only attribute specifies the vertical resolution
                of the computational grid. It corresponds to the length of
                the y_vector array and the first dimension of the mesh array.

                Type
                ----
                int
                    Number of y-direction grid points

                Notes
                -----
                For circular fibers, nx and ny are typically equal to maintain
                symmetry. For rectangular waveguides, they may differ based
                on the aspect ratio and required resolution.
            )pbdoc"
        )
        .def_readonly(
            "dx",
            &ModelParameters::dx,
            R"pbdoc(
                Grid spacing in the x-direction.

                This read-only attribute specifies the uniform spacing between
                adjacent grid points in the horizontal direction, in meters.

                Type
                ----
                float
                    Grid spacing in meters

                Notes
                -----
                The grid spacing determines the spatial resolution and affects
                the accuracy of finite difference approximations. It should be
                small enough to resolve the smallest features in the mode field,
                typically much smaller than the wavelength in the medium.

                Typical values: 10-100 nm for optical fibers at telecom wavelengths.
            )pbdoc")
        .def_readonly("dy", &ModelParameters::dy,
            R"pbdoc(
                Grid spacing in the y-direction.

                This read-only attribute specifies the uniform spacing between
                adjacent grid points in the vertical direction, in meters.

                Type
                ----
                float
                    Grid spacing in meters

                Notes
                -----
                For optimal finite difference accuracy, dx and dy should be
                approximately equal. The spacing must satisfy stability and
                accuracy requirements for the numerical scheme used.
            )pbdoc")
        .def_readwrite("wavelength", &ModelParameters::wavelength,
            R"pbdoc(
                Operating wavelength in vacuum.

                This read-only attribute stores the vacuum wavelength of the
                electromagnetic radiation being analyzed, in meters.

                Type
                ----
                float
                    Wavelength in meters (e.g., 1550e-9 for 1550 nm)

                Notes
                -----
                The wavelength determines:
                - Wave number k₀ = 2π/λ
                - Normalized frequency V = k₀ * a * NA
                - Mode cutoff conditions
                - Dispersion characteristics

                Common telecom wavelengths: 1310 nm, 1550 nm
            )pbdoc")
        .def_readonly("wavenumber", &ModelParameters::wavenumber,
            R"pbdoc(
                Free-space wave number.

                This read-only attribute provides the wave number in vacuum,
                calculated as k₀ = 2π/λ, where λ is the wavelength.

                Type
                ----
                float
                    Wave number in rad/m

                Notes
                -----
                The wave number is fundamental to all electromagnetic calculations:
                - Propagation constants: β = n_eff * k₀
                - Phase relationships: φ = β * z
                - Normalized parameters: V = k₀ * a * NA

                It's automatically computed from the wavelength parameter.
            )pbdoc"
        )
        .def_property_readonly(
            "mesh",
            [](const ModelParameters &self) {
                pybind11::array_t<double> mesh_array = NumpyInterface::eigen_to_ndarray<double>(self.mesh, {static_cast<unsigned long>(self.ny), static_cast<unsigned long>(self.nx)});
                return mesh_array;
            },
            R"pbdoc(
                Refractive index distribution mesh.

                This read-only attribute provides access to the 2D refractive
                index distribution that defines the waveguide structure and
                material properties across the computational domain.

                Type
                ----
                numpy.ndarray
                    2D array with shape (ny, nx) containing real refractive
                    index values

                Notes
                -----
                The mesh defines:
                - Core and cladding regions
                - Material boundaries and interfaces
                - Index contrast and numerical aperture
                - Waveguide geometry and structure

                Values typically range from 1.0 (air) to 1.5 (glass) for
                optical fibers, with core indices slightly higher than cladding.

                The mesh should be smooth to avoid numerical artifacts, with
                gradual transitions at material interfaces when possible.
            )pbdoc"
        )
        .def_property_readonly(
            "itr_list",
            [](const ModelParameters &self) {
                return NumpyInterface::eigen_to_ndarray<double>(self.itr_list, {static_cast<unsigned long>(self.itr_list.size())});
            },
            R"pbdoc(
                Iteration parameter list for adaptive computations.

                This read-only attribute provides access to the iteration
                parameters used by adaptive solvers and mesh refinement
                algorithms. These parameters control convergence behavior
                and computational accuracy.

                Type
                ----
                numpy.ndarray
                    1D array of iteration parameters

                Notes
                -----
                The iteration list is used by:
                - Adaptive eigenvalue solvers
                - Mesh refinement algorithms
                - Convergence acceleration schemes
                - Multi-grid methods

                Typical values are monotonically increasing: [1.0, 1.5, 2.0]
            )pbdoc"
        )
        .def_readonly(
            "x_vector",
            &ModelParameters::x_vector_py,
            R"pbdoc(
                X-coordinate vector for the computational grid.

                This read-only attribute provides the 1D array of x-coordinates
                corresponding to the horizontal grid points in the computational
                domain.

                Type
                ----
                numpy.ndarray
                    1D array of x-coordinates in meters, length nx

                Notes
                -----
                The x_vector defines the spatial grid in the horizontal direction.
                It should be uniformly spaced and centered around the waveguide
                structure for optimal accuracy.

                The spacing between points should match the dx parameter:
                dx = x_vector[1] - x_vector[0]
            )pbdoc"
        )
        .def_readonly(
            "y_vector",
            &ModelParameters::y_vector_py,
            R"pbdoc(
                Y-coordinate vector for the computational grid.

                This read-only attribute provides the 1D array of y-coordinates
                corresponding to the vertical grid points in the computational
                domain.

                Type
                ----
                numpy.ndarray
                    1D array of y-coordinates in meters, length ny

                Notes
                -----
                The y_vector defines the spatial grid in the vertical direction.
                It should be uniformly spaced and centered around the waveguide
                structure for optimal accuracy.

                The spacing between points should match the dy parameter:
                dy = y_vector[1] - y_vector[0]
            )pbdoc"
        )
        .def_readonly("mesh_gradient", &ModelParameters::mesh_gradient_py,
            R"pbdoc(
                Gradient of the refractive index mesh.

                This read-only attribute provides the computed gradient of the
                refractive index distribution, which is essential for accurate
                mode calculations in graded-index structures and for handling
                material interfaces properly.

                Type
                ----
                tuple of numpy.ndarray
                    Tuple containing (gradient_x, gradient_y) arrays, each
                    with shape (ny, nx)

                Notes
                -----
                The gradient information is used for:
                - Graded-index fiber analysis
                - Accurate boundary condition implementation
                - Interface treatment in finite difference schemes
                - Dispersion and nonlinearity calculations

                The gradient is typically computed using finite difference
                approximations with appropriate boundary handling. Sharp
                interfaces may require special treatment to maintain accuracy.

                For step-index structures, the gradient has non-zero values
                only at the core-cladding interface.
            )pbdoc")
        .def(
            pybind11::pickle(
                &ModelParameters::get_pickle,
                &ModelParameters::build_from_tuple
            ),
            R"pbdoc(
                Enable Python pickle serialization support.

                Provides serialization and deserialization capabilities for
                ModelParameters objects, allowing them to be saved to files,
                transmitted over networks, or used in parallel computing
                environments while preserving all parameter data.

                Notes
                -----
                Pickle support enables:
                - Saving parameter configurations to disk for reproducibility
                - Parallel processing with multiprocessing module
                - Caching expensive mesh computations
                - Parameter sharing between different analysis sessions
                - Version control of computational configurations

                The pickle protocol preserves all mesh data, coordinate arrays,
                physical parameters, and computational settings, allowing
                complete reconstruction of the ModelParameters object.

                Large meshes may result in significant file sizes when pickled.
                Consider using compressed formats for storage efficiency.

                Examples
                --------
                >>> import pickle
                >>> # Save ModelParameters object
                >>> with open('fiber_params.pkl', 'wb') as f:
                ...     pickle.dump(model_params, f)
                >>>
                >>> # Load ModelParameters object
                >>> with open('fiber_params.pkl', 'rb') as f:
                ...     loaded_params = pickle.load(f)
                >>>
                >>> # Verify parameters are identical
                >>> assert numpy.allclose(loaded_params.mesh, model_params.mesh)
            )pbdoc"
        )
        ;
}
