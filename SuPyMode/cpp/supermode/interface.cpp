#include <pybind11/pybind11.h>
#include "supermode.h"

namespace py = pybind11;

/**
 * @brief Python bindings for the SuperMode C++ module
 *
 * This module provides Python bindings for SuPyMode's supermode analysis,
 * enabling computation and analysis of coupled mode systems in multi-core
 * optical fibers and waveguide arrays.
 */
PYBIND11_MODULE(interface_supermode, module)
{
    module.doc() = R"pbdoc(
        SuPyMode SuperMode Module

        This module provides advanced supermode analysis capabilities for coupled
        waveguide systems. Supermodes are the normal modes of coupled waveguide
        structures, such as multi-core fibers, directional couplers, and waveguide
        arrays.

        The module enables:
        - Computation of supermode field distributions
        - Analysis of mode coupling and power transfer
        - Calculation of beating lengths and coupling coefficients
        - Overlap integral computation between different modes
        - Adiabatic coupling analysis for tapered structures

        Supermodes are essential for understanding light propagation in coupled
        systems and designing devices like mode multiplexers, beam splitters,
        and coupled resonators.
    )pbdoc";

    pybind11::class_<SuperMode, std::shared_ptr<SuperMode>>(module, "SUPERMODE",
        R"pbdoc(
            Supermode analysis class for coupled waveguide systems.

            This class represents a collection of supermodes in a coupled waveguide
            system and provides methods for analyzing their properties and interactions.
            Supermodes are the eigenmodes of the entire coupled system, each with
            its own propagation constant and field distribution.

            In a system of N coupled waveguides, there are N supermodes, each
            representing a different way energy can propagate through the system.
            The superposition of these modes determines the overall behavior of
            light in the coupled structure.

            Key Features:
            - Field distribution access and manipulation
            - Mode coupling analysis between different supermodes
            - Beating length calculation for power oscillations
            - Overlap integral computation for mode matching
            - Adiabatic coupling analysis for varying structures

        )pbdoc"
    )
    .def(
        // pybind11::init<
        //     const size_t,
        //     const pybind11::array_t<double> &,
        //     const pybind11::array_t<double> &,
        //     const pybind11::array_t<double> &,
        //     const pybind11::array_t<double> &,
        //     const ModelParameters &,
        //     const Boundaries&
        // >(),
        pybind11::init(
            [](
                const size_t mode_number,
                const pybind11::array_t<double> &fields,
                const pybind11::array_t<double> &index,
                const pybind11::array_t<double> &betas,
                const pybind11::array_t<double> &eigen_value,
                const ModelParameters &model_parameters,
                const Boundaries& boundaries
            ) {
                Eigen::VectorXd fields_eigen = NumpyInterface::convert_py_to_eigen(fields, model_parameters.nx * model_parameters.ny, model_parameters.n_slice);
                Eigen::VectorXd index_eigen = NumpyInterface::convert_py_to_eigen(index, model_parameters.n_slice);
                Eigen::VectorXd betas_eigen = NumpyInterface::convert_py_to_eigen(betas, model_parameters.n_slice);
                Eigen::VectorXd eigen_value_eigen = NumpyInterface::convert_py_to_eigen(eigen_value, model_parameters.n_slice);

                return std::make_shared<SuperMode>(
                    mode_number,
                    fields_eigen,
                    index_eigen,
                    betas_eigen,
                    eigen_value_eigen,
                    model_parameters,
                    boundaries
                );
            }
        ),
        pybind11::arg("mode_number"),
        pybind11::arg("fields"),
        pybind11::arg("index"),
        pybind11::arg("betas"),
        pybind11::arg("eigen_value"),
        pybind11::arg("model_parameters"),
        pybind11::arg("boundaries"),
        R"pbdoc(
            Constructor for SuperMode with full field data and parameters.

            Parameters
            ----------
            mode_number : int
                Identifier for the supermode (e.g., 0 for fundamental, 1 for first higher-order mode).
            fields : numpy.ndarray
                Complex field distribution of the supermode, typically flattened from a 2D grid.
            index : numpy.ndarray
                Effective index distribution associated with the supermode.
            betas : numpy.ndarray
                Propagation constants for the supermode.
            eigen_value : numpy.ndarray
                Eigenvalues from the mode solver corresponding to the supermode.
            model_parameters : ModelParameters
                Physical and computational parameters defining the coupled system.
            boundaries : Boundaries
                Boundary conditions applied in the mode solver.
        )pbdoc"
    )
    .def_readonly(
        "fields",
        &SuperMode::fields,
        R"pbdoc(
            Electric field distributions of all supermodes.

            This attribute provides read-only access to the computed electric
            field distributions for all supermodes in the coupled system.
            The fields are stored as complex arrays representing the spatial
            distribution of each supermode.

            Type
            ----
            numpy.ndarray
                Complex array with shape (n_points, n_modes) where n_points
                is the number of grid points and n_modes is the number of
                computed supermodes.

            Notes
            -----
            The field amplitudes are typically normalized according to power
            or energy conservation. The phase information is preserved and
            is crucial for understanding interference effects in coupled systems.

            Each column represents one supermode, and the spatial distribution
            shows how that particular supermode extends across all waveguides
            in the coupled system.
        )pbdoc"
    )
    .def_readwrite(
        "binding_number",
        &SuperMode::mode_number,
        R"pbdoc(
            Binding number or mode identifier for the supermode set.

            This attribute specifies which supermode or set of supermodes
            this object represents. It's used for indexing and identification
            purposes in multi-mode calculations.

            Type
            ----
            int
                Mode index or binding number. For single supermode analysis,
                this identifies which supermode (0 = fundamental, 1 = first
                higher-order, etc.). For multi-mode analysis, it may represent
                a mode group or calculation identifier.

            Notes
            -----
            The binding number is particularly important when dealing with
            multiple supermode calculations or when tracking specific modes
            through parameter variations (e.g., in bend or taper analysis).
        )pbdoc"
    )
    .def_readwrite(
        "label",
        &SuperMode::label,
        R"pbdoc(
            Label or identifier for the supermode.
        )pbdoc"
    )
    .def_property_readonly(
        "stylized_label",
        [](const SuperMode &self) {
            return self.label.empty() ? "Mode: " + std::to_string(self.mode_number) : "$" + self.label + "$";
        },
        R"pbdoc(
            Stylized label for the supermode, combining the mode number and label.
        )pbdoc"
    )
    .def_readwrite(
        "mode_number",
        &SuperMode::mode_number,
        R"pbdoc(
            Mode number or identifier for the supermode.
        )pbdoc"
    )
    .def_readwrite(
        "solver_number",
        &SuperMode::solver_number,
        R"pbdoc(
            Solver number or identifier for the supermode.
        )pbdoc"
    )
    .def_property_readonly(
        "ID",
        [](const SuperMode &self) {
            return py::make_tuple(py::int_(self.mode_number), py::int_(self.solver_number));
        },
        R"pbdoc(
            Unique identifier for the supermode instance.
        )pbdoc"
    )
    .def_readwrite(
        "boundaries",
        &SuperMode::boundaries,
        R"pbdoc(
            Boundary conditions for the supermode.
        )pbdoc"
    )
    .def(
        "__repr__",
        [](const SuperMode &self) {
            return "<SuperMode mode_number=" + std::to_string(self.mode_number) +
                   ", solver_number=" + std::to_string(self.solver_number) +
                   ", label='" + self.label + "'>";
        },
        R"pbdoc(
            String representation of the SuperMode instance for debugging and display.

            Returns
            -------
            str
                A string representation of the SuperMode instance.
        )pbdoc"
    )
    .def_readwrite(
        "model_parameters",
        &SuperMode::model_parameters,
        R"pbdoc(
            Physical and computational parameters for the supermode model.

            This attribute contains all the parameters that define the coupled
            waveguide system, including geometric dimensions, refractive index
            profiles, wavelength, and computational settings.

            Type
            ----
            ModelParameters
                Structure containing:
                - Wavelength and frequency information
                - Refractive index distributions
                - Geometric parameters (core sizes, separations)
                - Grid resolution and boundary conditions
                - Material properties and dispersion

            Notes
            -----
            Modifying these parameters requires recomputation of the supermodes.
            The parameters directly affect the coupling strength, effective
            indices, and field distributions of the supermodes.

            Key parameters for coupled systems include:
            - Core separation distance (affects coupling strength)
            - Individual core parameters (size, index contrast)
            - Operating wavelength (affects normalized frequency)
        )pbdoc"
    )
    .def(
        "get_fields",
        [](const SuperMode &self) {
            return NumpyInterface::eigen_to_ndarray<double>(self.fields, {self.model_parameters.n_slice, self.model_parameters.ny, self.model_parameters.nx});
        },
        R"pbdoc(
            Retrieve all supermode field distributions.

            Returns the complete set of electric field distributions for all
            computed supermodes in the coupled waveguide system. Each supermode
            represents a normal mode of the entire coupled structure.

            Returns
            -------
            numpy.ndarray
                Array containing electric field distributions with
                shape (n_slice, ny, nx). Each element represents
                the spatial field distribution of one supermode.

            Notes
            -----
            The returned fields are typically normalized to unit power or
            energy. The values include both amplitude and phase
            information, which is essential for understanding interference
            and beating effects in coupled systems.

        )pbdoc"
    )
    .def(
        "itr_list",
        [](const SuperMode &self) {
            return NumpyInterface::eigen_to_ndarray<double>(self.model_parameters.itr_list, {static_cast<unsigned long>(self.model_parameters.itr_list.size())});
        },
        R"pbdoc(
            Get the list of iteration parameters used in supermode computation.

            This method returns information about the iterative process used
            to compute the supermodes, including convergence history and
            computational parameters.

            Returns
            -------
            list
                List of iteration parameters, convergence criteria, and
                computational metrics from the eigenvalue solver.

        )pbdoc"
    )
    .def(
        "get_mesh_gradient",
        [](const SuperMode &self) {
            return NumpyInterface::eigen_to_ndarray<double>(self.model_parameters.mesh_gradient, {self.model_parameters.nx * self.model_parameters.ny, 2});
        },
        R"pbdoc(
            Get the mesh gradients used in the supermode computation.
            This method returns the spatial gradients of the computational mesh,
            which are essential for accurate finite-difference calculations of the supermodes.

            Returns
            -------
            numpy.ndarray
                Array of shape (n_grid_points, 2) containing the x and y gradients
                of the computational mesh. The first column corresponds to the x-gradient,
                and the second column corresponds to the y-gradient.

        )pbdoc"
    )
    .def(
        "get_norm",
        &SuperMode::get_norm,
        R"pbdoc(
            Calculate the normalization factors for all supermodes.

            Computes the normalization constants used to ensure proper
            scaling of the supermode fields. This is important for
            power conservation and accurate coupling calculations.

            Returns
            -------
            numpy.ndarray
                Array of real normalization factors, one for each supermode.

            Notes
            -----
            The normalization typically ensures unit power flow for each
            supermode, which is essential for accurate coupling coefficient
            calculations and power transfer analysis.
        )pbdoc"
    )
    .def(
        "get_index",
        [](const SuperMode &self) {
            return NumpyInterface::eigen_to_ndarray<double>(self.index, {self.model_parameters.n_slice});
        },
        R"pbdoc(
            Get the effective refractive indices of all supermodes.

            Returns the effective refractive indices (real part of the
            propagation constants divided by the free-space wave number)
            for all computed supermodes.

            Returns
            -------
            numpy.ndarray
                Real array containing the effective indices of all supermodes,
                typically sorted in descending order.

            Notes
            -----
            The effective index determines the phase velocity of each supermode.
            In coupled systems, the difference in effective indices between
            supermodes determines the beating length and power transfer
            characteristics.

            For weakly coupled systems, the supermode effective indices
            are close to those of the individual uncoupled waveguides.
        )pbdoc"
    )
    .def(
        "get_betas",
        [](const SuperMode &self) {
            return NumpyInterface::eigen_to_ndarray<double>(self.betas, {self.model_parameters.n_slice});
        },
        R"pbdoc(
            Get the propagation constants (beta values) of all supermodes.

            Returns the complex propagation constants for all computed
            supermodes. The propagation constant determines how the phase
            and amplitude of each supermode evolve along the propagation
            direction.

            Returns
            -------
            numpy.ndarray
                Complex array containing the propagation constants (β) in
                units of rad/m. The real part determines phase evolution,
                while the imaginary part represents attenuation.

            Notes
            -----
            The propagation constant is related to the effective index by:
            β = n_eff * k₀, where k₀ is the free-space wave number.

            The difference in propagation constants between supermodes
            determines the beating length: L_beat = 2π / |β₁ - β₂|.
        )pbdoc"
    )
    .def(
        "get_eigen_value",
        [](const SuperMode &self) {
            return NumpyInterface::eigen_to_ndarray<double>(self.eigen_value, {self.model_parameters.n_slice,});
        },
        R"pbdoc(
            Get the eigenvalues from the supermode eigenvalue problem.

            Returns the raw eigenvalues obtained from solving the wave
            equation eigenvalue problem. These are related to but not
            identical to the propagation constants.

            Returns
            -------
            numpy.ndarray
                Complex array of eigenvalues from the matrix eigenvalue
                problem. The relationship to physical propagation constants
                depends on the specific formulation used.

            Notes
            -----
            The eigenvalues are the direct output of the numerical eigenvalue
            solver and may require transformation to obtain physically
            meaningful propagation constants. The transformation depends
            on the specific wave equation formulation and discretization.
        )pbdoc"
    )
    .def(
        "get_normalized_coupling_with_mode",
        [](const SuperMode &self, const SuperMode &other) {
            Eigen::Matrix<complex128, Eigen::Dynamic, Eigen::Dynamic> coupling = self.get_normalized_coupling_with_mode(other);
            return NumpyInterface::eigen_to_ndarray<complex128>(coupling, {self.model_parameters.n_slice});
        },
        R"pbdoc(
            Calculate normalized coupling coefficients with another supermode set.

            Computes the normalized coupling coefficients between the supermodes
            in this object and those in another SuperMode object. This is
            essential for analyzing power transfer and mode evolution in
            varying coupled systems.

            Parameters
            ----------
            other_mode : SuperMode
                Another SuperMode object representing a different configuration
                or position along a varying structure.

            Returns
            -------
            numpy.ndarray
                Matrix of normalized coupling coefficients with shape
                (n_modes_self, n_modes_other). Element (i,j) represents
                the coupling between mode i in this object and mode j
                in the other object.

            Notes
            -----
            Normalized coupling coefficients are dimensionless and represent
            the fraction of power that can be transferred between modes.
            Values close to 1 indicate strong coupling, while values near
            0 indicate weak coupling.

            This is particularly useful for analyzing adiabatic tapers,
            where the coupling coefficients should remain close to unity
            for adiabatic mode evolution.

        )pbdoc"
    )
    .def(
        "get_overlap_integrals_with_mode",
        [](const SuperMode &self, const SuperMode &other) {
            Eigen::MatrixXd overlaps = self.get_overlap_integrals_with_mode(other);
            return NumpyInterface::eigen_to_ndarray<double>(overlaps, {self.model_parameters.n_slice});
        },
        pybind11::arg("other_mode"),
        R"pbdoc(
            Compute overlap integrals between supermodes and another mode set.

            Calculates the spatial overlap integrals between the field
            distributions of supermodes in this object and those in another
            SuperMode object. Overlap integrals quantify the similarity
            between mode field distributions.

            Parameters
            ----------
            other_mode : SuperMode
                Another SuperMode object with which to compute overlaps.

            Returns
            -------
            numpy.ndarray
                Matrix of complex overlap integrals with shape
                (n_modes_self, n_modes_other). The overlap integral
                between modes i and j is: ∫ E₁ᵢ* · E₂ⱼ dA

            Notes
            -----
            Overlap integrals are fundamental to understanding mode coupling
            and power transfer efficiency. They appear in coupled mode theory
            and are essential for calculating coupling coefficients.

            The overlap integral is defined as:
            S₁₂ = ∫∫ E₁*(x,y) · E₂(x,y) dx dy

            For identical modes, the overlap integral equals 1 (assuming
            proper normalization). For orthogonal modes, it equals 0.

            Examples
            --------
            >>> overlaps = supermode1.get_overlap_integrals_with_mode(supermode2)
            >>> diagonal_overlaps = numpy.diag(overlaps)  # Same mode overlaps
        )pbdoc"
    )
    .def(
        "get_overlap_integrals_with_mode",
        [](const SuperMode &self, const SuperMode &other) {
            Eigen::MatrixXd overlaps = self.get_overlap_integrals_with_mode(other);
            return NumpyInterface::eigen_to_ndarray<double>(overlaps, {self.model_parameters.n_slice});
        },
        pybind11::arg("other_mode"),
        R"pbdoc(
            Compute normalized overlap integrals between supermodes and another mode set.

            Calculates the normalized spatial overlap integrals between the field
            distributions of supermodes in this object and those in another SuperMode object. Normalized overlap integrals provide a dimensionless measure of mode similarity.


            Parameters
            ----------
            other_mode : SuperMode
                Another SuperMode object with which to compute normalized overlaps.


            Returns
            -------
            numpy.ndarray
                Matrix of complex normalized overlap integrals with shape (n_modes_self, n_modes_other).
        )pbdoc"
    )
    .def(
        "get_adiabatic_with_mode",
        [](const SuperMode &self, const SuperMode &other) {
            Eigen::MatrixXd adiabaticity = self.get_adiabatic_with_mode(other);
            return NumpyInterface::eigen_to_ndarray<double>(adiabaticity, {self.model_parameters.n_slice});
        },
        pybind11::arg("other_mode"),
        R"pbdoc(
            Calculate the adiabaticity parameter between supermodes and another mode set.
            Evaluates the adiabaticity parameter, which quantifies how closely the transition between two supermode configurations follows adiabatic evolution. A low adiabaticity parameter indicates that the system is evolving slowly enough to maintain mode purity, while a high value suggests significant non-adiabatic coupling and potential mode mixing.

            Parameters
            ----------
            other_mode : SuperMode
                Another SuperMode object representing a different configuration or position along a varying structure.

            Returns
            -------
            numpy.ndarray
                Matrix of adiabaticity parameters with shape (n_modes_self, n_modes_other).
            Notes
            -----
            The adiabaticity parameter is typically defined as:
            A_ij = |dγ_ij/dz| / |Δβ_ij|²
            where γ_ij are the coupling coefficients between modes i and j, and Δβ_ij are the differences in propagation constants. A value of A_ij << 1 indicates adiabatic evolution, while A_ij ~ 1 or greater indicates non-adiabatic behavior.
        )pbdoc"
    )
    .def(
        "get_beating_length_with_mode",
        [](const SuperMode &self, const SuperMode &other) {
            Eigen::MatrixXd beating_lengths = self.get_beating_length_with_mode(other);
            return NumpyInterface::eigen_to_ndarray<double>(beating_lengths, {self.model_parameters.n_slice});
        },
        R"pbdoc(
            Calculate beating lengths between supermodes and another mode set.

            Computes the beating lengths that characterize the periodic power
            exchange between supermodes. The beating length is the distance
            over which power oscillates between coupled modes.

            Parameters
            ----------
            other_mode : SuperMode
                Another SuperMode object representing a different supermode
                configuration (typically at a different position).

            Returns
            -------
            numpy.ndarray
                Matrix of beating lengths with shape (n_modes_self, n_modes_other).
                Element (i,j) is the beating length between mode i in this
                object and mode j in the other object, in units of length.

            Notes
            -----
            The beating length is calculated from the difference in propagation
            constants: L_beat = 2π / |β₁ - β₂|

            Short beating lengths indicate strong coupling and rapid power
            oscillation, while long beating lengths indicate weak coupling.

            In practical devices:
            - Directional couplers use specific lengths related to beating length
            - Multi-core fibers have beating lengths that determine crosstalk
            - Mode multiplexers require understanding of beating characteristics

            Examples
            --------
            >>> beating_lengths = supermode1.get_beating_length_with_mode(supermode2)
            >>> min_beating = numpy.min(beating_lengths[beating_lengths > 0])
        )pbdoc"
    )
    .def(
        "is_computation_compatible",
        &SuperMode::is_computation_compatible,
        R"pbdoc(
            Check computational compatibility with another SuperMode object.

            Verifies whether this SuperMode object is computationally compatible
            with another one, meaning they can be used together in coupled
            calculations. Compatibility requires matching grid sizes, coordinate
            systems, and computational parameters.

            Parameters
            ----------
            other_mode : SuperMode
                Another SuperMode object to check compatibility with.

            Returns
            -------
            bool
                True if the objects are compatible for joint calculations,
                False otherwise.

            Notes
            -----
            Compatibility is required for:
            - Coupling coefficient calculations
            - Overlap integral computation
            - Adiabatic analysis
            - Mode evolution tracking

            Incompatibility typically arises from:
            - Different grid resolutions
            - Mismatched coordinate systems
            - Different wavelengths or frequencies
            - Inconsistent boundary conditions

        )pbdoc"
    )
    .def(
        pybind11::pickle(
            [](const SuperMode &supermode) {
                return pybind11::make_tuple(
                    supermode.mode_number,
                    NumpyInterface::eigen_to_ndarray<double>(supermode.fields, {supermode.model_parameters.nx * supermode.model_parameters.ny, supermode.model_parameters.n_slice}),
                    NumpyInterface::eigen_to_ndarray<double>(supermode.index, {supermode.model_parameters.n_slice}),
                    NumpyInterface::eigen_to_ndarray<double>(supermode.betas, {supermode.model_parameters.n_slice}),
                    NumpyInterface::eigen_to_ndarray<double>(supermode.eigen_value, {supermode.model_parameters.n_slice}),
                    supermode.model_parameters,
                    supermode.boundaries
                );
            },
            [](const pybind11::tuple &tuple) {
                return SuperMode{
                    tuple[0].cast<size_t>(),                             // mode_number
                    tuple[1].cast<Eigen::VectorXd>(),                    // fields
                    tuple[2].cast<Eigen::VectorXd>(),                    // index
                    tuple[3].cast<Eigen::VectorXd>(),                    // betas
                    tuple[4].cast<Eigen::VectorXd>(),                    // eigen_values
                    tuple[5].cast<ModelParameters>(),                    // Model parameters
                    tuple[6].cast<Boundaries>()                          // boundaries
                }; // load
            }
        ),
        R"pbdoc(
            Enable Python pickle serialization support.

            Provides serialization and deserialization capabilities for SuperMode
            objects, allowing them to be saved to files, transmitted over networks,
            or used in parallel computing environments.

            Notes
            -----
            Pickle support enables:
            - Saving SuperMode objects to disk for later analysis
            - Parallel processing with multiprocessing module
            - Caching expensive computations
            - Inter-process communication in distributed computing

            The pickle protocol preserves all computed fields, parameters,
            and internal state, allowing full reconstruction of the object.

        )pbdoc"
    )
    .def_property_readonly(
        "beta",
        [](const SuperMode &self) {
            py::object Beta = py::module_::import("SuPyMode.representation").attr("Beta");
            return Beta(self);
        },
        R"pbdoc(
            Property to access the Beta representation (propagation constants) of the supermode.

            This property allows users to access the Beta representation, which
            provides a convenient interface for analyzing the propagation constants
            and related properties of the supermodes.

            Returns
            -------
            Beta
                An instance of the Beta class representing the propagation constants
                and related properties of the supermodes.

            Notes
            -----
            The Beta representation is useful for:
            - Analyzing effective indices and phase velocities
            - Calculating beating lengths and coupling characteristics
            - Visualizing mode dispersion and propagation behavior
        )pbdoc"
    )
    .def_property_readonly(
        "field",
        [](const SuperMode &self) {
            py::object Field = py::module_::import("SuPyMode.representation").attr("Field");
            return Field(self);
        },
        R"pbdoc(
            Property to access the Field representation of the supermode.

            This property allows users to access the Field representation, which
            provides a convenient interface for analyzing the electric field
            distributions and related properties of the supermodes.

            Returns
            -------
            Field
                An instance of the Field class representing the electric field
                distributions and related properties of the supermodes.

            Notes
            -----
            The Field representation is useful for:
            - Visualizing mode profiles and spatial distributions
            - Analyzing mode overlap and coupling characteristics
            - Calculating power flow and intensity distributions
        )pbdoc"
    )
    .def_property_readonly(
        "adiabatic",
        [](const SuperMode &self) {
            py::object Adiabatic = py::module_::import("SuPyMode.representation").attr("Adiabatic");
            return Adiabatic(self);
        },
        R"pbdoc(
            Property to access the Adiabatic representation of the supermode.

            This property allows users to access the Adiabatic representation, which
            provides a convenient interface for analyzing adiabatic coupling and mode evolution characteristics of the supermodes.

            Returns
            -------
            Adiabatic
                An instance of the Adiabatic class representing the adiabatic coupling and mode evolution properties of the supermodes.

            Notes
            -----
            The Adiabatic representation is useful for:
            - Analyzing adiabatic coupling conditions and mode evolution
            - Evaluating transition probabilities in varying structures
            - Designing adiabatic tapers and mode converters
        )pbdoc"
    )
    .def_property_readonly(
        "normalized_coupling",
        [](const SuperMode &self) {
            py::object NormalizedCoupling = py::module_::import("SuPyMode.representation").attr("NormalizedCoupling");
            return NormalizedCoupling(self);
        },
        R"pbdoc(
            Property to access the Normalized Coupling representation of the supermode.

            This property allows users to access the Normalized Coupling representation, which
            provides a convenient interface for analyzing the normalized coupling characteristics
            of the supermodes.

            Returns
            -------
            NormalizedCoupling
                An instance of the NormalizedCoupling class representing the normalized coupling
                characteristics of the supermodes.

            Notes
            -----
            The Normalized Coupling representation is useful for:
            - Analyzing mode overlap and coupling efficiency
            - Evaluating power transfer between modes
            - Designing mode converters and couplers
        )pbdoc"
    )
    .def_property_readonly(
        "index",
        [](const SuperMode &self) {
            py::object Index = py::module_::import("SuPyMode.representation").attr("Index");
            return Index(self);
        },
        R"pbdoc(
            Property to access the Index representation of the supermode.

            This property allows users to access the Index representation, which
            provides a convenient interface for analyzing the effective refractive indices
            and related properties of the supermodes.

            Returns
            -------
            Index
                An instance of the Index class representing the effective refractive indices
                and related properties of the supermodes.

            Notes
            -----
            The Index representation is useful for:
            - Analyzing effective indices and phase velocities
            - Calculating beating lengths and coupling characteristics
            - Visualizing mode dispersion and propagation behavior
        )pbdoc"
    )
    .def_property_readonly(
        "eigenvalue",
        [](const SuperMode &self) {
            py::object EigenValue = py::module_::import("SuPyMode.representation").attr("EigenValue");
            return EigenValue(self);
        },
        R"pbdoc(
            Property to access the EigenValue representation of the supermode.

            This property allows users to access the EigenValue representation, which
            provides a convenient interface for analyzing the raw eigenvalues obtained
            from the supermode eigenvalue problem.

            Returns
            -------
            EigenValue
                An instance of the EigenValue class representing the raw eigenvalues
                from the supermode eigenvalue problem.

            Notes
            -----
            The EigenValue representation is useful for:
            - Analyzing the raw output of the numerical eigenvalue solver
            - Understanding the relationship between eigenvalues and physical propagation constants
            - Diagnosing convergence and computational issues in the eigenvalue problem
        )pbdoc"
    )
    .def_property_readonly(
        "beating_length",
        [](const SuperMode &self) {
            py::object BeatingLength = py::module_::import("SuPyMode.representation").attr("BeatingLength");
            return BeatingLength(self);
        },
        R"pbdoc(
            Property to access the Beating Length representation of the supermode.
            This property allows users to access the Beating Length representation, which
            provides a convenient interface for analyzing the beating lengths that characterize the periodic power exchange between supermodes.

            Returns
            -------
            BeatingLength
                An instance of the BeatingLength class representing the beating lengths between supermodes.

            Notes
            -----
            The Beating Length representation is useful for:
            - Analyzing mode coupling and power transfer characteristics
            - Designing directional couplers and multi-core fibers
            - Understanding interference effects in coupled systems
        )pbdoc"
    )
    .def(
        "__hash__",
        [](const SuperMode &self) {
            // Create a hash based on the unique identifier of the supermode
            std::size_t h1 = std::hash<int>()(self.mode_number);
            std::size_t h2 = std::hash<int>()(self.solver_number);
            return h1 ^ (h2 << 1); // Combine hashes
        },
        R"pbdoc(
            Compute a hash value for the SuperMode instance.
            This allows SuperMode objects to be used in hash-based collections like sets and dictionaries.
        )pbdoc"
    )
    .def(
        "get_field_interpolation",
        [](const SuperMode &self, double itr) {
            Eigen::VectorXd interpolated_field = self.get_field_interpolation(itr);

            return NumpyInterface::eigen_to_ndarray<double>(interpolated_field, {self.model_parameters.ny, self.model_parameters.nx});
        },
        pybind11::arg("itr"),
        R"pbdoc(
            Get interpolated field data for a specific iteration.

            This method provides access to the interpolated field data for a given iteration parameter, allowing users to retrieve the mode field distribution at specific points in the iteration process.
            Parameters
            ----------
            itr : float
                The iteration value for which to retrieve the interpolated field data.

            Returns
            -------
            numpy.ndarray
                A 2D array of shape (ny, nx) containing the interpolated field data

            Notes
            -----
            This method is useful for:
            - Visualizing mode field evolution during iterative solvers
            - Analyzing convergence behavior of mode calculations
            - Extracting intermediate field distributions for diagnostics

        )pbdoc"
    )
    ;
}
