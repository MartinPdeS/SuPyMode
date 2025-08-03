#include <pybind11/pybind11.h>
#include "supermode.h"

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

    pybind11::class_<SuperMode>(module, "SUPERMODE",
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

            Examples
            --------
            >>> # Access supermode properties
            >>> fields = supermode.get_fields()
            >>> betas = supermode.get_betas()
            >>>
            >>> # Analyze coupling between modes
            >>> coupling = supermode.get_normalized_coupling_with_mode(other_mode)
            >>> beating_length = supermode.get_beating_length_with_mode(other_mode)
        )pbdoc")
    .def_readonly("fields", &SuperMode::fields,
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
        )pbdoc")
    .def_readwrite("binding_number", &SuperMode::mode_number,
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
        )pbdoc")
    .def_readwrite("model_parameters", &SuperMode::model_parameters,
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
        )pbdoc")

    .def("itr_list", &SuperMode::get_itr_list,
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

            Notes
            -----
            This information is useful for diagnosing convergence issues
            and understanding the computational cost of the supermode
            calculation. It can help optimize solver parameters for
            better performance or accuracy.
        )pbdoc")

    .def("get_fields", &SuperMode::get_fields_py,
        R"pbdoc(
            Retrieve all supermode field distributions.

            Returns the complete set of electric field distributions for all
            computed supermodes in the coupled waveguide system. Each supermode
            represents a normal mode of the entire coupled structure.

            Returns
            -------
            numpy.ndarray
                Complex array containing electric field distributions with
                shape (n_grid_points, n_supermodes). Each column represents
                the spatial field distribution of one supermode.

            Notes
            -----
            The returned fields are typically normalized to unit power or
            energy. The complex values include both amplitude and phase
            information, which is essential for understanding interference
            and beating effects in coupled systems.

            Examples
            --------
            >>> fields = supermode.get_fields()
            >>> fundamental_supermode = fields[:, 0]  # First supermode
            >>> field_intensity = numpy.abs(fields)**2
        )pbdoc")
    .def("get_field", &SuperMode::get_field_py, pybind11::arg("index"),
        R"pbdoc(
            Retrieve the field distribution of a specific supermode.

            Returns the electric field distribution for a single supermode
            identified by its index. This is more efficient than retrieving
            all fields when only one is needed.

            Parameters
            ----------
            index : int
                Index of the desired supermode (0-based). Must be less than
                the total number of computed supermodes.

            Returns
            -------
            numpy.ndarray
                Complex array containing the electric field distribution
                of the specified supermode with shape (n_grid_points,).

            Raises
            ------
            IndexError
                If the index is out of range for the computed supermodes.

            Examples
            --------
            >>> fundamental = supermode.get_field(0)  # Fundamental supermode
            >>> first_higher = supermode.get_field(1)  # First higher-order supermode
        )pbdoc")
    .def("get_norm", &SuperMode::get_norm,
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
        )pbdoc")
    .def("get_index", &SuperMode::get_index_py,
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
        )pbdoc")
    .def("get_betas", &SuperMode::get_betas_py,
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
        )pbdoc")
    .def("get_eigen_value", &SuperMode::get_eigen_value_py,
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
        )pbdoc")

    .def("get_normalized_coupling_with_mode", &SuperMode::get_normalized_coupling_with_mode_py,
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

            Examples
            --------
            >>> coupling_matrix = supermode1.get_normalized_coupling_with_mode(supermode2)
            >>> max_coupling = numpy.max(numpy.abs(coupling_matrix))
        )pbdoc")
    .def("get_adiabatic_with_mode", &SuperMode::get_adiabatic_with_mode_py,
        R"pbdoc(
            Analyze adiabatic coupling conditions with another supermode set.

            Evaluates whether the transition between this supermode configuration
            and another can be considered adiabatic. Adiabatic transitions
            preserve mode identity and minimize unwanted mode coupling.

            Parameters
            ----------
            other_mode : SuperMode
                Another SuperMode object representing a nearby configuration
                in parameter space (e.g., different position in a taper).

            Returns
            -------
            dict or numpy.ndarray
                Adiabatic coupling metrics, which may include:
                - Adiabaticity parameter values
                - Mode correlation coefficients
                - Transition probabilities
                - Local coupling strengths

            Notes
            -----
            Adiabatic evolution occurs when the rate of change of system
            parameters is much slower than the characteristic time scales
            of the system. This analysis helps determine if a taper or
            varying structure will maintain mode purity.

            The adiabatic criterion is typically expressed as:
            |dγ/dz| << |Δβ|², where γ are coupling coefficients and
            Δβ are propagation constant differences.
        )pbdoc")
    .def("get_overlap_integrals_with_mode", &SuperMode::get_overlap_integrals_with_mode_py,
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
        )pbdoc")
    .def("get_beating_length_with_mode", &SuperMode::get_beating_length_with_mode_py,
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
        )pbdoc")
    .def("is_computation_compatible", &SuperMode::is_computation_compatible,
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

            Examples
            --------
            >>> if supermode1.is_computation_compatible(supermode2):
            ...     coupling = supermode1.get_normalized_coupling_with_mode(supermode2)
            ... else:
            ...     print("Supermodes are not compatible for coupling analysis")
        )pbdoc")
    .def(pybind11::pickle(&SuperMode::get_pickle, &SuperMode::build_from_tuple),
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

            Examples
            --------
            >>> import pickle
            >>> # Save SuperMode object
            >>> with open('supermode.pkl', 'wb') as f:
            ...     pickle.dump(supermode, f)
            >>>
            >>> # Load SuperMode object
            >>> with open('supermode.pkl', 'rb') as f:
            ...     loaded_supermode = pickle.load(f)
        )pbdoc")
    ;
}

