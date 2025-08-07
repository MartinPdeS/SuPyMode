#pragma once

#include <complex>
#include "../utils/utils.h"
#include "../utils/numpy_interface.h"
#include "../model_parameters/model_parameters.h"
#include <Eigen/Core>
#include <pybind11/pybind11.h>

typedef std::complex<double> complex128;

/**
 * @brief SuperMode class for optical waveguide mode analysis and computations.
 *
 * This class represents an optical supermode in a waveguide system, containing
 * the electromagnetic field distributions, propagation constants, and various
 * methods for mode analysis including coupling calculations, overlap integrals,
 * and normalization operations.
 *
 * The SuperMode class supports:
 * - Field normalization and arrangement
 * - Overlap integral calculations between modes
 * - Coupling coefficient computations
 * - Beating length and adiabatic parameter calculations
 * - Python interface through pybind11
 */
class SuperMode {
public:
    size_t mode_number;          ///< Index/number identifier for this mode
    Eigen::MatrixXd fields;      ///< Electromagnetic field distribution matrix
    Eigen::VectorXd index;       ///< Refractive index profile along propagation direction
    Eigen::VectorXd betas;       ///< Propagation constants (beta values) along the mode
    Eigen::VectorXd eigen_value; ///< Eigenvalues from the mode solver
    ModelParameters model_parameters; ///< Physical and numerical parameters for the model

    std::string left_boundary;   ///< Boundary condition specification for left edge
    std::string right_boundary;  ///< Boundary condition specification for right edge
    std::string top_boundary;    ///< Boundary condition specification for top edge
    std::string bottom_boundary; ///< Boundary condition specification for bottom edge

    /**
     * @brief Default constructor for SuperMode.
     */
    SuperMode() = default;
    /**
     * @brief Constructor with full field data and parameters.
     *
     * @param mode_number Index/identifier for this mode
     * @param fields_py Python array containing electromagnetic field data
     * @param index_py Python array containing refractive index profile
     * @param betas_py Python array containing propagation constants
     * @param eigen_value_py Python array containing eigenvalues
     * @param model_parameters Physical and numerical model parameters
     * @param left_boundary Boundary condition for left edge
     * @param right_boundary Boundary condition for right edge
     * @param top_boundary Boundary condition for top edge
     * @param bottom_boundary Boundary condition for bottom edge
     */
    SuperMode(
        const size_t mode_number,
        const pybind11::array_t<double> &fields_py,
        const pybind11::array_t<double> &index_py,
        const pybind11::array_t<double> &betas_py,
        const pybind11::array_t<double> &eigen_value_py,
        const ModelParameters &model_parameters,
        const std::string &left_boundary,
        const std::string &right_boundary,
        const std::string &top_boundary,
        const std::string &bottom_boundary
    );

    /**
     * @brief Constructor with minimal parameters (fields computed separately).
     *
     * @param mode_number Index/identifier for this mode
     * @param model_parameters Physical and numerical model parameters
     * @param left_boundary Boundary condition for left edge
     * @param right_boundary Boundary condition for right edge
     * @param top_boundary Boundary condition for top edge
     * @param bottom_boundary Boundary condition for bottom edge
     */
    SuperMode(
        const size_t mode_number,
        const ModelParameters &model_parameters,
        const std::string& left_boundary,
        const std::string& right_boundary,
        const std::string& top_boundary,
        const std::string& bottom_boundary
    );

    /**
     * @brief Calculate field normalization for a specific slice.
     *
     * @param slice The slice index along the propagation direction
     * @param normalization_type Type of normalization ("l2", "max", "cmt", etc.)
     * @return Normalization factor for the specified slice
     */
    double get_norm(const size_t &slice, const std::string &normalization_type) const;

    /**
     * @brief Get normalization array for all slices.
     *
     * @param normalization_type Type of normalization to apply
     * @return Vector containing normalization factors for all slices
     */
    Eigen::VectorXd get_norm_array(const std::string &normalization_type) const;

    /**
     * @brief Calculate CMT (Coupled Mode Theory) normalization for a slice.
     *
     * @param slice The slice index
     * @return CMT normalization factor
     */
    double get_norm_cmt(const size_t &slice) const;

    /**
     * @brief Calculate maximum field normalization for a slice.
     *
     * @param slice The slice index
     * @return Maximum field normalization factor
     */
    double get_norm_max(const size_t &slice) const;

    /**
     * @brief Calculate L2 normalization for a slice.
     *
     * @param slice The slice index
     * @return L2 normalization factor
     */
    double get_norm_l2(const size_t &slice) const;

    /**
     * @brief Get CMT normalization array for all slices.
     *
     * @return Vector of CMT normalization factors
     */
    Eigen::VectorXd get_norm_cmt_array() const;

    /**
     * @brief Get scalar coupling normalization array for all slices.
     *
     * @return Vector of scalar coupling normalization factors
     */
    Eigen::VectorXd get_norm_scalar_coupling_array() const;

    /**
     * @brief Get L2 normalization array for all slices.
     *
     * @return Vector of L2 normalization factors
     */
    Eigen::VectorXd get_norm_l2_array() const;

    /**
     * @brief Get maximum field normalization array for all slices.
     *
     * @return Vector of maximum normalization factors
     */
    Eigen::VectorXd get_norm_max_array() const;

    /**
     * @brief Calculate overlap integral with another SuperMode at a specific slice.
     *
     * @param other_supermode The other SuperMode to calculate overlap with
     * @param slice The slice index for the calculation
     * @return Overlap integral value
     */
    double get_overlap_integral(const SuperMode &other_supermode, size_t slice) const;

    /**
     * @brief Calculate overlap integral between two slices of different SuperModes.
     *
     * @param other_supermode The other SuperMode
     * @param slice_0 Slice index from this SuperMode
     * @param slice_1 Slice index from the other SuperMode
     * @return Overlap integral value
     */
    double get_overlap_integral(const SuperMode &other_supermode, size_t slice_0, size_t slice_1) const;

    /**
     * @brief Calculate overlap integrals matrix with another SuperMode.
     *
     * @param other_supermode The other SuperMode
     * @return Matrix of overlap integrals for all slice combinations
     */
    Eigen::MatrixXd get_overlap_integrals_with_mode(const SuperMode& other_supermode) const;

    /**
     * @brief Normalize the electromagnetic fields of this SuperMode.
     *
     * Applies field normalization to ensure proper scaling and physical consistency.
     */
    void normalize_fields();

    /**
     * @brief Arrange and organize the field data.
     *
     * Performs field arrangement operations for consistent data organization.
     */
    void arrange_fields();

    /**
     * @brief Calculate L2 norm of field slice.
     *
     * @param slice The slice index
     * @return L2 norm of the field at the specified slice
     */
    double get_field_slice_norm_l2(const size_t slice) const;

    /**
     * @brief Perform trapezoidal integration over a 2D mesh.
     *
     * @param mesh The 2D mesh data for integration
     * @param dx Grid spacing in x-direction
     * @param dy Grid spacing in y-direction
     * @return Integrated value over the mesh
     */
    double get_trapez_integral(const Eigen::VectorXd &mesh, const double &dx, const double &dy) const;

    /**
     * @brief Calculate normalized coupling coefficients with another SuperMode.
     *
     * @param other_supermode The other SuperMode to calculate coupling with
     * @return Complex vector of normalized coupling coefficients
     */
    Eigen::Matrix<complex128, Eigen::Dynamic, 1> get_normalized_coupling_with_mode(const SuperMode &other_supermode) const;

    /**
     * @brief Calculate beating length between this mode and another.
     *
     * @param other_supermode The other SuperMode
     * @return Vector of beating lengths along the propagation direction
     */
    Eigen::VectorXd get_beating_length_with_mode(const SuperMode &other_supermode) const;

    /**
     * @brief Calculate adiabatic parameter with another SuperMode.
     *
     * @param other_supermode The other SuperMode
     * @return Vector of adiabatic parameters along the propagation direction
     */
    Eigen::VectorXd get_adiabatic_with_mode(const SuperMode &other_supermode) const;

    /**
     * @brief Calculate gradient field overlap with another SuperMode.
     *
     * @param other_supermode The other SuperMode
     * @return Vector of gradient field overlap values
     */
    Eigen::VectorXd get_gradient_field_overlap(const SuperMode &other_supermode) const;

    /**
     * @brief Perform trapezoidal integration over a 2D matrix.
     *
     * @param mesh The 2D matrix data for integration
     * @param dx Grid spacing in x-direction
     * @param dy Grid spacing in y-direction
     * @return Integrated value over the matrix
     */
    double get_trapz_integral(const Eigen::MatrixXd& mesh, double dx, double dy) const;

    /**
     * @brief Check if this SuperMode has the same symmetry as another.
     *
     * @param other_supermode The other SuperMode to compare with
     * @return True if both modes have the same boundary conditions
     */
    bool is_same_symmetry(const SuperMode &other_supermode) const;

    /**
     * @brief Check if this SuperMode is computation compatible with another.
     *
     * Two modes are computation compatible if they have the same boundary
     * conditions but are different modes (different mode numbers).
     *
     * @param other_supermode The other SuperMode to check compatibility with
     * @return True if modes are compatible for computations
     */
    bool is_computation_compatible(const SuperMode &other_supermode) const {
        // Return True if both mode have same boundary condisitons but are different
        return (this->mode_number != other_supermode.mode_number) && (this->is_same_symmetry(other_supermode));
    }

    // Python interface methods using pybind11

    /**
     * @brief Get overlap integrals with another mode (Python interface).
     *
     * @param supermode The other SuperMode
     * @return Python numpy array of overlap integrals
     */
    pybind11::array_t<double> get_overlap_integrals_with_mode_py(const SuperMode& supermode) const {
        return numy_interface::eigen_to_ndarray<double>(this->get_overlap_integrals_with_mode(supermode), {model_parameters.n_slice});
    }

    /**
     * @brief Get normalized coupling coefficients (Python interface).
     *
     * @param supermode The other SuperMode
     * @return Python numpy array of complex coupling coefficients
     */
    pybind11::array_t<std::complex<double>> get_normalized_coupling_with_mode_py(const SuperMode& supermode) const {
        return numy_interface::eigen_to_ndarray<complex128>(this->get_normalized_coupling_with_mode(supermode), {model_parameters.n_slice});
    }

    /**
     * @brief Get adiabatic parameters with another mode (Python interface).
     *
     * @param supermode The other SuperMode
     * @return Python numpy array of adiabatic parameters
     */
    pybind11::array_t<double> get_adiabatic_with_mode_py(const SuperMode& supermode) const {
        return numy_interface::eigen_to_ndarray<double>(this->get_adiabatic_with_mode(supermode), {model_parameters.n_slice});
    }

    /**
     * @brief Get beating lengths with another mode (Python interface).
     *
     * @param supermode The other SuperMode
     * @return Python numpy array of beating lengths
     */
    pybind11::array_t<double> get_beating_length_with_mode_py(const SuperMode& supermode) const {
        return numy_interface::eigen_to_ndarray<double>(this->get_beating_length_with_mode(supermode), {model_parameters.n_slice});
    }

    /**
     * @brief Get electromagnetic fields (Python interface).
     *
     * @return Python numpy array of field data with shape [n_slice, ny, nx]
     */
    pybind11::array_t<double> get_fields_py() const {
        return numy_interface::eigen_to_ndarray<double>(this->fields, {model_parameters.n_slice, model_parameters.ny, model_parameters.nx});
    }

    /**
     * @brief Get refractive index profile (Python interface).
     *
     * @return Python numpy array of refractive indices with shape [n_slice]
     */
    pybind11::array_t<double> get_index_py() const {
        return numy_interface::eigen_to_ndarray<double>(this->index, {model_parameters.n_slice});
    }

    /**
     * @brief Get eigenvalues (Python interface).
     *
     * @return Python numpy array of eigenvalues with shape [n_slice]
     */
    pybind11::array_t<double> get_eigen_value_py() const {
        return numy_interface::eigen_to_ndarray<double>(this->eigen_value, {model_parameters.n_slice});
    }

    /**
     * @brief Get propagation constants (Python interface).
     *
     * @return Python numpy array of beta values with shape [n_slice]
     */
    pybind11::array_t<double> get_betas_py() const {
        return numy_interface::eigen_to_ndarray<double>(this->betas, {model_parameters.n_slice});
    }

    /**
     * @brief Get iteration list from model parameters (Python interface).
     *
     * @return Python numpy array of iteration values
     */
    pybind11::array_t<double> get_itr_list() const {
        return numy_interface::eigen_to_ndarray<double>(this->model_parameters.itr_list, {model_parameters.n_slice});
    }

    /**
     * @brief Get mesh gradient (Python interface).
     *
     * @return Python numpy array of mesh gradient with shape [nx, ny]
     */
    pybind11::array_t<double> get_mesh_gradient() const {
        return numy_interface::eigen_to_ndarray<double>(this->model_parameters.mesh_gradient, {this->model_parameters.nx, this->model_parameters.ny});
    }

    /**
     * @brief Get field data for a specific slice (Python interface).
     *
     * @param index The slice index to retrieve
     * @return Python numpy array of field data with shape [nx, ny]
     */
    pybind11::array_t<double> get_field_py(const size_t index) const {
        Eigen::MatrixXd output = this->fields.col(index);
        return numy_interface::eigen_to_ndarray<double>(output, {this->model_parameters.nx, this->model_parameters.ny});
    }

    // Serialization methods for Python pickling

    /**
     * @brief Create a pickle tuple for Python serialization.
     *
     * @param supermode The SuperMode to serialize
     * @return Python tuple containing all serializable data
     */
    static pybind11::tuple get_pickle(SuperMode &supermode);

    /**
     * @brief Reconstruct SuperMode from pickle tuple.
     *
     * @param tuple Python tuple containing serialized SuperMode data
     * @return Reconstructed SuperMode object
     */
    static SuperMode build_from_tuple(pybind11::tuple tuple);

};


