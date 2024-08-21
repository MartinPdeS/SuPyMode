#pragma once


#include "model_parameters.cpp"
#include "numpy_interface.cpp"
#include <Eigen/Core>

typedef std::complex<double> complex128;

class SuperMode
{
public:
    size_t mode_number;
    Eigen::MatrixXd fields;
    Eigen::VectorXd index, betas, eigen_value;
    ModelParameters model_parameters;

    std::string left_boundary;
    std::string right_boundary;
    std::string top_boundary;
    std::string bottom_boundary;

    SuperMode() = default;
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
        const std::string &bottom_boundary) : mode_number(mode_number), model_parameters(model_parameters)
    {
        fields = convert_py_to_eigen(fields_py, this->model_parameters.nx * this->model_parameters.ny, this->model_parameters.n_slice);
        index = convert_py_to_eigen(index_py, this->model_parameters.n_slice, 1);
        betas = convert_py_to_eigen(betas_py, this->model_parameters.n_slice, 1);
        eigen_value = convert_py_to_eigen(eigen_value_py, this->model_parameters.n_slice, 1);
    }

    SuperMode(
        const size_t mode_number,
        const ModelParameters &model_parameters,
        const std::string& left_boundary,
        const std::string& right_boundary,
        const std::string& top_boundary,
        const std::string& bottom_boundary) : mode_number(mode_number)
    {
        this->model_parameters = model_parameters;
        this->fields = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(this->model_parameters.nx * this->model_parameters.ny, model_parameters.n_slice);
        this->eigen_value = Eigen::Matrix<double, Eigen::Dynamic, 1>(model_parameters.n_slice);
        this->betas = Eigen::Matrix<double, Eigen::Dynamic, 1>(model_parameters.n_slice);
        this->index = Eigen::Matrix<double, Eigen::Dynamic, 1>(model_parameters.n_slice);
    }

    double get_norm(const size_t &slice, const std::string &normalization_type) const;

    Eigen::VectorXd get_norm_array(const std::string &normalization_type) const;

    double get_norm_cmt(const size_t &slice) const;
    double get_norm_max(const size_t &slice) const;
    double get_norm_l2(const size_t &slice) const;

    Eigen::VectorXd get_norm_cmt_array() const;
    Eigen::VectorXd get_norm_scalar_coupling_array() const;
    Eigen::VectorXd get_norm_l2_array() const;
    Eigen::VectorXd get_norm_max_array() const;

    double get_overlap_integral(const SuperMode &other_supermode, size_t slice) const;
    double get_overlap_integral(const SuperMode &other_supermode, size_t slice_0, size_t slice_1) const;

    Eigen::MatrixXd get_overlap_integrals_with_mode(const SuperMode& other_supermode) const;

    void normalize_fields();
    void arrange_fields();

    double get_field_slice_norm_l2(const size_t slice) const;

    double get_trapez_integral(const Eigen::VectorXd &mesh, const double &dx, const double &dy) const;

    Eigen::Matrix<complex128, Eigen::Dynamic, 1> get_normalized_coupling_with_mode(const SuperMode &other_supermode) const;

    Eigen::VectorXd get_beating_length_with_mode(const SuperMode &other_supermode) const;

    Eigen::VectorXd get_adiabatic_with_mode(const SuperMode &other_supermode) const;

    Eigen::VectorXd get_gradient_field_overlap(const SuperMode &other_supermode) const;

    double get_trapz_integral(const Eigen::MatrixXd& mesh, double dx, double dy) const;

    bool is_same_symmetry(const SuperMode &other_supermode) const {
        // Return True if all boundaries condition are the same else False
        return
            (this->left_boundary == other_supermode.left_boundary) &&
            (this->right_boundary == other_supermode.right_boundary) &&
            (this->top_boundary == other_supermode.top_boundary) &&
            (this->bottom_boundary == other_supermode.bottom_boundary);
    }

    bool is_computation_compatible(const SuperMode &other_supermode) const {
        // Return True if both mode have same boundary condisitons but are different
        return (this->mode_number != other_supermode.mode_number) && (this->is_same_symmetry(other_supermode));
    }

    pybind11::array_t<double> get_overlap_integrals_with_mode_py(const SuperMode& supermode) const {
        return eigen_to_ndarray<double>(this->get_overlap_integrals_with_mode(supermode), {model_parameters.n_slice});
    }

    pybind11::array_t<std::complex<double>> get_normalized_coupling_with_mode_py(const SuperMode& supermode) const {
        return eigen_to_ndarray<complex128>(this->get_normalized_coupling_with_mode(supermode), {model_parameters.n_slice});
    }

    pybind11::array_t<double> get_adiabatic_with_mode_py(const SuperMode& supermode) const {
        return eigen_to_ndarray<double>(this->get_adiabatic_with_mode(supermode), {model_parameters.n_slice});
    }

    pybind11::array_t<double> get_beating_length_with_mode_py(const SuperMode& supermode) const {
        return eigen_to_ndarray<double>(this->get_beating_length_with_mode(supermode), {model_parameters.n_slice});
    }

    pybind11::array_t<double> get_fields_py() const {
        return eigen_to_ndarray<double>(this->fields, {model_parameters.n_slice, model_parameters.ny, model_parameters.nx});
    }

    pybind11::array_t<double> get_index_py() const {
        return eigen_to_ndarray<double>(this->index, {model_parameters.n_slice});
    }

    pybind11::array_t<double> get_eigen_value_py() const {
        return eigen_to_ndarray<double>(this->eigen_value, {model_parameters.n_slice});
    }

    pybind11::array_t<double> get_betas_py() const {
        return eigen_to_ndarray<double>(this->betas, {model_parameters.n_slice});
    }

    pybind11::array_t<double> get_itr_list() const {
        return eigen_to_ndarray<double>(this->model_parameters.itr_list, {model_parameters.n_slice});
    }

    pybind11::array_t<double> get_mesh_gradient() const {
        return eigen_to_ndarray<double>(this->model_parameters.mesh_gradient, {this->model_parameters.nx, this->model_parameters.ny});
    }

    pybind11::array_t<double> get_field_py(const size_t index) const {
        Eigen::MatrixXd output = this->fields.col(index);
        return eigen_to_ndarray<double>(output, {this->model_parameters.nx, this->model_parameters.ny});
    }

    static pybind11::tuple get_pickle(SuperMode &supermode) {
        return pybind11::make_tuple(
            supermode.mode_number,
            supermode.get_fields_py(),
            supermode.get_index_py(),
            supermode.get_betas_py(),
            supermode.get_eigen_value_py(),
            supermode.model_parameters,
            supermode.left_boundary,
            supermode.right_boundary,
            supermode.top_boundary,
            supermode.bottom_boundary
        );
    }

    static SuperMode build_from_tuple(pybind11::tuple tuple) {
        return SuperMode{
            tuple[0].cast<size_t>(),                             // mode_number
            tuple[1].cast<pybind11::array_t<double>>(),          // fields
            tuple[2].cast<pybind11::array_t<double>>(),          // index
            tuple[3].cast<pybind11::array_t<double>>(),          // betas
            tuple[4].cast<pybind11::array_t<double>>(),          // eigen_values
            tuple[5].cast<ModelParameters>(),                    // Model parameters
            tuple[6].cast<std::string>(),                        // left_boundary
            tuple[7].cast<std::string>(),                        // right_boundary
            tuple[8].cast<std::string>(),                        // top_boundary
            tuple[9].cast<std::string>()                         // bottom_boundary
        }; // load
    }
};


