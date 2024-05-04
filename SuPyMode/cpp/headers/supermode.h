#pragma once


#include "model_parameters.cpp"
#include "numpy_interface.cpp"
#include <eigen/Eigen>
#include "equations.cpp"

typedef std::complex<double> complex128;

class SuperMode
{
public:
    size_t mode_number;
    Eigen::MatrixXd fields;
    Eigen::VectorXd index, betas, eigen_value;
    ModelParameters model_parameters;

    SuperMode() = default;
    SuperMode(size_t mode_number, pybind11::array_t<double> fields_py, pybind11::array_t<double> index_py,
    pybind11::array_t<double> betas_py, pybind11::array_t<double> eigen_value_py, ModelParameters model_parameters) :
    mode_number(mode_number), model_parameters(model_parameters)
    {
        fields = convert_py_to_eigen(fields_py, this->model_parameters.field_size, this->model_parameters.n_slice);
        index = convert_py_to_eigen(index_py, this->model_parameters.n_slice, 1);
        betas = convert_py_to_eigen(betas_py, this->model_parameters.n_slice, 1);
        eigen_value = convert_py_to_eigen(eigen_value_py, this->model_parameters.n_slice, 1);
    }

    SuperMode(const size_t mode_number, const ModelParameters &model_parameters) :
    mode_number(mode_number)
    {
        this->model_parameters = model_parameters;
        this->fields = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(this->model_parameters.field_size, model_parameters.n_slice);
        this->eigen_value = Eigen::Matrix<double, Eigen::Dynamic, 1>(model_parameters.n_slice);
        this->betas = Eigen::Matrix<double, Eigen::Dynamic, 1>(model_parameters.n_slice);
        this->index = Eigen::Matrix<double, Eigen::Dynamic, 1>(model_parameters.n_slice);
    }

    Eigen::VectorXd convert_py_to_eigen(pybind11::array_t<double> array_py, size_t rows, size_t cols) const {
        auto info = array_py.request();
        double* ptr = static_cast<double*>(info.ptr);
        return Eigen::Map<Eigen::VectorXd>(ptr, rows, cols);
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

    double get_overlap_integral(const SuperMode& other_supermode, size_t slice) const;
    double get_overlap_integral(const SuperMode& other_supermode, size_t slice_0, size_t slice_1) const;

    Eigen::MatrixXd get_overlap_integrals_with_mode(const SuperMode& other_supermode) const;

    void normalize_fields();
    void arrange_fields();
    void normalize_field_slice_l2(const size_t slice);

    double get_field_slice_norm_custom(const size_t slice) const;

    double get_trapez_integral(const Eigen::VectorXd &mesh, const double &dx, const double &dy) const;

    Eigen::Matrix<complex128, Eigen::Dynamic, 1> get_normalized_coupling_with_mode(const SuperMode &other_supermode) const;

    Eigen::VectorXd get_beating_length_with_mode(const SuperMode &other_supermode) const;

    Eigen::VectorXd get_adiabatic_with_mode(const SuperMode &other_supermode) const;

    Eigen::VectorXd get_gradient_field_overlap(const SuperMode &other_supermode) const;




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
        return eigen_to_ndarray<double>(this->fields, {model_parameters.n_slice, model_parameters.nx, model_parameters.ny});
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

    SuperMode get_supermode_from_tuple(pybind11::tuple tuple) const;
    pybind11::tuple get_state();
};


