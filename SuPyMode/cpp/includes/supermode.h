#ifndef SUPERMODE_H
#define SUPERMODE_H

#include "definitions.cpp"
#include "numpy_interface.cpp"

class SuperMode
{
public:
    size_t mode_number;
    Eigen::MatrixXd fields;
    Eigen::VectorXd index, betas, eigen_value;
    ModelParameters model_parameters;

    SuperMode();
    SuperMode(
        size_t mode_number,
        pybind11::array_t<double> fields_py,
        pybind11::array_t<double> index_py,
        pybind11::array_t<double> betas_py,
        pybind11::array_t<double> eigen_value_py,
        ModelParameters model_parameters);


    SuperMode(const size_t mode_number, const ModelParameters &model_parameters);


    pybind11::tuple get_state();


    double get_norm(const size_t &slice, const std::string &normalization_type) const;


    Eigen::VectorXd get_norm_array(const std::string &normalization_type) const;


    double get_norm_cmt(const size_t &slice) const;


    double get_norm_max(const size_t &slice) const;


    double get_norm_l2(const size_t &slice) const;


    Eigen::VectorXd get_norm_cmt_array() const;


    Eigen::VectorXd get_norm_scalar_coupling_array() const;


    Eigen::VectorXd get_norm_l2_array() const;


    Eigen::VectorXd get_norm_max_array() const;


    double get_overlap_integral(const SuperMode& other_supermode, const size_t &slice) const;


    double get_overlap_integral(const SuperMode& other_supermode, const size_t &slice_0, const size_t &slice_1) const;


    Eigen::MatrixXd get_overlap_integrals_with_mode(const SuperMode& other_supermode) const;


    void normalize_fields();


    void arrange_fields();


    void normalize_field_slice_l2(const size_t slice);


    double get_field_slice_norm_custom(const size_t slice) const;


    double get_trapez_integral(const Eigen::VectorXd &mesh, const double &dx, const double &dy) const;


    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> get_normalized_coupling_with_mode(const SuperMode& other_supermode) const;


    Eigen::VectorXd get_beating_length_with_mode(const SuperMode& other_supermode) const;


    Eigen::VectorXd get_adiabatic_with_mode(const SuperMode& other_supermode) const;


    pybind11::array_t<double> get_overlap_integrals_with_mode_py(const SuperMode& supermode) const;


    pybind11::array_t<std::complex<double>> get_normalized_coupling_with_mode_py(const SuperMode& supermode) const;


    pybind11::array_t<double> get_adiabatic_with_mode_py(const SuperMode& supermode) const;


    pybind11::array_t<double> get_beating_length_with_mode_py(const SuperMode& supermode) const;


    pybind11::array_t<double> get_fields_py();


    pybind11::array_t<double> get_index_py();


    pybind11::array_t<double> get_eigen_value_py();


    pybind11::array_t<double> get_betas_py();


    pybind11::array_t<double> get_itr_list();


    pybind11::array_t<double> get_mesh_gradient();


    SuperMode get_supermode_from_tuple(pybind11::tuple tuple);
};

#endif

