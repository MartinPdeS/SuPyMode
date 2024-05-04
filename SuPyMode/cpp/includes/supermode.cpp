#pragma once

#include "definitions.cpp"
#include "utils.cpp"
#include "supermode.h"


typedef std::complex<double> complex128;

pybind11::tuple SuperMode::get_state()
{
    return pybind11::make_tuple(
        this->mode_number,
        this->get_fields_py(),
        this->get_index_py(),
        this->get_betas_py(),
        this->get_eigen_value_py(),
        this->model_parameters
        );
}

SuperMode SuperMode::get_supermode_from_tuple(pybind11::tuple tuple) const
{
    return SuperMode{
        tuple[0].cast<size_t>(),                             // mode_number
        tuple[1].cast<pybind11::array_t<double>>(),          // fields
        tuple[2].cast<pybind11::array_t<double>>(),          // index
        tuple[3].cast<pybind11::array_t<double>>(),          // betas
        tuple[4].cast<pybind11::array_t<double>>(),          // eigen_values
        tuple[5].cast<ModelParameters>()                     // model parameters
    }; // load
}

double SuperMode::get_norm(const size_t &slice, const std::string &normalization_type) const
{
    if (normalization_type == "max")
        return this->get_norm_max(slice);
    if (normalization_type == "l2")
        return this->get_norm_l2(slice);
    if (normalization_type == "cmt")
        return this->get_norm_cmt(slice);
}

Eigen::VectorXd SuperMode::get_beating_length_with_mode(const SuperMode &other_supermode) const {

    Eigen::VectorXd beta_0 = this->betas, beta_1 = other_supermode.betas;

    return (beta_0 - beta_1).cwiseAbs().cwiseInverse() * (2 * PI);
}


Eigen::VectorXd SuperMode::get_adiabatic_with_mode(const SuperMode &other_supermode) const {
    Eigen::VectorXd delta_beta = (this->betas - other_supermode.betas).cwiseAbs();

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>
        coupling = this->get_normalized_coupling_with_mode(other_supermode).cwiseAbs();

    return delta_beta.cwiseProduct(coupling.cwiseInverse()).cwiseAbs();
}

Eigen::Matrix<complex128, Eigen::Dynamic, 1> SuperMode::get_normalized_coupling_with_mode(const SuperMode &other_supermode) const {
    Eigen::VectorXd
        integral = this->get_gradient_field_overlap(other_supermode),
        beta_0 = this->betas,
        beta_1 = other_supermode.betas,
        delta_beta = (beta_0 - beta_1).cwiseInverse(),
        beta_prod = (beta_0.cwiseProduct(beta_1)).cwiseSqrt().cwiseInverse(),
        denominator = delta_beta.cwiseProduct(beta_prod);

    double k2 = pow(model_parameters.wavenumber, 2);

    std::complex<double> scalar = -J * k2 / 2.0;

    return scalar * denominator.cwiseProduct(integral);
}

Eigen::VectorXd SuperMode::get_gradient_field_overlap(const SuperMode &other_supermode) const
{
    Eigen::VectorXd output(model_parameters.n_slice);

    Eigen::MatrixXd
        mesh_2D = model_parameters.mesh_gradient,
        field_0,
        field_1,
        overlap;

    double gradient_overlap;
    for (size_t slice = 0; slice < model_parameters.n_slice; ++slice)
    {
        field_0 = this->fields.col(slice);
        field_0.resize(model_parameters.nx, model_parameters.ny);

        field_1 = other_supermode.fields.col(slice);
        field_1.resize(model_parameters.nx, model_parameters.ny);

        overlap = mesh_2D.cwiseProduct(field_0).cwiseProduct(field_1);

        gradient_overlap = get_trapz_integral(
            overlap,
            model_parameters.dx_scaled[slice],
            model_parameters.dy_scaled[slice]
        );

        output[slice] = gradient_overlap;
    }

    return output;
}

double SuperMode::get_norm_cmt(const size_t &slice) const // Equation 7.35 from Bures
{
    double
        itr = model_parameters.itr_list[slice],
        factor = 0.5 * this->model_parameters.dx * itr * this->model_parameters.dy * itr;

    return this->fields.col(slice).cwiseAbs2().sum() * factor;
}

double SuperMode::get_norm_max(const size_t &slice) const
{
    return this->fields.col(slice).cwiseAbs().maxCoeff();
}

double SuperMode::get_norm_l2(const size_t &slice) const
{
    return sqrt(this->fields.col(slice).cwiseAbs2().sum());
}

Eigen::VectorXd SuperMode::get_norm_cmt_array() const
{
    Eigen::VectorXd
        term_0 = 0.5 * this->fields.cwiseAbs2().colwise().sum(),
        term_1 = model_parameters.itr_list.cwiseProduct(model_parameters.itr_list) * this->model_parameters.dx * this->model_parameters.dy,
        norm_array = term_0.cwiseProduct(term_1);

    return norm_array;
}

Eigen::VectorXd SuperMode::get_norm_scalar_coupling_array() const
{
    return this->fields.cwiseAbs2().colwise().sum() * (2 * PI);
}

Eigen::VectorXd SuperMode::get_norm_max_array() const
{
    return this->fields.colwise().maxCoeff();
}

double SuperMode::get_overlap_integral(const SuperMode& other_supermode, size_t slice) const
{
    return this->get_overlap_integral(other_supermode, slice, slice);
}

double SuperMode::get_overlap_integral(const SuperMode& other_supermode, size_t slice_0, size_t slice_1) const
{
    return this->fields.col(slice_0).transpose() * other_supermode.fields.col(slice_0);
}

Eigen::MatrixXd SuperMode::get_overlap_integrals_with_mode(const SuperMode& other_supermode) const
{
    Eigen::VectorXd
        overlap_integrals = this->fields.cwiseProduct(other_supermode.fields).colwise().sum().cwiseAbs();

    return overlap_integrals;
}

void SuperMode::normalize_fields()
{
    for(size_t slice = 0; slice < this->fields.cols(); ++ slice)
        this->normalize_field_slice_l2(slice);
}

void SuperMode::arrange_fields()
{
    for(size_t slice = 0; slice < this->fields.cols() - 1; ++ slice)
    {
        Eigen::VectorXd
            field_0 = this->fields.col(slice),
            field_1 = this->fields.col(slice + 1);


        double overlap = field_0.cwiseProduct(field_1).sum();
        if (overlap < 0)
            this->fields.col(slice + 1) *= -1;
    }
}

void SuperMode::normalize_field_slice_l2(const size_t slice)
{
    double norm = this->get_field_slice_norm_custom(slice);

    this->fields.col(slice) /= sqrt(norm);
}

double SuperMode::get_field_slice_norm_custom(const size_t slice) const
{
    Eigen::VectorXd field = this->fields.col(slice);

    double
        itr = this->model_parameters.itr_list[slice],
        dx = this->model_parameters.dx * itr,
        dy = this->model_parameters.dy * itr,
        dA = dx * dy,
        norm = field.cwiseAbs2().sum() * dA;

    return norm;
}





