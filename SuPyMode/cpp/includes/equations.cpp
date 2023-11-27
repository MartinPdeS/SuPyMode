#pragma once

#include "definitions.cpp"
#include "supermode.h"

double
get_trapz_integral(const Eigen::MatrixXd &mesh, const double &dx, const double &dy)
{
    Eigen::VectorXd
        vector_integral(mesh.rows());

    for (size_t row = 0; row < mesh.rows(); ++row)
    {
        double integral_cols = 0.0;
        for (size_t col = 0; col < mesh.cols() - 1; ++col)
            integral_cols += 0.5 * (mesh(row, col) + mesh(row, col + 1)) * dx;

        vector_integral[row] = integral_cols;
    }

    double integral = 0.0;
    for (size_t row = 0; row< mesh.rows() - 1; ++row)
        integral += 0.5 * (vector_integral[row] + vector_integral[row + 1]) * dy;

    return integral;
}

double
get_basic_integral(const Eigen::MatrixXd &mesh, const double &dx, const double &dy)
{
    double
        dA = dx * dy,
        integral = mesh.sum() * dA;

    return integral;
}



Eigen::VectorXd
get_gradient_field_overlap(const SuperMode& mode_0, const SuperMode& mode_1, const ModelParameters &model_parameters)
{
    Eigen::VectorXd
        output(model_parameters.n_slice);

    Eigen::MatrixXd
        mesh_2D = model_parameters.mesh_gradient,
        field_0,
        field_1,
        overlap;

    double gradient_overlap;
    for (size_t slice = 0; slice < model_parameters.n_slice; ++slice)
    {
        field_0 = mode_0.fields.col(slice);
        field_0.resize(model_parameters.nx, model_parameters.ny);

        field_1 = mode_1.fields.col(slice);
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

Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>
get_mode_normalized_coupling(const SuperMode& mode_0, const SuperMode& mode_1, const ModelParameters &model_parameters)
{
    Eigen::VectorXd
        integral = get_gradient_field_overlap(mode_0, mode_1, model_parameters),
        beta_0 = mode_0.betas,
        beta_1 = mode_1.betas,
        delta_beta = (beta_0 - beta_1).cwiseInverse(),
        beta_prod = (beta_0.cwiseProduct(beta_1)).cwiseSqrt().cwiseInverse(),
        denominator = delta_beta.cwiseProduct(beta_prod);

    double k2 = pow(model_parameters.wavenumber, 2);

    std::complex<double> scalar = -J * k2 / 2.0;

    return scalar * denominator.cwiseProduct(integral);
}

Eigen::VectorXd
get_mode_adiabatic_with_mode(const SuperMode& mode_0, const SuperMode& mode_1, const ModelParameters &model_parameters)
{
    Eigen::VectorXd delta_beta = (mode_1.betas - mode_0.betas).cwiseAbs();

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>
        coupling = get_mode_normalized_coupling(mode_0, mode_1, model_parameters).cwiseAbs();

    return delta_beta.cwiseProduct(coupling.cwiseInverse()).cwiseAbs();
}

Eigen::VectorXd
get_mode_beating_length_with_mode(const SuperMode& mode_0, const SuperMode& mode_1, const ModelParameters &model_parameters)
{
    Eigen::VectorXd
        beta_0 = mode_0.betas,
        beta_1 = mode_1.betas;

    return (beta_0 - beta_1).cwiseAbs().cwiseInverse() * (2 * PI);
}