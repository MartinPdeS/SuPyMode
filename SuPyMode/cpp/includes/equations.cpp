#pragma once

#include <eigen/Eigen>

double get_trapz_integral(const Eigen::MatrixXd &mesh, double dx, double dy)
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

double get_basic_integral(const Eigen::MatrixXd &mesh, double dx, double dy)
{
    double
        dA = dx * dy,
        integral = mesh.sum() * dA;

    return integral;
}