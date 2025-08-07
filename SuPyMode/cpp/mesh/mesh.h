#pragma once

#include <pybind11/pybind11.h>
#include <utility>
#include <Eigen/Dense>
#include "../utils/numpy_interface.h"
#include <unsupported/Eigen/MatrixFunctions>
#include <cmath>


std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
compute_gradient_2p(const Eigen::MatrixXd& image, double dx, double dy);

// Rows is y-axis -- Cols is x--axis
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
compute_gradient_5p(const Eigen::MatrixXd& image, const double dx, const double dy);

std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
compute_gradient_7p(const Eigen::MatrixXd& image, double dx, double dy);

Eigen::MatrixXd
get_rho_gradient_time_rho(const Eigen::MatrixXd &mesh, const Eigen::VectorXd &y_vector, const Eigen::VectorXd &x_vector);

Eigen::MatrixXd
get_rho_gradient(const Eigen::MatrixXd &mesh, const Eigen::VectorXd &y_vector, const Eigen::VectorXd &x_vector);

pybind11::array_t<double>
get_rho_gradient_py(const pybind11::array_t<double> &mesh_py, const pybind11::array_t<double> &x_vector_py, const pybind11::array_t<double> &y_vector_py);

pybind11::array_t<double>
get_rho_gradient_time_rho_py(const pybind11::array_t<double> &mesh_py, const pybind11::array_t<double> &x_vector_py, const pybind11::array_t<double> &y_vector_py);


class Geometry {
public:
    Geometry() = default;

    Geometry(const pybind11::array_t<double>& mesh, const pybind11::array_t<double>& x, const pybind11::array_t<double>& y)
    {
        this->mesh = NumpyInterface::convert_py_to_eigen_matrix(mesh);
        this->x = NumpyInterface::convert_py_to_eigen_vector(x);
        this->y = NumpyInterface::convert_py_to_eigen_vector(y);
        dx = this->x[1] - this->x[0];
        dy = this->y[1] - this->y[0];
        dx = std::abs(dx);
        dy = std::abs(dy);
    }

    /**
     * @brief Computes the gradient of the mesh using a 2-point finite difference method.
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_gradient_2p();

    /**
     * @brief Computes the gradient of the mesh using a 5-point finite difference method.
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_gradient_5p();

    /**
     * @brief Computes the gradient of the mesh using a 7-point finite difference method.
     */
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_gradient_7p();

    /**
     * @brief Computes the rho gradient of the mesh.
     */
    Eigen::MatrixXd get_rho_gradient_time_rho();

    /**
     * @brief Computes the rho gradient of the mesh.
     */
    Eigen::MatrixXd get_rho_gradient();


private:
    Eigen::MatrixXd mesh;
    Eigen::VectorXd x;
    Eigen::VectorXd y;
    double dx;
    double dy;
};



