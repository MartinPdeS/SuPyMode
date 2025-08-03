#pragma once

#include <pybind11/pybind11.h>
#include <utility>
#include <Eigen/Dense>
#include "../utils/numpy_interface.h"
#include <unsupported/Eigen/MatrixFunctions>


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






