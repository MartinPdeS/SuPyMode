#pragma once

#include <Eigen/Sparse>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "numpy_interface.cpp"
#include "mesh.cpp"

#define PI 3.1415926535897932384626f


class ModelParameters {
public:
    double wavelength;
    double wavenumber;
    double dx;
    double dy;
    double ditr;

    size_t nx;
    size_t ny;
    size_t n_slice;

    pybind11::array_t<double> mesh_py;
    pybind11::array_t<double> x_vector_py;
    pybind11::array_t<double> y_vector_py;
    pybind11::array_t<double> itr_list_py;
    pybind11::array_t<double> mesh_gradient_py;

    Eigen::MatrixXd mesh;
    Eigen::VectorXd x_vector;
    Eigen::VectorXd y_vector;
    Eigen::MatrixXd mesh_gradient;
    Eigen::VectorXd itr_list;
    Eigen::VectorXd dx_scaled;
    Eigen::VectorXd dy_scaled;
    Eigen::VectorXd wavenumber_scaled;

    int debug_mode = 0;

    ModelParameters() = default;

    ModelParameters(
        const pybind11::array_t<double> &mesh_py,
        const pybind11::array_t<double> &x_vector_py,
        const pybind11::array_t<double> &y_vector_py,
        const pybind11::array_t<double> &itr_list_py,
        const double wavelength,
        const double dx,
        const double dy,
        const int debug_mode = 0)
        : wavelength(wavelength), dx(dx), dy(dy),
          debug_mode(debug_mode),
          mesh_py(mesh_py), x_vector_py(x_vector_py),
          y_vector_py(y_vector_py), itr_list_py(itr_list_py)
    {
        init();
    }

    void compute_scaled_parameters() {
        dx_scaled = itr_list * dx;
        dy_scaled = itr_list * dy;
        wavenumber_scaled = itr_list * wavenumber;
    }

    static pybind11::tuple get_pickle(const ModelParameters &model_parameters) {
        return pybind11::make_tuple(
            model_parameters.mesh_py,
            model_parameters.x_vector_py,
            model_parameters.y_vector_py,
            model_parameters.itr_list_py,
            model_parameters.wavelength,
            model_parameters.dx,
            model_parameters.dy,
            model_parameters.debug_mode
        );
    }

    static ModelParameters build_from_tuple(const pybind11::tuple &tuple) {
        return ModelParameters{
            tuple[0].cast<pybind11::array_t<double>>(),
            tuple[1].cast<pybind11::array_t<double>>(),
            tuple[2].cast<pybind11::array_t<double>>(),
            tuple[3].cast<pybind11::array_t<double>>(),
            tuple[4].cast<double>(),
            tuple[5].cast<double>(),
            tuple[6].cast<double>(),
            tuple[7].cast<int>()
        };
    }

    friend std::ostream &operator<<(std::ostream &os, const ModelParameters &params) {
        os << "dx: " << params.dx
           << " dy: " << params.dy
           << " nx: " << params.nx
           << " ny: " << params.ny
           << " n_slice: " << params.n_slice;
        return os;
    }

private:
    void init() {
        nx = mesh_py.request().shape[0];
        ny = mesh_py.request().shape[1];
        n_slice = itr_list_py.request().size;
        wavenumber = 2 * PI / wavelength;

        mesh = convert_py_to_eigen<double>(mesh_py, nx, ny);
        x_vector = convert_py_to_eigen<double>(x_vector_py, nx);
        y_vector = convert_py_to_eigen<double>(y_vector_py, ny);
        itr_list = convert_py_to_eigen<double>(itr_list_py, n_slice);
        mesh_gradient = get_rho_gradient_time_rho(mesh.cwiseProduct(mesh), x_vector, y_vector);
        mesh_gradient_py = eigen_to_ndarray<double>(mesh_gradient, {nx, ny});

        compute_scaled_parameters();
        ditr = std::abs(itr_list[1] - itr_list[0]);
    }
};


// -
