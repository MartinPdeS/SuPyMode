#pragma once

#include "Eigen/Sparse"
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "numpy_interface.cpp"


#define PI 3.1415926535897932384626f


class ModelParameters
{
    public:
        double wavelength;
        double wavenumber;
        double dx;
        double dy;
        double ditr;

        size_t nx;
        size_t ny;
        size_t n_slice;
        size_t field_size;

        pybind11::array_t<double> itr_list_py;
        pybind11::array_t<double> mesh_gradient_py;

        Eigen::VectorXd mesh_gradient;
        Eigen::VectorXd itr_list;
        Eigen::VectorXd dx_scaled;
        Eigen::VectorXd dy_scaled;
        Eigen::VectorXd wavenumber_scaled;

        int debug_mode;

        ModelParameters() = default;

        ModelParameters(double wavelength, pybind11::array_t<double> mesh_gradient_py, pybind11::array_t<double> itr_list_py, double dx, double dy, int debug_mode = 0)
        : wavelength(wavelength), itr_list_py(itr_list_py), mesh_gradient_py(mesh_gradient_py), dx(dx), dy(dy), debug_mode(debug_mode)
        {
            this->nx = mesh_gradient_py.request().shape[0];
            this->ny = mesh_gradient_py.request().shape[1];

            this->wavenumber = 2 * PI / this->wavelength;
            this->field_size = nx * ny;
            this->n_slice = itr_list_py.size();

            this->itr_list = convert_py_to_eigen<double>(itr_list_py, n_slice);
            this->mesh_gradient = convert_py_to_eigen<double>(mesh_gradient_py, nx * ny);

            this->ditr = abs(itr_list[1] - itr_list[0]);

            this->compute_scaled_parameters();
        }

        void compute_scaled_parameters()
        {
            this->dx_scaled = itr_list * dx;
            this->dy_scaled = itr_list * dy;
            this->wavenumber_scaled = itr_list * this->wavenumber;
        }

        static pybind11::tuple get_pickle(ModelParameters &model_parameters)
        {
            return pybind11::make_tuple(
                model_parameters.wavelength,
                model_parameters.mesh_gradient_py,
                model_parameters.itr_list_py,
                model_parameters.dx,
                model_parameters.dy,
                model_parameters.n_slice
            );
        }

        static ModelParameters build_from_tuple(pybind11::tuple tuple) {
            return ModelParameters{
                tuple[0].cast<double>(),                             // wavelength
                tuple[1].cast<pybind11::array_t<double>>(),          // mesh_gradient_py,
                tuple[2].cast<pybind11::array_t<double>>(),          // itr_list_py,
                tuple[3].cast<double>(),                             // dx
                tuple[4].cast<double>()                              // dy
            }; // load
        }

        std::ostream &operator<<(std::ostream &os)
        {
            return os << 'dx: '<< this->dx
                      << 'dy: '<< this->dy
                      << 'nx: '<< this->nx
                      << 'ny: '<< this->ny
                      << 'n_slice: '<< this->n_slice;
        }

};