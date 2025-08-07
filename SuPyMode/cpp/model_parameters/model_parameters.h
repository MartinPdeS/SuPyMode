#pragma once

#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "../mesh/mesh.h"
#include "../utils/numpy_interface.h"


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
        const int debug_mode = 0);

    void compute_scaled_parameters();

    /**
     * \brief Returns a tuple that can be used to pickle the ModelParameters object.
     * \return A tuple containing the mesh, x_vector, y_vector, itr_list,
     * wavelength, dx, dy, and debug_mode.
     */
    static pybind11::tuple get_pickle(const ModelParameters &model_parameters);

    /**
     * \brief Builds a ModelParameters object from a tuple.
     * \param tuple A tuple containing the mesh, x_vector, y_vector, itr_list,
     * wavelength, dx, dy, and debug_mode.
     * \return A ModelParameters object constructed from the tuple.
     */
    static ModelParameters build_from_tuple(const pybind11::tuple &tuple);

    /**
     * \brief Overloads the output stream operator for ModelParameters.
     * \param os The output stream.
     * \param params The ModelParameters object.
     * \return The output stream.
     */
    friend std::ostream &operator<<(std::ostream &os, const ModelParameters &params);

private:
    void init();
};


// -
