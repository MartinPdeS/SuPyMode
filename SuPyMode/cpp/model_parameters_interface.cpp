#include <pybind11/pybind11.h>
#include "includes/model_parameters.cpp"

PYBIND11_MODULE(ModelParameters, module) {
    module.doc() = "A C++ wrapper class for ModelParameters";

    pybind11::class_<ModelParameters>(module, "ModelParameters")
        .def(
            pybind11::init<
                const pybind11::array_t<double>&,
                const pybind11::array_t<double>&,
                const pybind11::array_t<double>&,
                const pybind11::array_t<double>&,
                const double,
                const double,
                const double,
                const int>(),
            pybind11::arg("mesh"),
            pybind11::arg("x_vector"),
            pybind11::arg("y_vector"),
            pybind11::arg("itr_list"),
            pybind11::arg("wavelength"),
            pybind11::arg("dx"),
            pybind11::arg("dy"),
            pybind11::arg("debug_mode") = 0
        )
        .def_readonly("n_slice", &ModelParameters::n_slice)
        .def_readonly("nx", &ModelParameters::nx)
        .def_readonly("ny", &ModelParameters::ny)
        .def_readonly("dx", &ModelParameters::dx)
        .def_readonly("dy", &ModelParameters::dy)
        .def_readonly("wavelength", &ModelParameters::wavelength)
        .def_readonly("wavenumber", &ModelParameters::wavenumber)
        .def_readonly("itr_list", &ModelParameters::itr_list_py)
        .def_readonly("mesh", &ModelParameters::mesh_py)
        .def_readonly("mesh_gradient", &ModelParameters::mesh_gradient_py)
        .def(
            pybind11::pickle(
                &ModelParameters::get_pickle,
                &ModelParameters::build_from_tuple
            )
        );
}
