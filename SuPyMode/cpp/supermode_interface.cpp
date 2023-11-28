#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include "definitions.cpp"
#include "includes/supermode.h"
#include "includes/supermode.cpp"


PYBIND11_MODULE(SuperMode, module)
{
    module.doc() = "A c++ wrapper class for SuperMode";

    pybind11::class_<ModelParameters>(module, "ModelParameters")
    .def(
        pybind11::pickle(
            [](ModelParameters& model_parameter)
            {
                return model_parameter.get_state();  // dump
            },
            [](pybind11::tuple t)
            {
                return ModelParameters{
                    t[0].cast<double>(),                             // wavelength
                    t[1].cast<pybind11::array_t<double>>(),          // mesh_gradient_py,
                    t[2].cast<pybind11::array_t<double>>(),          // itr_list_py,
                    t[3].cast<double>(),                             // dx
                    t[4].cast<double>()                              // dy
                }; // load
            }
        )


    );

    pybind11::class_<SuperMode>(module, "SuperMode")

    .def_readwrite("binding_number", &SuperMode::mode_number)
    .def_readwrite("model_parameters", &SuperMode::model_parameters)

    .def("itr_list", &SuperMode::get_itr_list)

    .def("get_fields", &SuperMode::get_fields_py)
    .def("get_norm", &SuperMode::get_norm)
    .def("get_index", &SuperMode::get_index_py)
    .def("get_betas", &SuperMode::get_betas_py)
    .def("get_eigen_value", &SuperMode::get_eigen_value_py)

    .def("get_normalized_coupling_with_mode", &SuperMode::get_normalized_coupling_with_mode_py)
    .def("get_adiabatic_with_mode", &SuperMode::get_adiabatic_with_mode_py)
    .def("get_overlap_integrals_with_mode", &SuperMode::get_overlap_integrals_with_mode_py)
    .def("get_beating_length_with_mode", &SuperMode::get_beating_length_with_mode_py)
    .def(
        pybind11::pickle(
            [](SuperMode& a)
            {
                return  a.get_state();  // dump
            },
            [](pybind11::tuple t)
            {
                return SuperMode{
                    t[0].cast<size_t>(),                             // mode_number
                    t[1].cast<pybind11::array_t<double>>(),          // fields
                    t[2].cast<pybind11::array_t<double>>(),          // index
                    t[3].cast<pybind11::array_t<double>>(),          // betas
                    t[4].cast<pybind11::array_t<double>>(),          // eigen_values
                    t[5].cast<ModelParameters>()                     // Model parameters
                }; // load
            }
        )


    );
}

