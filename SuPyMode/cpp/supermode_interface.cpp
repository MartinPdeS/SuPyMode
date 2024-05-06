#include <pybind11/pybind11.h>
#include "includes/supermode.cpp"


PYBIND11_MODULE(SuperMode, module)
{
    module.doc() = "A c++ wrapper class for SuperMode";

    pybind11::class_<SuperMode>(module, "SuperMode")
    .def_readonly("fields", &SuperMode::fields)
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
    .def(pybind11::pickle(&SuperMode::get_pickle, &SuperMode::build_from_tuple))
    .def("get_field", &SuperMode::get_field_py, pybind11::arg("index"));
}

