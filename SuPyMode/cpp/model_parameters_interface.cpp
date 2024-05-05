#include <pybind11/pybind11.h>
#include "includes/model_parameters.cpp"


PYBIND11_MODULE(ModelParameters, module)
{
    module.doc() = "A c++ wrapper class for ModelParameters";

    pybind11::class_<ModelParameters>(module, "ModelParameters")

    .def_readwrite("nx", &ModelParameters::nx)
    .def_readwrite("ny", &ModelParameters::ny)
    .def_readwrite("dx", &ModelParameters::dx)
    .def_readwrite("dy", &ModelParameters::dy)
    .def_readwrite("wavelength", &ModelParameters::wavelength)
    .def_readwrite("wavenumber", &ModelParameters::wavenumber)
    .def_readwrite("itr_list", &ModelParameters::itr_list_py)
    .def(pybind11::pickle(&ModelParameters::get_pickle, &ModelParameters::build_from_tuple));
}

