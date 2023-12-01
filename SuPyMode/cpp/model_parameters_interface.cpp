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

}

