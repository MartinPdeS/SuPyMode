#include "mesh.h"

PYBIND11_MODULE(Example, module) {
    module.def(
        "get_rho_gradient_5p",
        &get_rho_gradient_py,
        pybind11::arg("mesh"),
        pybind11::arg("x_vector"),
        pybind11::arg("y_vector")
    );
    module.def(
        "get_rho_gradient_time_rho_5p",
        &get_rho_gradient_time_rho_py,
        pybind11::arg("mesh"),
        pybind11::arg("x_vector"),
        pybind11::arg("y_vector")
    );
}
