#include "mesh.h"

PYBIND11_MODULE(interface_mesh, module) {
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

    pybind11::class_<Geometry>(module, "GEOMETRY")
        .def(
            pybind11::init<const pybind11::array_t<double>&, const pybind11::array_t<double>&, const pybind11::array_t<double>&>(),
            pybind11::arg("mesh"),
            pybind11::arg("x"),
            pybind11::arg("y")
        )
        .def("compute_gradient_2p", &Geometry::compute_gradient_2p)
        .def("compute_gradient_5p", &Geometry::compute_gradient_5p)
        .def("compute_gradient_7p", &Geometry::compute_gradient_7p)
        .def("get_rho_gradient_time_rho", &Geometry::get_rho_gradient_time_rho)
        .def("get_rho_gradient", &Geometry::get_rho_gradient);
}
