#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "boundaries.h"

namespace py = pybind11;


PYBIND11_MODULE(interface_boundaries, module) {

    py::enum_<BoundaryCondition>(module, "BoundaryValue")
        .value("Zero", BoundaryCondition::Zero)
        .value("Symmetric", BoundaryCondition::Symmetric)
        .value("AntiSymmetric", BoundaryCondition::AntiSymmetric)
        .export_values();

    py::class_<Boundaries, std::shared_ptr<Boundaries>>(module, "Boundaries")
        .def(
            py::init<
                const BoundaryCondition&,
                const BoundaryCondition&,
                const BoundaryCondition&,
                const BoundaryCondition&
            >(),
            py::arg("left") = BoundaryCondition::Zero,
            py::arg("right") = BoundaryCondition::Zero,
            py::arg("top") = BoundaryCondition::Zero,
            py::arg("bottom") = BoundaryCondition::Zero
        )
        .def_readonly(
            "left",
            &Boundaries::left
        )
        .def_readonly(
            "right",
            &Boundaries::right
        )
        .def_readonly(
            "top",
            &Boundaries::top
        )
        .def_readonly(
            "bottom",
            &Boundaries::bottom
        )
        .def("__eq__", &Boundaries::operator==, py::is_operator());

}
