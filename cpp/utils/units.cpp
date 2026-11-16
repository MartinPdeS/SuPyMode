#include "units.h"
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

py::object numerical_input(py::handle value, const std::string &unit, const std::string &name,
                           int rank, bool positive, bool finite) {
    py::module_ typedunit = py::module_::import("TypedUnit");
    py::object numeric = py::reinterpret_borrow<py::object>(value);
    if (!unit.empty() && py::isinstance(value, typedunit.attr("Quantity"))) {
        typedunit.attr("validate_quantity")(
            value, py::arg("unit") = unit, py::arg("name") = name,
            py::arg("registry") = typedunit.attr("ureg"), py::arg("positive") = positive,
            py::arg("finite") = finite, py::arg("scalar") = rank == 0);
        numeric = numeric.attr("to")(unit).attr("magnitude");
    }
    if (!unit.empty() && py::module_::import("numpy").attr("iscomplexobj")(numeric).cast<bool>())
        throw py::value_error(name + " must contain real physical coordinates or lengths.");
    std::vector<py::object> markers;
    if (rank >= 0) {
        const char *shapes[] = {"Scalar", "Vector", "Matrix"};
        markers.push_back(typedunit.attr(shapes[rank]));
    }
    if (finite)
        markers.push_back(typedunit.attr("Finite"));
    if (positive)
        markers.push_back(typedunit.attr("Positive"));
    py::tuple marker_tuple(markers.size());
    for (size_t index = 0; index < markers.size(); ++index)
        marker_tuple[index] = markers[index];
    py::object annotation = markers.empty()
                                ? typedunit.attr("NoUnit")
                                : typedunit.attr("NoUnit").attr("__class_getitem__")(marker_tuple);
    typedunit.attr("validate_annotation")(numeric, py::arg("annotation") = annotation,
                                          py::arg("name") = name);
    return numeric;
}
