#include <pybind11/pybind11.h>
#include "includes/eigensolver.cpp"



PYBIND11_MODULE(CppSolver, module)
{
    pybind11::class_<CppSolver>(module, "CppSolver")
    .def(
        pybind11::init<
            const ModelParameters&,
            const pybind11::array_t<double>&,
            const size_t,
            const size_t,
            const size_t,
            const double,
            const std::string&,
            const std::string&,
            const std::string&,
            const std::string&>(),
        pybind11::arg("model_parameters"),
        pybind11::arg("finit_matrix"),
        pybind11::arg("n_computed_mode"),
        pybind11::arg("n_sorted_mode"),
        pybind11::arg("max_iter"),
        pybind11::arg("tolerance"),
        pybind11::arg("left_boundary"),
        pybind11::arg("right_boundary"),
        pybind11::arg("top_boundary"),
        pybind11::arg("bottom_boundary")
    )

    .def("loop_over_itr", &CppSolver::loop_over_itr, pybind11::arg("extrapolation_order"), pybind11::arg("alpha"))
    .def("compute_laplacian", &CppSolver::compute_laplacian)
    .def("get_mode", &CppSolver::get_sorted_mode)

    .def_readwrite("n_sorted_mode", &CppSolver::n_sorted_mode)
    .def_readwrite("n_computed_mode", &CppSolver::n_computed_mode)
    .def_readwrite("alpha_vector", &CppSolver::alpha_vector)

    .def_readonly("left_boundary", &CppSolver::left_boundary)
    .def_readonly("right_boundary", &CppSolver::right_boundary)
    .def_readonly("top_boundary", &CppSolver::top_boundary)
    .def_readonly("bottom_boundary", &CppSolver::bottom_boundary);
}

