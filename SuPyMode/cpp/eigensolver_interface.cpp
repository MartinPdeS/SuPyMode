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
            const double
        >(),
        pybind11::arg("model_parameters"),
        pybind11::arg("finit_matrix"),
        pybind11::arg("n_computed_mode"),
        pybind11::arg("n_sorted_mode"),
        pybind11::arg("max_iter"),
        pybind11::arg("tolerance")
    )

     .def("loop_over_itr", &CppSolver::loop_over_itr, pybind11::arg("extrapolation_order"), pybind11::arg("alpha"))
     .def("compute_laplacian", &CppSolver::compute_laplacian)
     .def("get_mode", &CppSolver::get_sorted_mode)

     .def_readwrite("n_sorted_mode", &CppSolver::n_sorted_mode)
     .def_readwrite("n_computed_mode", &CppSolver::n_computed_mode)
     .def_readwrite("alpha_vector", &CppSolver::alpha_vector);
}

