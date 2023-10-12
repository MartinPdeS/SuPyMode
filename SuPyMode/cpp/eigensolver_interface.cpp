#include <pybind11/pybind11.h>
#include "includes/eigensolver.cpp"



PYBIND11_MODULE(CppSolver, module)
{
    pybind11::class_<CppSolver>(module, "CppSolver")
    .def(pybind11::init<ndarray&, ndarray&, Vecf2D&, ndarray&, size_t, size_t, size_t, double, double, double, double, bool, bool>(),
         pybind11::arg("mesh"),
         pybind11::arg("gradient"),
         pybind11::arg("finit_matrix"),
         pybind11::arg("itr_list"),
         pybind11::arg("n_computed_mode"),
         pybind11::arg("n_sorted_mode"),
         pybind11::arg("max_iter"),
         pybind11::arg("tolerance"),
         pybind11::arg("wavelength"),
         pybind11::arg("dx"),
         pybind11::arg("dy"),
         pybind11::arg("show_iteration") = false,
         pybind11::arg("show_eigenvalues") = false
       )

     .def("loop_over_itr",             &CppSolver::loop_over_itr, pybind11::arg("extrapolation_order"), pybind11::arg("alpha"))
     .def("compute_laplacian",         &CppSolver::compute_laplacian)
     .def("get_mode",                  &CppSolver::get_sorted_mode)

     .def_readwrite("n_sorted_mode",   &CppSolver::n_sorted_mode)
     .def_readwrite("n_computed_mode", &CppSolver::n_computed_mode)
     .def_readwrite("alpha_vector",    &CppSolver::alpha_vector)
     .def_readwrite("wavelength",      &CppSolver::wavelength);
}

