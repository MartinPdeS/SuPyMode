#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Sparse>
#include "/home/martth/temporary/gnuplot-iostream/gnuplot-iostream.h"
#include <Spectra/SymEigsShiftSolver.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <cstdio>
#include <ctime>
#include <vector>

using namespace Eigen;
using namespace std;
using namespace Spectra;
namespace py = pybind11;

typedef std::complex<float>              fComplex;
typedef py::array_t<float>               ndarray;
typedef py::array_t<fComplex>            Cndarray;
typedef py::buffer_info                  info;
typedef unsigned int                     uint;
typedef SparseMatrix<float, ColMajor>    MSparse;
typedef vector<float>                    Vecf1D;
typedef vector<vector<float>>            Vecf2D;
#include "utils.cpp"
#include "class.cpp"



PYBIND11_MODULE(EigenSolver, module) {
    module.doc() = "A c++ solver for EigenPairs";

    py::class_<EigenSolving>(module, "EigenSolving")
    .def(py::init<ndarray&, uint, uint, float>(),
         py::arg("Mesh"),
         py::arg("nMode"),
         py::arg("MaxIter"),
         py::arg("Tolerance") )

     .def("GetMatrix", &EigenSolving::GetMatrix)

     .def("LoopOverITR", &EigenSolving::LoopOverITR)

     .def("ComputingOverlap", &EigenSolving::ComputingOverlap)

     .def("ComputingCoupling", &EigenSolving::ComputingCoupling)

     .def("SortModesFields", &EigenSolving::SortModesFields)

     .def("SortModesIndex", &EigenSolving::SortModesIndex)

     .def("ComputingIndices", &EigenSolving::ComputingIndices)

     .def("GetSlice", &EigenSolving::GetSlice, py::arg("slice")  = 0)

     .def_property("dx", &EigenSolving::Getdx, &EigenSolving::Setdx)

     .def_property("dy", &EigenSolving::Getdy, &EigenSolving::Setdy)

     .def_property("Lambda", &EigenSolving::Getlambda, &EigenSolving::Setlambda);

}
