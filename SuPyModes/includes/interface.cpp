#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Sparse>
#include "/home/martth/temporary/gnuplot-iostream/gnuplot-iostream.h"
#include <Spectra/GenEigsRealShiftSolver.h>

#include <Spectra/MatOp/SparseSymShiftSolve.h>

#include <Spectra/MatOp/SparseGenRealShiftSolve.h>
#include <cstdio>
#include <ctime>
#include <vector>

using namespace Eigen;
using namespace std;
using namespace Spectra;
namespace py = pybind11;

#define PI 3.1415926535897932384626f

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
    .def(py::init<ndarray&, ndarray&, uint, uint, float, bool>(),
         py::arg("Mesh"),
         py::arg("Gradient"),
         py::arg("nMode"),
         py::arg("MaxIter"),
         py::arg("Tolerance"),
         py::arg("debug") )

     .def("LoopOverITR", &EigenSolving::LoopOverITR, py::arg("ITR"), py::arg("alpha"), py::arg("ExtrapolationOrder"))

     .def("ComputingOverlap", &EigenSolving::ComputingOverlap)

     .def("ComputingAdiabatic", &EigenSolving::ComputingAdiabatic)

     .def("ComputingCoupling", &EigenSolving::ComputingCoupling)

     .def("SortModesFields", &EigenSolving::SortModesFields)

     .def("SortModesIndex", &EigenSolving::SortModesIndex)

     .def("ComputingIndices", &EigenSolving::ComputingIndices)

     .def("ComputeLaplacian", &EigenSolving::ComputeLaplacian)

     .def("PringLaplacian", &EigenSolving::PringLaplacian)

     .def("GetSlice", &EigenSolving::GetSlice, py::arg("slice")  = 0)

     .def_property("LeftSymmetry", &EigenSolving::GetLeftSymmetry, &EigenSolving::SetLeftSymmetry)
     .def_property("RightSymmetry", &EigenSolving::GetRightSymmetry, &EigenSolving::SetRightSymmetry)
     .def_property("TopSymmetry", &EigenSolving::GetTopSymmetry, &EigenSolving::SetTopSymmetry)
     .def_property("BottomSymmetry", &EigenSolving::GetBottomSymmetry, &EigenSolving::SetBottomSymmetry)

     .def_property("dx", &EigenSolving::Getdx, &EigenSolving::Setdx)

     .def_property("dy", &EigenSolving::Getdy, &EigenSolving::Setdy)

     .def_property("Lambda", &EigenSolving::Getlambda, &EigenSolving::Setlambda);

}
