#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>
#include "../../extern/eigen/Eigen/Eigenvalues"
#include "../../extern/eigen/Eigen/Sparse"
#include "../../extern/eigen/Eigen/Core"
#include "../../extern/spectra/include/Spectra/GenEigsRealShiftSolver.h"
#include "../../extern/spectra/include/Spectra/MatOp/SparseGenRealShiftSolve.h"
#include <vector>
#include <iostream>



using namespace Eigen;
using namespace std;
using namespace Spectra;
namespace py = pybind11;

#define PI 3.1415926535897932384626f

typedef double                                 ScalarType;
typedef std::complex<ScalarType>               ComplexScalarType;
typedef Matrix<ScalarType, Dynamic, 1>         VectorType;
typedef Matrix<ComplexScalarType, Dynamic, 1>  ComplexVectorType;
typedef Matrix<ScalarType, Dynamic, Dynamic>   MatrixType;
typedef py::array_t<ScalarType>                ndarray;
typedef py::array_t<ComplexScalarType>         Cndarray;
typedef py::buffer_info                        info;
typedef SparseMatrix<ScalarType, ColMajor>     MSparse;
typedef vector<ScalarType>                     Vecf1D;
typedef vector<vector<ScalarType>>             Vecf2D;
typedef Eigen::Triplet<ScalarType> T;

ScalarType inf = numeric_limits<ScalarType>::infinity();

#include "utils.cpp"
#include "class.cpp"



PYBIND11_MODULE(EigenSolver, module) {
    module.doc() = "A c++ solver for EigenPairs";

    py::class_<EigenSolving>(module, "EigenSolving")
    .def(py::init<ndarray&, ndarray&, size_t, size_t, size_t, ScalarType, bool>(),
         py::arg("Mesh"),
         py::arg("Gradient"),
         py::arg("nMode"),
         py::arg("sMode"),
         py::arg("MaxIter"),
         py::arg("Tolerance"),
         py::arg("Debug")=false
       )

     .def("LoopOverITR", &EigenSolving::LoopOverITR, py::arg("ITR"), py::arg("ExtrapolationOrder"))

     .def("LoopOverITR_", &EigenSolving::LoopOverITR_, py::arg("ITR"), py::arg("ExtrapolationOrder"))

     .def("ComputingOverlap", &EigenSolving::ComputingOverlap)

     .def("ComputingAdiabatic", &EigenSolving::ComputingAdiabatic)

     .def("ComputingCoupling", &EigenSolving::ComputingCoupling)

     .def("SortModesFields", &EigenSolving::SortModesFields)

     .def("SortModesIndex", &EigenSolving::SortModesIndex)

     .def("ComputeLaplacian", &EigenSolving::ComputeLaplacian, py::arg("Order"))

     .def("GetSlice", &EigenSolving::GetSlice, py::arg("slice"))
     .def("GetFields", &EigenSolving::GetFields, py::arg("slice"))
     .def("GetIndices", &EigenSolving::GetIndices)
     .def("GetBetas", &EigenSolving::GetBetas)

     .def_property("LeftSymmetry", &EigenSolving::GetLeftSymmetry, &EigenSolving::SetLeftSymmetry)
     .def_property("RightSymmetry", &EigenSolving::GetRightSymmetry, &EigenSolving::SetRightSymmetry)
     .def_property("TopSymmetry", &EigenSolving::GetTopSymmetry, &EigenSolving::SetTopSymmetry)
     .def_property("BottomSymmetry", &EigenSolving::GetBottomSymmetry, &EigenSolving::SetBottomSymmetry)

     .def_property("dx", &EigenSolving::Getdx, &EigenSolving::Setdx)

     .def_property("dy", &EigenSolving::Getdy, &EigenSolving::Setdy)

     .def_property("Lambda", &EigenSolving::Getlambda, &EigenSolving::Setlambda);

}
