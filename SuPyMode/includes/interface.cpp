#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
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

    py::class_<SuperMode>(module, "SuperMode")
    .def("GetFields",                &SuperMode::GetFields)
    .def("GetIndex",                 &SuperMode::GetIndex)
    .def("GetBetas",                 &SuperMode::GetBetas)
    .def("GetCoupling",              &SuperMode::GetCoupling)
    .def("GetAdiabatic",             &SuperMode::GetAdiabatic)
    .def("GetAdiabaticSpecific",     &SuperMode::GetAdiabaticSpecific)
    .def("GetCouplingSpecific",      &SuperMode::GetCouplingSpecific)
    .def_readwrite("LeftSymmetry",   &SuperMode::LeftSymmetry)
    .def_readwrite("RightSymmetry",  &SuperMode::RightSymmetry)
    .def_readwrite("TopSymmetry",    &SuperMode::TopSymmetry)
    .def_readwrite("BottomSymmetry", &SuperMode::BottomSymmetry)
    .def_readwrite("BindingNumber",  &SuperMode::ModeNumber);


    py::class_<EigenSolving>(module, "EigenSolving")
    .def(py::init<ndarray&, ndarray&, size_t, size_t, size_t, ScalarType, ScalarType, ScalarType, ScalarType, bool>(),
         py::arg("Mesh"),
         py::arg("Gradient"),
         py::arg("nMode"),
         py::arg("sMode"),
         py::arg("MaxIter"),
         py::arg("Tolerance"),
         py::arg("dx"),
         py::arg("dy"),
         py::arg("Wavelength"),
         py::arg("Debug") = false
       )

     .def("LoopOverITR",                 &EigenSolving::LoopOverITR, py::arg("ITR"), py::arg("ExtrapolationOrder"))
     .def("ComputeAdiabatic",            &EigenSolving::ComputeAdiabatic)
     .def("ComputeCouplingAdiabatic",    &EigenSolving::ComputeCouplingAdiabatic)
     .def("ComputeCoupling",             &EigenSolving::ComputeCoupling)
     .def("SortModes",                   &EigenSolving::SortModes, py::arg("Sorting"))
     .def("ComputeLaplacian",            &EigenSolving::ComputeLaplacian, py::arg("Order"))
     .def("GetSlice",                    &EigenSolving::GetSlice, py::arg("slice"))


     .def("GetMode",    &EigenSolving::GetMode)


     .def_readwrite("LeftSymmetry",   &EigenSolving::LeftSymmetry)
     .def_readwrite("RightSymmetry",  &EigenSolving::RightSymmetry)
     .def_readwrite("TopSymmetry",    &EigenSolving::TopSymmetry)
     .def_readwrite("BottomSymmetry", &EigenSolving::BottomSymmetry)


     .def_readwrite("sMode",          &EigenSolving::sMode)
     .def_readwrite("nMode",          &EigenSolving::nMode)
     .def_readwrite("ITRLength",      &EigenSolving::ITRLength)
     .def_readwrite("ITRList",        &EigenSolving::ITRList)

     .def_readwrite("Wavelength",     &EigenSolving::lambda);

}
