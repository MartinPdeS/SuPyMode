#pragma once

#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "../../../extern/eigen/Eigen/Sparse"
#include "../../../extern/spectra/include/Spectra/GenEigsRealShiftSolver.h"
#include "../../../extern/spectra/include/Spectra/MatOp/SparseGenRealShiftSolve.h"
#include <vector>
#include <iostream>

#define PYBIND11_DETAILED_ERROR_MESSAGES

typedef Eigen::Matrix<double, Eigen::Dynamic, 1>                VectorType;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>  ComplexVectorType;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>   MatrixType;
typedef pybind11::array_t<double>                               ndarray;
typedef pybind11::array_t<std::complex<double> >                Cndarray;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor>            MSparse;
typedef std::vector<double>                                     Vecf1D;
typedef std::vector<std::vector<double>>                        Vecf2D;
typedef Eigen::Triplet<double>                                  fTriplet;

double inf = std::numeric_limits<double>::infinity();

std::complex<double> J(0.0, 1.0);

#define PI 3.1415926535897932384626f
