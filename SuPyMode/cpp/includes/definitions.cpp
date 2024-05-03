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

std::complex<double> J(0.0, 1.0);

#define PI 3.1415926535897932384626f
