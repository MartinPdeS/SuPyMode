#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "../../../extern/eigen/Eigen/Sparse"
#include "../../../extern/spectra/include/Spectra/GenEigsRealShiftSolver.h"
#include "../../../extern/spectra/include/Spectra/MatOp/SparseGenRealShiftSolve.h"
#include <vector>
#include <iostream>

#define PYBIND11_DETAILED_ERROR_MESSAGES
typedef double                                                      ScalarType;
typedef std::complex<double>                                    ComplexScalarType;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1>                VectorType;
typedef Eigen::Matrix<ComplexScalarType, Eigen::Dynamic, 1>         ComplexVectorType;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>   MatrixType;
typedef pybind11::array_t<double>                               ndarray;
typedef pybind11::array_t<ComplexScalarType>                        Cndarray;
typedef Eigen::SparseMatrix<double, Eigen::ColMajor>            MSparse;
typedef std::vector<double>                                     Vecf1D;
typedef std::vector<std::vector<double>>                        Vecf2D;
typedef Eigen::Triplet<double>                                  fTriplet;

double inf = std::numeric_limits<double>::infinity();

std::complex<double> J(0.0, 1.0);

#define PI 3.1415926535897932384626f


class ModelParameters
{
    public:

        double dx, dy;
        size_t nx, ny, n_slice, field_size;
        Eigen::VectorXd itr_list, dx_scaled, dy_scaled;
        pybind11::array_t<double> itr_list_py;

        ModelParameters(){}
        ModelParameters(pybind11::array_t<double> itr_list_py, double dx, double dy, size_t nx, size_t ny)
        : dx(dx), dy(dy), nx(nx), ny(ny), itr_list_py(itr_list_py)
        {
            this->field_size = nx * ny;
            this->n_slice = itr_list_py.size();

            double *itr_list_ptr = (double*) itr_list_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_itr_list(itr_list_ptr, n_slice);
            this->itr_list = mapping_itr_list;

            this->dx_scaled = itr_list * dx;
            this->dy_scaled = itr_list * dy;
        }

    pybind11::tuple get_state()
    {
        return pybind11::make_tuple(
            this->itr_list_py,
            this->dx,
            this->dy,
            this->nx,
            this->ny,
            this->n_slice
            );
    }
};