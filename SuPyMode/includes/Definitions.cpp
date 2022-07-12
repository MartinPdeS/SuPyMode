#pragma once


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

std::complex<ScalarType> J(0.0, 1.0);
