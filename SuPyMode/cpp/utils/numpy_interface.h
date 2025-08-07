#pragma once

#include <vector>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <Eigen/Core>

namespace numy_interface {

template<typename T>
std::vector<size_t> get_stride_from_dimension(std::vector<size_t> dimension) {
  std::reverse(dimension.begin(), dimension.end());

  std::vector<size_t> stride;
  stride.push_back( sizeof(T) );

  for (size_t i=0; i<dimension.size()-1; ++i)
      stride.push_back( stride[i] * dimension[i] );

  std::reverse(stride.begin(), stride.end());

  return stride;
}

template<typename T, typename MatrixType>
pybind11::array_t<T> eigen_to_ndarray(const MatrixType &eigen_matrix, const std::vector<size_t> &dimension) {
    // Create a dynamic copy of the input Eigen matrix
    MatrixType* matrix_pointer = new MatrixType;
    *matrix_pointer = eigen_matrix;

    // Calculate the stride from the given dimensions
    std::vector<size_t> stride = get_stride_from_dimension<T>(dimension);

    // Define a capsule to manage memory. This lambda ensures the allocated matrix is properly deleted
    pybind11::capsule free_when_done(
        matrix_pointer->data(), [](void *f)
        {
            T* data_pointer = reinterpret_cast<T*>(f);
            delete[] data_pointer; // Make sure to delete the matrix data, not the matrix object itself
        }
    );

    // Create a numpy array using the dimensions, stride, and data pointer
    pybind11::array_t<T> numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        matrix_pointer->data(),
        free_when_done
    );

    return numpy_array;
}

template<typename T>
Eigen::MatrixXd convert_py_to_eigen(const pybind11::array_t<T> &array_py) {
    auto info = array_py.request();
    T* ptr = static_cast<T*>(info.ptr);
    return Eigen::Map<Eigen::MatrixXd>(ptr, info.shape[0], info.shape[1]);
}

template<typename T>
Eigen::MatrixXd convert_py_to_eigen(const pybind11::array_t<T> &array_py, const size_t rows, const size_t cols) {
    auto info = array_py.request();
    T* ptr = static_cast<T*>(info.ptr);
    return Eigen::Map<Eigen::MatrixXd>(ptr, rows, cols);
}

template<typename T>
Eigen::VectorXd convert_py_to_eigen(const pybind11::array_t<T> &array_py, const size_t rows) {
    auto info = array_py.request();
    T* ptr = static_cast<T*>(info.ptr);
    return Eigen::Map<Eigen::VectorXd>(ptr, rows);
}

}