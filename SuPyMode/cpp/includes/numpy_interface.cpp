#pragma once

#include "definitions.cpp"


template<typename T>
std::vector<size_t>
get_stride_from_dimension(std::vector<size_t> dimension)
{
  std::reverse(dimension.begin(), dimension.end());

  std::vector<size_t> stride;
  stride.push_back( sizeof(T) );

  for (size_t i=0; i<dimension.size()-1; ++i)
      stride.push_back( stride[i] * dimension[i] );

  std::reverse(stride.begin(), stride.end());

  return stride;
}




template<typename T>
pybind11::array_t<T> templated_eigen_to_ndarray(
    Eigen::Matrix<T, Eigen::Dynamic, 1> &&eigen3_vector,
    std::vector<size_t> dimension)
{
    Eigen::Matrix<T, Eigen::Dynamic, 1> * vector_pointer = new Eigen::Matrix<T, Eigen::Dynamic, 1>;
    (*vector_pointer) = eigen3_vector;

    std::vector<size_t> stride = get_stride_from_dimension<T>(dimension);

    pybind11::capsule free_when_done(
        vector_pointer->data(), [](void *f)
        {
            T *foo = reinterpret_cast<T *>(f);
            delete []foo;
        }
    );

    pybind11::array_t<T> numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        vector_pointer->data(),
        free_when_done
    );

    return numpy_array;
}


template<typename T>
pybind11::array_t<T> templated_eigen_to_ndarray_copy(
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &eigen3_vector,
    std::vector<size_t> dimension)
{
    Eigen::Matrix<T, Eigen::Dynamic, 1> * vector_pointer = new Eigen::Matrix<T, Eigen::Dynamic, 1>;
    (*vector_pointer) = eigen3_vector;

    std::vector<size_t> stride = get_stride_from_dimension<T>(dimension);

    pybind11::capsule free_when_done(
        vector_pointer->data(), [](void *f)
        {
            T *foo = reinterpret_cast<T *>(f);
            delete []foo;
        }
    );

    pybind11::array_t<T> numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        vector_pointer->data(),
        free_when_done
    );

    return numpy_array;
}

template<typename T>
pybind11::array_t<T> templated_eigen_to_ndarray(
    Eigen::Matrix<T, Eigen::Dynamic, 1> &eigen3_vector,
    std::vector<size_t> dimension)
{
    Eigen::Matrix<T, Eigen::Dynamic, 1> * vector_pointer = new Eigen::Matrix<T, Eigen::Dynamic, 1>;
    (*vector_pointer) = eigen3_vector;

    std::vector<size_t> stride = get_stride_from_dimension<T>(dimension);

    pybind11::capsule free_when_done(
        vector_pointer->data(), [](void *f)
        {
            T *foo = reinterpret_cast<T *>(f);
            delete []foo;
        }
    );

    pybind11::array_t<T> numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        vector_pointer->data(),
        free_when_done
    );

    return numpy_array;
}


template<typename T>
pybind11::array_t<T> templated_eigen_to_ndarray(
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &&eigen3_vector,
    std::vector<size_t> dimension)
{
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> * vector_pointer = new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    (*vector_pointer) = eigen3_vector;

    std::vector<size_t> stride = get_stride_from_dimension<T>(dimension);

    pybind11::capsule free_when_done(
        vector_pointer->data(), [](void *f)
        {
            T *foo = reinterpret_cast<T *>(f);
            delete []foo;
        }
    );

    pybind11::array_t<T> numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        vector_pointer->data(),
        free_when_done
    );

    return numpy_array;
}

template<typename T>
pybind11::array_t<T> templated_eigen_to_ndarray(
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> &eigen3_vector,
    std::vector<size_t> dimension)
{
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> * vector_pointer = new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    (*vector_pointer) = eigen3_vector;

    std::vector<size_t> stride = get_stride_from_dimension<T>(dimension);

    pybind11::capsule free_when_done(
        vector_pointer->data(), [](void *f)
        {
            T *foo = reinterpret_cast<T *>(f);
            delete []foo;
        }
    );

    pybind11::array_t<T> numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        vector_pointer->data(),
        free_when_done
    );

    return numpy_array;
}



template<typename T>
Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>
templated_ndarray_to_eigen_copy(pybind11::array_t<T> array)
{
    double *p = (double*) array.request().ptr;

    size_t nx = array.request().shape[0],
           ny = array.request().shape[1];

    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> mapping(p, ny, nx);

    return mapping;
}

template<typename T>
Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>
templated_ndarray_to_eigen(pybind11::array_t<T> &array)
{
    double *p = (double*) array.request().ptr;

    size_t nx = array.request().shape[0],
           ny = array.request().shape[1];

    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> mapping(p, ny, nx);

    return mapping;
}


template<typename T>
Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>
templated_ndarray_to_eigen_vector(pybind11::array_t<T> &array)
{
    double *p = (double*) array.request().ptr;

    size_t nx = array.request().shape[0];

    Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>> mapping(p, nx);

    return mapping;
}


template<typename T>
pybind11::array_t<T> templated_eigen_to_ndarray_copy(
    const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigen3_vector,
    std::vector<size_t> dimension)
{
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> * vector_pointer = new Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    (*vector_pointer) = eigen3_vector;

    std::vector<size_t> stride = get_stride_from_dimension<T>(dimension);

    pybind11::capsule free_when_done(
        vector_pointer->data(), [](void *f)
        {
            T *foo = reinterpret_cast<T *>(f);
            delete []foo;
        }
    );

    pybind11::array_t<T> numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        vector_pointer->data(),
        free_when_done
    );

    return numpy_array;
}