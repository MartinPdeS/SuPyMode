#pragma once

#include "definitions.cpp"

class BaseLaplacian{

public:
    pybind11::array_t<double> mesh_py;
    size_t n_x, n_y, size;
    double k_taper, *meshPtr;
    Eigen::SparseMatrix<double, Eigen::ColMajor> Laplacian, Identity, eigen_matrix;
    std::vector<std::vector<double>> finit_difference_matrix;

    BaseLaplacian(ndarray&  mesh_py, std::vector<std::vector<double>> &finit_difference_matrix)
    : mesh_py(mesh_py), finit_difference_matrix(finit_difference_matrix)
    {
        this->n_y = mesh_py.request().shape[0];
        this->n_x = mesh_py.request().shape[1];
        this->size = mesh_py.request().size;
        this->Laplacian = Eigen::SparseMatrix<double, Eigen::ColMajor>(this->size, this->size);
        this->meshPtr = (double*) mesh_py.request().ptr;
    }


    void FromTriplets()
    {
        std::vector<double> Row  = finit_difference_matrix[0],
                            Col  = finit_difference_matrix[1],
                            Data = finit_difference_matrix[2];

        std::vector<fTriplet> Tri;
        Tri.reserve(Row.size());

        for (int i=0; i<Row.size(); i++)
            Tri.push_back(fTriplet(Col[i], Row[i], Data[i]));


        Laplacian.setFromTriplets(Tri.begin(), Tri.end());
    }

    void compute_laplacian()
    {
        Identity = MSparse(size,size); Identity.setIdentity();

        FromTriplets();
    }


    void compute_finit_diff_matrix()
    {
        this->eigen_matrix = Laplacian;

        for (size_t index=0; index<this->size; ++index)
            Identity.coeffRef(index, index) = + pow(meshPtr[index] * k_taper, 2);

        this->eigen_matrix += Identity;

        this->eigen_matrix *= -1.0;
    }

};
