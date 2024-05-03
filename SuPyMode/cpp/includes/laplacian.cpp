#pragma once

#include "definitions.cpp"

class BaseLaplacian
{

    public:
        pybind11::array_t<double> mesh_py;
        size_t n_x;
        size_t n_y;
        size_t size;
        double k_taper;
        double *mesh_ptr;
        Eigen::SparseMatrix<double, Eigen::ColMajor> Laplacian;
        Eigen::SparseMatrix<double, Eigen::ColMajor> identity_matrix;
        Eigen::SparseMatrix<double, Eigen::ColMajor> eigen_matrix;
        std::vector<std::vector<double>> finit_difference_matrix;

        BaseLaplacian(pybind11::array_t<double>&  mesh_py, std::vector<std::vector<double>> &finit_difference_matrix)
        : mesh_py(mesh_py), finit_difference_matrix(finit_difference_matrix)
        {
            this->n_y = mesh_py.request().shape[0];
            this->n_x = mesh_py.request().shape[1];
            this->size = mesh_py.request().size;
            this->Laplacian = Eigen::SparseMatrix<double, Eigen::ColMajor>(this->size, this->size);
            this->mesh_ptr = (double*) mesh_py.request().ptr;
        }


        void build_from_triplet()
        {
            std::vector<double>
                row  = finit_difference_matrix[0],
                column  = finit_difference_matrix[1],
                data = finit_difference_matrix[2];

            std::vector<Eigen::Triplet<double>> triplet;
            triplet.reserve(row.size());

            for (int i = 0; i < row.size(); i++)
                triplet.push_back(Eigen::Triplet<double>(column[i], row[i], data[i]));


            Laplacian.setFromTriplets(triplet.begin(), triplet.end());
        }

        void compute_laplacian()
        {
            this->identity_matrix = Eigen::SparseMatrix<double, Eigen::ColMajor>(size,size);
            this->identity_matrix.setIdentity();

            build_from_triplet();
        }


        void compute_finit_diff_matrix()
        {
            this->eigen_matrix = Laplacian;

            for (size_t index=0; index < this->size; ++index)
                identity_matrix.coeffRef(index, index) = + pow(mesh_ptr[index] * k_taper, 2);

            this->eigen_matrix += identity_matrix;

            this->eigen_matrix *= -1.0;
        }

};
