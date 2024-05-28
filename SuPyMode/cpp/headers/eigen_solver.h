#pragma once

#include "power_shift_iteration.cpp"

Eigen::VectorXd
pybind11_to_eigen_vector_mapping(const pybind11::array_t<double> &pybind11_array)
{
    size_t nx = pybind11_array.request().shape[0],
           ny = pybind11_array.request().shape[1];

    double *pybind11_array_ptr = (double*) pybind11_array.request().ptr;

    Eigen::Map<Eigen::VectorXd> mapping_gradient(pybind11_array_ptr, nx * ny);

    return mapping_gradient;
}

class EigenSolver
{
public:
    double itr_initial;
    double itr_final;
    double wavenumber;
    double tolerance;
    double initial_value_guess;
    double dx;
    double dy;
    size_t n_slice;
    size_t n_computed_mode;
    size_t n_sorted_mode;
    size_t max_iteration;
    Eigen::VectorXd mesh, gradient, itr_list;
    Eigen::MatrixXd eigen_values, eigen_vectors;
    Eigen::SparseMatrix<double, Eigen::ColMajor> laplacian, eigen_matrix;

    EigenSolver(double itr_initial, double itr_final, double wavenumber, double tolerance, double initial_value_guess,
        size_t n_slice, size_t n_computed_mode, size_t n_sorted_mode, size_t max_iteration,
        const std::vector<std::vector<double>> &finit_difference_matrix,
        const pybind11::array_t<double> &mesh_py,
        const pybind11::array_t<double> &gradient_py,
        const double dx, const double dy)
        : itr_initial(itr_initial), itr_final(itr_final),
          wavenumber(wavenumber), tolerance(tolerance),
          initial_value_guess(initial_value_guess), n_slice(n_slice),
          n_computed_mode(n_computed_mode),
          n_sorted_mode(n_sorted_mode), max_iteration(max_iteration), dx(dx), dy(dy)
         {
            this->mesh = pybind11_to_eigen_vector_mapping(mesh_py);
            this->gradient = pybind11_to_eigen_vector_mapping(gradient_py);
            this->itr_list = Eigen::VectorXd::LinSpaced(n_slice, itr_initial, itr_final);
            this->eigen_values.resize(n_sorted_mode, n_slice);
            this->eigen_vectors.resize(n_sorted_mode, mesh.size());

            this->generate_laplacian_from_triplet(finit_difference_matrix);

            this->compute_first_slice();
         }


    void compute_first_slice();

    void generate_laplacian_from_triplet(const std::vector<std::vector<double>> &finit_difference_matrix);

    Eigen::SparseMatrix<double, Eigen::ColMajor> generate_eigen_matrix(double wavenumber) const;

    void iterate_over_same_slice(const Eigen::VectorXd &initial_vector_guess, double initial_value_guess);
};


