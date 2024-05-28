#pragma once


#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "supermode.cpp"
#include "extrapolate.cpp"
#include "progressbar.cpp"
#include "utils.cpp"
#include "numpy_interface.cpp"
#include <iostream>

class CppSolver
{
    public:
        size_t n_computed_mode;
        size_t n_sorted_mode;
        size_t max_iteration;
        size_t iteration;

        double tolerance;

        std::vector<SuperMode> computed_supermodes;
        std::vector<SuperMode> sorted_supermodes;

        Eigen::MatrixXd previous_eigen_vectors;
        Eigen::VectorXd previous_eigen_values;

        std::vector<double> alpha_vector;

        ModelParameters model_parameters;

        Eigen::SparseMatrix<double, Eigen::ColMajor> laplacian_matrix;
        Eigen::SparseMatrix<double, Eigen::ColMajor> identity_matrix;
        Eigen::SparseMatrix<double, Eigen::ColMajor> eigen_matrix;
        Eigen::MatrixXd finit_difference_triplets;

        double k_taper;

    CppSolver(
        const ModelParameters &model_parameters,
        const pybind11::array_t<double> &finit_difference_triplets_py,
        size_t n_computed_mode,
        size_t n_sorted_mode,
        size_t max_iteration,
        double tolerance
    )
        : model_parameters(model_parameters),
          n_computed_mode(n_computed_mode),
          n_sorted_mode(n_sorted_mode),
          max_iteration(max_iteration),
          tolerance(tolerance),
          iteration(iteration),
          k_taper(k_taper)
    {
        this->finit_difference_triplets = convert_py_to_eigen<double>(
            finit_difference_triplets_py,
            finit_difference_triplets_py.request().shape[0],
            finit_difference_triplets_py.request().shape[1]
        );

        this->generate_mode_set();
    }

    void compute_laplacian();

    void compute_finit_diff_matrix();

    SuperMode &get_sorted_mode(size_t mode_number){ return sorted_supermodes[mode_number]; }

    void generate_mode_set();

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> compute_eigen_solution(const double alpha);

    void sort_eigen_decomposition(Eigen::MatrixXd &eigen_vectors, Eigen::VectorXd &eigen_values);

    void sort_eigen_with_fields(Eigen::MatrixXd &eigen_vectors, Eigen::VectorXd &eigen_values);

    void populate_sorted_supermodes(size_t slice, Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values);

    void loop_over_itr(size_t extrapolation_order = 1, double alpha = 0);

    void sort_mode_per_last_propagation_constant();

    void arrange_mode_field();

    void normalize_mode_field();

    double compute_max_index();
};
