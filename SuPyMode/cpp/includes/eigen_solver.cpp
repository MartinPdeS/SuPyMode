#pragma once

#include "eigen_solver.h"


void EigenSolver::generate_laplacian_from_triplet(const std::vector<std::vector<double>> &finit_difference_matrix)
{
    std::vector<double>
        row  = finit_difference_matrix[0],
        col  = finit_difference_matrix[1],
        data = finit_difference_matrix[2];

    std::vector<Eigen::Triplet<double>> triplet;
    triplet.reserve(row.size());

    for (int i=0; i<row.size(); i++)
        triplet.push_back(
            Eigen::Triplet<double>(col[i], row[i], data[i])
        );


    laplacian.setFromTriplets(triplet.begin(), triplet.end());
}

Eigen::SparseMatrix<double, Eigen::ColMajor> EigenSolver::generate_eigen_matrix(double wavenumber) const
{
    Eigen::SparseMatrix<double, Eigen::ColMajor> eigen_matrix = laplacian;

    Eigen::VectorXd mesh_squared = (mesh * wavenumber).cwiseAbs2();

    eigen_matrix.diagonal() += mesh_squared;

    eigen_matrix *= -1.0;

    return eigen_matrix;
}

void EigenSolver::compute_first_slice()
{
    eigen_matrix = generate_eigen_matrix(wavenumber);

    Spectra::SparseGenRealShiftSolve<double> solver(eigen_matrix);

    Spectra::GenEigsRealShiftSolver<Spectra::SparseGenRealShiftSolve<double>> solution(
        solver,
        n_computed_mode,
        n_computed_mode + 3 ,
        initial_value_guess
    );

    solution.init();

    int nconv = solution.compute(
        Spectra::SortRule::LargestMagn,
        this->max_iteration,
        this->tolerance,
        Spectra::SortRule::LargestMagn
    );

    Eigen::MatrixXd slice_eigen_vectors = solution.eigenvectors().real();
    Eigen::VectorXd slice_eigen_values = solution.eigenvalues().real();

    this->eigen_values.col(0) = slice_eigen_values;
    this->eigen_vectors.col(0) = slice_eigen_vectors;
}


void EigenSolver::iterate_over_same_slice(const Eigen::VectorXd &initial_vector_guess, double initial_value_guess)
{
    double eigen_value;
    Eigen::VectorXd eigen_vector;

    eigen_matrix = generate_eigen_matrix(wavenumber);

    std::tie(eigen_value, eigen_vector) = inverse_power_method(
        eigen_matrix,
        initial_vector_guess,
        initial_value_guess,
        tolerance,
        max_iteration
    );

    return ;
}



