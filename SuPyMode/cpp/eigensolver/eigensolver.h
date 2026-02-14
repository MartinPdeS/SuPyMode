#pragma once

#include "Spectra/GenEigsRealShiftSolver.h"
#include "Spectra/MatOp/SparseGenRealShiftSolve.h"
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "../mesh/mesh.h"
#include "../supermode/supermode.h"
#include "../utils/utils.h"
#include "../utils/extrapolator.h"
#include "../utils/progress_bar.h"
#include "../utils/numpy_interface.h"
#include "boundaries/boundaries.h"
#include <iostream>

class EigenSolver
{
    public:
        ModelParameters model_parameters;
        size_t n_computed_mode;
        size_t n_sorted_mode;
        size_t max_iteration;
        double tolerance;
        size_t accuracy;
        Boundaries boundaries;
        size_t iteration;
        double wavelength;


        std::vector<SuperMode> computed_supermodes;
        std::vector<SuperMode> sorted_supermodes;

        Eigen::MatrixXd previous_eigen_vectors;
        Eigen::VectorXd previous_eigen_values;

        std::vector<double> alpha_vector;

    public:
        Eigen::SparseMatrix<double, Eigen::ColMajor> laplacian_matrix;
        Eigen::SparseMatrix<double, Eigen::ColMajor> identity_matrix;
        Eigen::SparseMatrix<double, Eigen::ColMajor> eigen_matrix;
        Eigen::MatrixXd finit_difference_triplets;

        double k_taper;

    EigenSolver(
        const size_t max_iteration,
        const double tolerance,
        const size_t accuracy
    );

    void initialize(
        const ModelParameters &model_parameters,
        const pybind11::array_t<double> &finit_difference_triplets_py,
        const size_t n_computed_mode,
        const size_t n_sorted_mode);

    void reset_solver() {
        if (this->model_parameters.debug_mode > 0)
            std::cout << "Resetting the solver to empty state." << std::endl;

        iteration = 0;
        computed_supermodes.clear();
        sorted_supermodes.clear();
        previous_eigen_vectors.setZero();
        previous_eigen_values.setZero();
        alpha_vector.clear();
    }


    /**
     * \brief Sets up the boundary conditions for the eigenvalue solver.
     *
     * \param left The boundary condition for the left edge.
     * \param right The boundary condition for the right edge.
     * \param top The boundary condition for the top edge.
     * \param bottom The boundary condition for the bottom edge.
     */
    void setup_boundaries(const std::string &left, const std::string &right, const std::string &top, const std::string &bottom);


    /**
     * \brief Computes the Laplacian matrix for the given model parameters.
     */
    void compute_laplacian();

    /**
     * \brief Computes the finite difference matrix for the given model parameters.
     */
    void compute_finit_diff_matrix();

    /**
     * \brief Computes the eigenvalues and eigenvectors of the Laplacian matrix.
     *
     * \param alpha The alpha parameter for the eigenvalue problem.
     * \return A tuple containing the eigenvalues and eigenvectors.
     */
    SuperMode &get_sorted_mode(size_t mode_number){ return sorted_supermodes[mode_number]; }

    /**
     * \brief Generates the mode set based on the model parameters.
     */
    void generate_mode_set();

    /**
     * \brief Computes the eigenvalues and eigenvectors of the Laplacian matrix.
     *
     * \param alpha The alpha parameter for the eigenvalue problem.
     * \return A tuple containing the eigenvalues and eigenvectors.
     */
    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> compute_eigen_solution(const double alpha);

    /**
     * \brief Computes the eigenvalues and eigenvectors of the Laplacian matrix.
     *
     * \param alpha The alpha parameter for the eigenvalue problem.
     * \return A tuple containing the eigenvalues and eigenvectors.
     */
    void sort_eigen_decomposition(Eigen::MatrixXd &eigen_vectors, Eigen::VectorXd &eigen_values);

    /**
     * \brief Computes the eigenvalues and eigenvectors of the Laplacian matrix.
     *
     * \param alpha The alpha parameter for the eigenvalue problem.
     * \return A tuple containing the eigenvalues and eigenvectors.
     */
    void sort_eigen_with_fields(Eigen::MatrixXd &eigen_vectors, Eigen::VectorXd &eigen_values);

    /**
     * \brief Populates the sorted supermodes with the computed eigenvalues and eigenvectors.
     *
     * \param slice The slice index for the supermode.
     * \param eigen_vectors The computed eigenvectors.
     * \param eigen_values The computed eigenvalues.
     */
    void populate_sorted_supermodes(size_t slice, Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values);

    /**
     * \brief Loops over the iterations to compute the eigenvalues and eigenvectors.
     *
     * \param extrapolation_order The order of extrapolation for the eigenvalue problem.
     * \param alpha The alpha parameter for the eigenvalue problem.
     */
    void loop_over_itr(size_t extrapolation_order = 1, double alpha = 0);

    /**
     * \brief Sorts the modes based on the last propagation constant.
     */
    void sort_mode_per_last_propagation_constant();

    /**
     * \brief Arranges the mode field.
     */
    void arrange_mode_field();

    /**
     * \brief Normalizes the mode field.
     */
    void normalize_mode_field();

    /**
     * \brief Computes the maximum index for the eigenvalue problem.
     *
     * \return The maximum index.
     */
    double compute_max_index();
};
