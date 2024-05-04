#pragma once

#include "supermode.cpp"
#include "extrapolate.cpp"
#include "progressbar.cpp"
#include "laplacian.cpp"
#include "definitions.cpp"
#include "utils.cpp"
#include "numpy_interface.cpp"





class CppSolver : public BaseLaplacian
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

    CppSolver(
        pybind11::array_t<double> &mesh_py,
        pybind11::array_t<double> &mesh_gradient_py,
        std::vector<std::vector<double>> &finit_difference_matrix,
        pybind11::array_t<double> &itr_list_py,
        size_t n_computed_mode,
        size_t n_sorted_mode,
        size_t max_iteration,
        double tolerance,
        double wavelength,
        double dx,
        double dy,
        int debug_mode)
        : BaseLaplacian(mesh_py, finit_difference_matrix),
          max_iteration(max_iteration),
          tolerance(tolerance),
          n_sorted_mode(n_sorted_mode),
          n_computed_mode(n_computed_mode)
    {
        this->model_parameters = ModelParameters(
            wavelength,
            mesh_gradient_py,
            itr_list_py,
            dx,
            dy
        );

        this->model_parameters.debug_mode = debug_mode;

        this->generate_mode_set();
    }

    SuperMode &get_sorted_mode(size_t mode_number){ return sorted_supermodes[mode_number]; }

    void generate_mode_set();

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd> compute_eigen_solution(const double &alpha);

    void sort_eigen_decomposition(Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values);

    void sort_eigen_with_fields(Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values);

    void populate_sorted_supermodes(size_t slice, Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values);

    void loop_over_itr(size_t extrapolation_order = 1, double alpha = 0);

    void sort_mode_per_last_propagation_constant();

    void arrange_mode_field();

    void normalize_mode_field();

    double compute_max_index();
};
