#pragma once

    #include "power_shift_iteration.cpp"
    // #include "ritz_iteration.cpp"

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
        double itr_initial, itr_final, wavenumber, tolerance, initial_value_guess, dx, dy;
        size_t n_slice, n_computed_mode, n_sorted_mode, max_iteration;
        Eigen::VectorXd mesh, gradient, itr_list;
        Eigen::MatrixXd eigen_values, eigen_vectors;
        Eigen::SparseMatrix<double, Eigen::ColMajor> laplacian, eigen_matrix;

        EigenSolver(
            const double &itr_initial,
            const double &itr_final,
            const double &wavenumber,
            const double &tolerance,
            const double initial_value_guess,
            const size_t &n_slice,
            const size_t n_computed_mode,
            const size_t n_sorted_mode,
            const size_t max_iteration,
            const std::vector<std::vector<double>> &finit_difference_matrix,
            const pybind11::array_t<double> &mesh_py,
            const pybind11::array_t<double> &gradient_py,
            const double dx,
            const double dy)
            : itr_initial(itr_initial),
              itr_final(itr_final),
              wavenumber(wavenumber),
              tolerance(tolerance),
              initial_value_guess(initial_value_guess),
              n_slice(n_slice),
              n_computed_mode(n_computed_mode),
              n_sorted_mode(n_sorted_mode),
              max_iteration(max_iteration),
              dx(dx),
              dy(dy)
             {
                this->mesh = pybind11_to_eigen_vector_mapping(mesh_py);
                this->gradient = pybind11_to_eigen_vector_mapping(gradient_py);
                this->itr_list = Eigen::VectorXd::LinSpaced(n_slice, itr_initial, itr_final);
                this->eigen_values.resize(n_sorted_mode, n_slice);
                this->eigen_vectors.resize(n_sorted_mode, mesh.size());

                this->generate_laplacian_from_triplet(finit_difference_matrix);

                this->compute_first_slice();
             }


        void generate_laplacian_from_triplet(const std::vector<std::vector<double>> &finit_difference_matrix)
        {
            std::vector<double> row  = finit_difference_matrix[0],
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

        Eigen::SparseMatrix<double, Eigen::ColMajor> generate_eigen_matrix(const double &wavenumber)
        {
            Eigen::SparseMatrix<double, Eigen::ColMajor> eigen_matrix = laplacian;

            Eigen::VectorXd mesh_squared = (mesh * wavenumber).cwiseAbs2();

            eigen_matrix.diagonal() += mesh_squared;

            eigen_matrix *= -1.0;

            return eigen_matrix;
        }

        void compute_first_slice()
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


        void iterate_over_same_slice(const Eigen::VectorXd &initial_vector_guess, const double &initial_value_guess)
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
    };


