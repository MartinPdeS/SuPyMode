#pragma once

    // https://netlib.org/utk/people/JackDongarra/etemplates/node96.html
    std::tuple<double, Eigen::VectorXd>
    inverse_power_method(const Eigen::SparseMatrix<double, Eigen::ColMajor> &eigen_matrix,
                         const Eigen::VectorXd &initial_vector_guess,
                         const double &initial_value_guess,
                         const double &tolerance,
                         const size_t &max_iteration)

    {
        Eigen::SparseMatrix<double, Eigen::ColMajor> shifted_matrix(eigen_matrix.cols(), eigen_matrix.rows()),
                                                     identity(eigen_matrix.cols(), eigen_matrix.rows());

        Eigen::VectorXd y(initial_vector_guess.size()),
                        v(initial_vector_guess.size()),
                        eigen_vector(initial_vector_guess.size());

        Eigen::BiCGSTAB<Eigen::SparseMatrix<double>> solver;

        double theta, quotient, eigen_value;

        y = initial_vector_guess / initial_vector_guess.squaredNorm();

        identity.setIdentity();

        for (size_t iteration=0; iteration<max_iteration; ++iteration)
        {
            shifted_matrix = eigen_matrix - initial_value_guess * identity;

            solver.compute(shifted_matrix);

            v = solver.solve(y);

            theta = v.transpose() * y;
            quotient = (y - theta * v).squaredNorm() / abs(theta);
            eigen_value = initial_value_guess + 1 / theta;

            eigen_vector = y / theta;

            if (quotient < tolerance)
                return std::make_tuple(eigen_value, eigen_vector);
        }
    }



    std::tuple<double, Eigen::VectorXd>
    inverse_power_method(const Eigen::SparseMatrix<double, Eigen::ColMajor> &eigen_matrix,
                         const double &initial_value_guess,
                         const double &tolerance,
                         const size_t &max_iteration)

    {
        Eigen::VectorXd initial_vector_guess = Eigen::VectorXd::Random(eigen_matrix.cols());

        return inverse_power_method(eigen_matrix, initial_vector_guess, initial_value_guess, tolerance, max_iteration);
    }



