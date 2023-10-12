#pragma once

    #include <eigen/Eigen>

    class RitzIterator  //https://en.wikipedia.org/wiki/Rayleigh%E2%80%93Ritz_method
    {
        Eigen::SparseMatrix<double, Eigen::ColMajor> A;
        size_t n_ritz;

        RitzIterator(const size_t &n_ritz): n_ritz(n_ritz)
        {

        }

        void set_matrix(const Eigen::SparseMatrix<double, Eigen::ColMajor> &A)
        {
            this->A = A;
        }

        void solve()
        {
            Eigen::MatrixXd V = Eigen::MatrixXd::Random(A.cols(), n_ritz);

            Eigen::HouseholderQR<Eigen::MatrixXd> qr(V);

            V = qr.householderQ();

            // Eigen::MatrixXd M = V.transpose() * A * V


        }



    };



