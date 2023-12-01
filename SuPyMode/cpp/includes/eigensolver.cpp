#pragma once

#include "supermode.cpp"
#include "extrapolate.cpp"
#include "progressbar.cpp"
#include "laplacian.cpp"
#include "definitions.cpp"
#include "utils.cpp"
#include "numpy_interface.cpp"
#include "eigen_solver.cpp"





class CppSolver : public BaseLaplacian
{
    public:
        size_t
            n_computed_mode,
            n_sorted_mode,
            max_iteration,
            iteration;

        double
            tolerance;

        std::vector<SuperMode>
            computed_supermodes,
            sorted_supermodes;

        Eigen::MatrixXd
            previous_eigen_vectors;

        Eigen::VectorXd
            previous_eigen_values;

        std::vector<double>
            alpha_vector;

        ModelParameters
            model_parameters;

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

   void generate_mode_set()
   {
     for (int mode_number=0; mode_number<n_computed_mode; ++mode_number)
     {
        SuperMode supermode = SuperMode(mode_number, this->model_parameters);
        computed_supermodes.push_back(supermode);
     }


     for (int mode_number=0; mode_number<n_sorted_mode; ++mode_number)
     {
        SuperMode supermode = SuperMode(mode_number, this->model_parameters);
        sorted_supermodes.push_back(supermode);
     }

   }

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd>
    compute_eigen_solution(const double &alpha)
    {
        this->compute_finit_diff_matrix();

        Spectra::SparseGenRealShiftSolve<double> op(eigen_matrix);

        Spectra::GenEigsRealShiftSolver<Spectra::SparseGenRealShiftSolve<double>> eigs(
            op,
            n_computed_mode,
            6 * n_computed_mode + 1, // See https://spectralib.org/doc/classspectra_1_1geneigsrealshiftsolver
            alpha
        );

        eigs.init();

        int nconv = eigs.compute(
            Spectra::SortRule::LargestMagn,
            this->max_iteration,
            this->tolerance,
            Spectra::SortRule::LargestMagn
        );

        eigen_matrix.resize(0, 0);

        Eigen::MatrixXd eigen_vectors = eigs.eigenvectors().real();
        Eigen::VectorXd eigen_values = eigs.eigenvalues().real();

        auto norm = eigen_vectors.col(0).cwiseProduct(eigen_vectors.col(1)).sum();

        this->sort_eigen_decomposition(eigen_vectors, eigen_values);

        return std::make_tuple(eigen_vectors, eigen_values);
    }

    void sort_eigen_decomposition(Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values)
    {
        if (this->iteration != 0)
            this->sort_eigen_with_fields(eigen_vectors, eigen_values);

        this->previous_eigen_values = eigen_values.block(0, 0, n_sorted_mode, 1);
        this->previous_eigen_vectors = eigen_vectors.block(0, 0, this->model_parameters.field_size, n_sorted_mode);
    }

    void sort_eigen_with_fields(Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values)
    {
        Eigen::Index max_index;

        Eigen::MatrixXd
            overlap_matrix = (previous_eigen_vectors.transpose() * eigen_vectors).cwiseAbs(),
            temporary_vector = eigen_vectors;

        Eigen::VectorXd temporary_value = eigen_values;

        for (size_t mode_0=0; mode_0<n_sorted_mode; ++mode_0)
        {
            double best_overlap = overlap_matrix.row(mode_0).maxCoeff(&max_index);

            if (best_overlap < 0.5)
                std::cout<<"Warning: bad mode field matching"<<"\n"
                         <<"\titeration: "<<this->iteration<<"\n"
                         <<"\tbest overlap: "<<best_overlap<<"\n"
                         <<"\tmodes index: "<<mode_0<<"\n"
                         <<"\tmax_index: "<<max_index<<"\n"
                         <<std::endl;

            eigen_vectors.col(mode_0) = temporary_vector.col(max_index);
            eigen_values.row(mode_0) = temporary_value.row(max_index);

            overlap_matrix.col(max_index) *= 0;
        }
    }

    void populate_sorted_supermodes(size_t slice, Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values)
    {
        for (SuperMode& mode : sorted_supermodes)
        {
            mode.eigen_value[slice] = eigen_values[mode.mode_number];

            mode.fields.col(slice) << eigen_vectors.col(mode.mode_number);


            mode.betas[slice] = sqrt( abs(mode.eigen_value[slice]) ) / this->model_parameters.itr_list[slice];
            mode.index[slice] = mode.betas[slice] / this->model_parameters.wavenumber;
        }
    }


    void loop_over_itr(size_t extrapolation_order = 1, double alpha = 0)
    {
        this->iteration = 0;

        Eigen::MatrixXd eigen_vectors;
        Eigen::VectorXd eigen_values;

        Extrapolator extrapolator = Extrapolator(this->model_parameters.ditr, extrapolation_order);

        ProgressBar progress_bar = ProgressBar(this->model_parameters.n_slice, 70, true, true);

        alpha_vector.reserve(this->model_parameters.n_slice);

        if (alpha == 0)
            alpha = -pow(this->model_parameters.wavenumber * this->compute_max_index(), 2);

        for (size_t slice=0; slice<this->model_parameters.n_slice; ++slice)
        {
            double itr_value = this->model_parameters.itr_list[slice];

            if (this->model_parameters.debug_mode > 0)
                progress_bar.show_next(this->model_parameters.itr_list[slice]);

            k_taper = this->model_parameters.wavenumber * this->model_parameters.itr_list[slice];



            std::tie(eigen_vectors, eigen_values) = this->compute_eigen_solution(alpha);

            if (this->model_parameters.debug_mode > 2)
            {
                std::cout.precision(16);
                std::cout<<" alpha guess: "<<alpha<<"\t"
                <<eigen_values.transpose()
                <<std::endl;
            }

            this->populate_sorted_supermodes(slice, eigen_vectors, eigen_values);

            alpha_vector.insert(alpha_vector.begin(), eigen_values[0]);

            alpha = extrapolator.extrapolate_next(alpha_vector);

            this->iteration++;
        }

        this->sort_mode_per_last_propagation_constant();

        this->normalize_mode_field();

        this->arrange_mode_field();
    }

    void sort_mode_per_last_propagation_constant()
    {
        if (this->model_parameters.debug_mode > 1)
            std::cout<<"Sorting supermode with propgation constant"<<std::endl;

        std::vector<double> list_index;
        for (SuperMode &supermode: sorted_supermodes)
        {
            double value = supermode.index.tail<1>().value();
            list_index.push_back(value);
        }

        std::vector<size_t> sorted_index = sort_indexes(list_index);

        inplace_reorder_vector(sorted_supermodes, sorted_index);

    }

    void arrange_mode_field()
    {
        if (this->model_parameters.debug_mode > 1)
            std::cout<<"Arranging fields"<<std::endl;

        for (SuperMode &supermode: sorted_supermodes)
            supermode.arrange_fields();
    }

    void normalize_mode_field()
    {
        if (this->model_parameters.debug_mode > 1)
            std::cout<<"Normalizing fields"<<std::endl;

        for (SuperMode &supermode: sorted_supermodes)
            supermode.normalize_fields();
    }



    double compute_max_index()
    {
        double max_index = 0.0;
        for (size_t i=0; i<size; ++i)
            if (meshPtr[i] > max_index)
                max_index = meshPtr[i];

        return max_index;
    }
};
