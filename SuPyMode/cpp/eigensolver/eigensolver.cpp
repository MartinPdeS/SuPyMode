#include "eigensolver.h"


EigenSolver::EigenSolver(
    const ModelParameters &model_parameters,
    const pybind11::array_t<double> &finit_difference_triplets_py,
    const size_t n_computed_mode,
    const size_t n_sorted_mode,
    const size_t max_iteration,
    const double tolerance,
    const std::string &left_boundary,
    const std::string &right_boundary,
    const std::string &top_boundary,
    const std::string &bottom_boundary
)
    :
    model_parameters(model_parameters),
    n_computed_mode(n_computed_mode),
    n_sorted_mode(n_sorted_mode),
    max_iteration(max_iteration),
    tolerance(tolerance),
    left_boundary(left_boundary),
    right_boundary(right_boundary),
    top_boundary(top_boundary),
    bottom_boundary(bottom_boundary)
{
    this->finit_difference_triplets = numy_interface::convert_py_to_eigen<double>(
        finit_difference_triplets_py,
        finit_difference_triplets_py.request().shape[0],
        finit_difference_triplets_py.request().shape[1]
    );

    this->generate_mode_set();
}

void EigenSolver::compute_laplacian() {
    size_t number_of_triplet = finit_difference_triplets.cols();

    this->identity_matrix = Eigen::SparseMatrix<double, Eigen::ColMajor>(model_parameters.mesh.size(), model_parameters.mesh.size());
    this->laplacian_matrix = Eigen::SparseMatrix<double, Eigen::ColMajor>(model_parameters.mesh.size(), model_parameters.mesh.size());

    this->identity_matrix.setIdentity();

    this->laplacian_matrix.reserve(number_of_triplet);

    // Populate the sparse matrix directly
    for (size_t i = 0; i < number_of_triplet; ++i) {
        int row = static_cast<int>(finit_difference_triplets(0, i));
        int col = static_cast<int>(finit_difference_triplets(1, i));
        double value = finit_difference_triplets(2, i);
        this->laplacian_matrix.insert(col, row) = value;
    }

    this->eigen_matrix = -1.0 * laplacian_matrix;
}

void EigenSolver::compute_finit_diff_matrix() {
    for (long index = 0; index < model_parameters.mesh.size(); ++index)
        this->eigen_matrix.coeffRef(index, index) = - this->laplacian_matrix.coeffRef(index, index) - pow(model_parameters.mesh(index) * k_taper, 2);
}

void EigenSolver::generate_mode_set()
{
    for (size_t mode_number=0; mode_number < n_computed_mode; ++mode_number)
    {
        SuperMode supermode = SuperMode(mode_number, this->model_parameters, this->left_boundary, this->right_boundary, this->top_boundary, this->bottom_boundary);
        computed_supermodes.push_back(supermode);
    }

    for (size_t mode_number=0; mode_number<n_sorted_mode; ++mode_number)
    {
        SuperMode supermode = SuperMode(mode_number, this->model_parameters, this->left_boundary, this->right_boundary, this->top_boundary, this->bottom_boundary);
        sorted_supermodes.push_back(supermode);
    }

}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd> EigenSolver::compute_eigen_solution(const double alpha)
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

    eigs.compute(
        Spectra::SortRule::LargestMagn,
        this->max_iteration,
        this->tolerance,
        Spectra::SortRule::LargestMagn
    );

    Eigen::MatrixXd eigen_vectors = eigs.eigenvectors().real();
    Eigen::VectorXd eigen_values = eigs.eigenvalues().real();

    this->sort_eigen_decomposition(eigen_vectors, eigen_values);

    return std::make_tuple(eigen_vectors, eigen_values);
}

void EigenSolver::sort_eigen_decomposition(Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values)
{
    if (this->iteration != 0)
        this->sort_eigen_with_fields(eigen_vectors, eigen_values);

    this->previous_eigen_values = eigen_values.block(0, 0, n_sorted_mode, 1);
    this->previous_eigen_vectors = eigen_vectors.block(0, 0, this->model_parameters.nx * this->model_parameters.ny, n_sorted_mode);
}

void EigenSolver::sort_eigen_with_fields(Eigen::MatrixXd& eigen_vectors, Eigen::VectorXd& eigen_values)
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
            printf("Warning: bad mode field matching\n\titeration: %zu\n\tbest overlap: %f\n\tmodes index: %zu\n\tmax_index: %zu\n", this->iteration, best_overlap, mode_0, max_index);


        eigen_vectors.col(mode_0) = temporary_vector.col(max_index);
        eigen_values.row(mode_0) = temporary_value.row(max_index);

        overlap_matrix.col(max_index) *= 0;
    }
}

void EigenSolver::populate_sorted_supermodes(size_t slice, Eigen::MatrixXd &eigen_vectors, Eigen::VectorXd &eigen_values)
{
    for (SuperMode& mode : sorted_supermodes)
    {
        mode.eigen_value[slice] = eigen_values[mode.mode_number];

        mode.fields.col(slice) << eigen_vectors.col(mode.mode_number);

        mode.betas[slice] = sqrt( abs(mode.eigen_value[slice]) ) / this->model_parameters.itr_list[slice];
        mode.index[slice] = mode.betas[slice] / this->model_parameters.wavenumber;
    }
}

void EigenSolver::loop_over_itr(size_t extrapolation_order, double alpha)
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

void EigenSolver::sort_mode_per_last_propagation_constant()
{
    if (this->model_parameters.debug_mode > 1)
        printf("Sorting supermode with propgation constant\n");

    std::vector<double> list_index;
    for (SuperMode &supermode: sorted_supermodes)
    {
        double value = supermode.index.tail<1>().value();
        list_index.push_back(value);
    }

    std::vector<size_t> sorted_index = sort_indexes(list_index);

    inplace_reorder_vector(sorted_supermodes, sorted_index);

}

void EigenSolver::arrange_mode_field()
{
    if (this->model_parameters.debug_mode > 1)
        printf("Arranging fields\n");

    for (SuperMode &supermode: sorted_supermodes)
        supermode.arrange_fields();
}

void EigenSolver::normalize_mode_field()
{
    if (this->model_parameters.debug_mode > 1)
        printf("Normalizing fields\n");

    for (SuperMode &supermode: sorted_supermodes)
        supermode.normalize_fields();
}


double EigenSolver::compute_max_index()
{
    double max_index = 0.0;
    for (long i=0; i<model_parameters.mesh.size(); ++i)
        if (model_parameters.mesh(i) > max_index)
            max_index = model_parameters.mesh(i);

    return max_index;
}
