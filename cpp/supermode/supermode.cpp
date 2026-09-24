#include "supermode.h"

double SuperMode::get_norm(const size_t &slice, const std::string &normalization_type) const {
    if (normalization_type == "max")
        return this->get_norm_max(slice);
    if (normalization_type == "l2")
        return this->get_norm_l2(slice);
    if (normalization_type == "cmt")
        return this->get_norm_cmt(slice);

    return 0.0;
}


Eigen::VectorXd SuperMode::get_beating_length_with_mode(const SuperMode &other_supermode) const {

    Eigen::VectorXd beta_0 = this->betas, beta_1 = other_supermode.betas;

    return (beta_0 - beta_1).cwiseAbs().cwiseInverse() * (2 * PI);
}


Eigen::VectorXd SuperMode::get_adiabatic_with_mode(const SuperMode &other_supermode) const {
    Eigen::VectorXd delta_beta = (this->betas - other_supermode.betas).cwiseAbs();

    Eigen::Matrix<complex128, Eigen::Dynamic, 1>
        coupling = this->get_normalized_coupling_with_mode(other_supermode).cwiseAbs();

    return delta_beta.cwiseProduct(coupling.cwiseInverse()).cwiseAbs();
}

Eigen::Matrix<complex128, Eigen::Dynamic, 1> SuperMode::get_normalized_coupling_with_mode(const SuperMode &other_supermode) const {
    Eigen::VectorXd
        integral = this->get_gradient_field_overlap(other_supermode),
        beta_0 = this->betas,
        beta_1 = other_supermode.betas,
        delta_beta = (beta_0 - beta_1).cwiseInverse(),
        beta_prod = (beta_0.cwiseProduct(beta_1)).cwiseSqrt().cwiseInverse(),
        denominator = delta_beta.cwiseProduct(beta_prod);

    double k2 = pow(model_parameters.wavenumber, 2);

    std::complex<double> scalar = -complex128(0, 1) * k2 / 2.0;

    return scalar * denominator.cwiseProduct(integral);
}

Eigen::VectorXd SuperMode::get_gradient_field_overlap(const SuperMode &other_supermode) const {
    Eigen::VectorXd output(model_parameters.n_slice);

    for (size_t slice = 0; slice < model_parameters.n_slice; ++slice)
    {
        Eigen::Map<const Eigen::MatrixXd> field_0(this->fields.col(slice).data(), model_parameters.nx, model_parameters.ny);
        Eigen::Map<const Eigen::MatrixXd> field_1(other_supermode.fields.col(slice).data(), model_parameters.nx, model_parameters.ny);

        Eigen::MatrixXd overlap = model_parameters.mesh_gradient.cwiseProduct(field_0).cwiseProduct(field_1);

        double gradient_overlap = this->get_trapz_integral(
            overlap,
            model_parameters.dx_scaled[slice],
            model_parameters.dy_scaled[slice]
        );

        output[slice] = gradient_overlap;
    }

    return output;
}

double SuperMode::get_norm_cmt(const size_t &slice) const {
    // Equation 7.35 from Bures
    double
        itr = model_parameters.itr_list[slice],
        factor = 0.5 * this->model_parameters.dx * itr * this->model_parameters.dy * itr;

    return this->fields.col(slice).cwiseAbs2().sum() * factor;
}

double SuperMode::get_norm_max(const size_t &slice) const {
    return this->fields.col(slice).cwiseAbs().maxCoeff();
}

double SuperMode::get_norm_l2(const size_t &slice) const {
    return sqrt(this->fields.col(slice).cwiseAbs2().sum());
}

Eigen::VectorXd SuperMode::get_norm_cmt_array() const {
    Eigen::VectorXd
        term_0 = 0.5 * this->fields.cwiseAbs2().colwise().sum(),
        term_1 = model_parameters.itr_list.cwiseProduct(model_parameters.itr_list) * this->model_parameters.dx * this->model_parameters.dy,
        norm_array = term_0.cwiseProduct(term_1);

    return norm_array;
}

Eigen::VectorXd SuperMode::get_norm_scalar_coupling_array() const {
    return this->fields.cwiseAbs2().colwise().sum() * (2 * PI);
}

Eigen::VectorXd SuperMode::get_norm_max_array() const {
    return this->fields.colwise().maxCoeff();
}

double SuperMode::get_overlap_integral(const SuperMode& other_supermode, size_t slice) const {
    return this->get_overlap_integral(other_supermode, slice, slice);
}

double SuperMode::get_overlap_integral(const SuperMode& other_supermode, size_t slice_0, size_t slice_1) const {
    return this->fields.col(slice_0).transpose() * other_supermode.fields.col(slice_1);
}

Eigen::MatrixXd SuperMode::get_overlap_integrals_with_mode(const SuperMode& other_supermode) const {
    Eigen::VectorXd
        overlap_integrals = this->fields.cwiseProduct(other_supermode.fields).colwise().sum().cwiseAbs();

    return overlap_integrals;
}

void SuperMode::normalize_fields()
{
    for(int slice = 0; slice < this->fields.cols(); ++ slice) {
        double
            itr = this->model_parameters.itr_list[slice],
            dx = this->model_parameters.dx * itr,
            dy = this->model_parameters.dy * itr,
            norm = this->fields.col(slice).cwiseAbs2().sum() * dx * dy;

        this->fields.col(slice) /= sqrt(norm);
    }
}

void SuperMode::arrange_fields()
{
    for(int slice = 0; slice < this->fields.cols() - 1; ++ slice) {
        Eigen::VectorXd
            field_0 = this->fields.col(slice),
            field_1 = this->fields.col(slice + 1);


        double overlap = field_0.cwiseProduct(field_1).sum();
        if (overlap < 0)
            this->fields.col(slice + 1) *= -1;
    }
}

double SuperMode::get_trapz_integral(const Eigen::MatrixXd& mesh, double dx, double dy) const {
    // Precompute the scaling factors for the integral to simplify the loop calculations.
    const double dx_half = 0.5 * dx;
    const double dy_half = 0.5 * dy;

    // Initialize vector to store the integral of each row.
    Eigen::VectorXd row_integrals(mesh.rows());

    // Compute the integral for each row using the trapezoidal rule.
    for (int row = 0; row < mesh.rows(); ++row)
    {
        double row_integral = 0.0;
        for (int col = 0; col < mesh.cols() - 1; ++col)
            row_integral += (mesh(row, col) + mesh(row, col + 1)) * dx_half;

        row_integrals(row) = row_integral;
    }

    // Compute the total integral using the trapezoidal rule across the rows' integrals.
    double total_integral = 0.0;
    for (int row = 0; row < mesh.rows() - 1; ++row)
    {
        total_integral += (row_integrals(row) + row_integrals(row + 1)) * dy_half;
    }

    return total_integral;
}


bool SuperMode::is_same_symmetry(const SuperMode &other_supermode) const {
    // Return True if all boundaries condition are the same else False
    return (this->boundaries == other_supermode.boundaries);
}

Eigen::VectorXd SuperMode::get_field_interpolation(const double itr) const {
    // Create an output array to hold the interpolated field data
    Eigen::VectorXd interpolated_field(model_parameters.nx * model_parameters.ny);

    // Map the output array to an Eigen matrix for easier manipulation
    Eigen::Map<Eigen::MatrixXd> interpolated_field_map(
        interpolated_field.data(),
        model_parameters.nx,
        model_parameters.ny

    );

    // Perform linear interpolation for each point in the field
    for (size_t i = 0; i < model_parameters.nx; ++i) {
        for (size_t j = 0; j < model_parameters.ny; ++j) {
            // Find the two slices that bracket the desired iteration value
            size_t slice_0 = 0;
            while (slice_0 < model_parameters.n_slice - 1 && model_parameters.itr_list[slice_0 + 1] < itr) {
                ++slice_0;
            }
            size_t slice_1 = slice_0 + 1;

            // Get the field values at the two slices
            double field_0 = this->fields(i + j * model_parameters.nx, slice_0);
            double field_1 = this->fields(i + j * model_parameters.nx, slice_1);

            // Get the iteration values at the two slices
            double itr_0 = model_parameters.itr_list[slice_0];
            double itr_1 = model_parameters.itr_list[slice_1];

            // Perform linear interpolation
            if (itr_1 != itr_0) {
                interpolated_field_map(i, j) = field_0 + (field_1 - field_0) * (itr - itr_0) / (itr_1 - itr_0);
            } else {
                interpolated_field_map(i, j) = field_0; // or field_1, they are the same in this case
            }
        }
    }

    return interpolated_field;
}
