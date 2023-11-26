#pragma once

#include "definitions.cpp"
#include "utils.cpp"
#include "numpy_interface.cpp"

struct FieldStructure
{
    public:
        Eigen::MatrixXd fields;
        Eigen::VectorXd itr_list, betas;
        double wavenumber, wavelength, dx_initial, dy_initial;

        FieldStructure(
            const Eigen::MatrixXd& fields,
            const Eigen::VectorXd &itr_list,
            const Eigen::VectorXd &betas,
            const double wavenumber,
            const double dx_initial,
            const double dy_initial
        )
        : fields(fields), itr_list(itr_list)
        {
            this->wavelength = 2. * PI / wavenumber;
        }

        double get_norm(const size_t &slice, const string &normalization_type)
        {
            if (normalization_type == "max")
                return this->get_norm_max(slice);
            if (normalization_type == "l2")
                return this->get_norm_l2(slice);
            if (normalization_type == "cmt")
                return this->get_norm_cmt(slice);
        }

        double get_norm_cmt(const size_t &slice)
        {
            double itr = itr_list[slice];

            return 0.5 * this->fields.col(slice).cwiseAbs2().sum() * dx_initial * itr * dy_initial * itr;
        }

        double get_norm_max(const size_t &slice)
        {
            return this->fields.col(slice).cwiseAbs().maxCoeff();
        }

        double get_norm_l2(const size_t &slice)
        {
            return sqrt(this->fields.col(slice).cwiseAbs2().sum());
        }

        Eigen::VectorXd get_norm_cmt_array()
        {
            Eigen::VectorXd term_0, term_1, norm_array;

            term_0 = 0.5 * this->fields.cwiseAbs2().colwise().sum();

            term_1 = itr_list.cwiseProduct(itr_list) * dx_initial * dy_initial;

            norm_array = term_0.cwiseProduct(term_1);

            return norm_array;
        }
};

struct SuperMode
{
    size_t mode_number;
    double wavenumber, wavelength;
    pybind11::array_t<double> mesh_gradient_py;
    Eigen::MatrixXd fields;
    Eigen::VectorXd index, betas, eigen_value, mesh_gradient;
    ModelParameters model_parameters;

    SuperMode(){}
    SuperMode(
        size_t mode_number,
        double wavenumber,
        pybind11::array_t<double> mesh_gradient_py,
        pybind11::array_t<double> fields_py,
        pybind11::array_t<double> index_py,
        pybind11::array_t<double> betas_py,
        pybind11::array_t<double> eigen_value_py,
        ModelParameters model_parameters
        )
        : mode_number(mode_number), wavenumber(wavenumber)
        {
            this->model_parameters = model_parameters;
            this->wavelength = 2. * PI / wavenumber;

            double *field_ptr = (double*) fields_py.request().ptr;
            Eigen::Map<Eigen::MatrixXd> mapping_fields(field_ptr, this->model_parameters.field_size, model_parameters.n_slice);
            this->fields = mapping_fields;

            double *index_ptr = (double*) index_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_index(index_ptr, model_parameters.n_slice, 1);
            this->index = mapping_index;

            double *beta_ptr = (double*) betas_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_betas(beta_ptr, model_parameters.n_slice, 1);
            this->betas = mapping_betas;

            double *eigen_value_ptr = (double*) eigen_value_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_eigen_value(eigen_value_ptr, model_parameters.n_slice, 1);
            this->eigen_value = mapping_eigen_value;

            double *mesh_gradient_ptr = (double*) mesh_gradient_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_mesh_gradient(mesh_gradient_ptr, this->model_parameters.field_size);
            this->mesh_gradient = mapping_mesh_gradient;
        }

    SuperMode(
        const size_t mode_number,
        const double wavenumber,
        const Eigen::VectorXd &mesh_gradient,
        const ModelParameters &model_parameters
        )
        : mode_number(mode_number),
          wavenumber(wavenumber),
          mesh_gradient(mesh_gradient)
        {
            this->model_parameters = model_parameters;
            this->wavelength = 2. * PI / wavenumber;
            this->fields = MatrixType(this->model_parameters.field_size, model_parameters.n_slice);
            this->eigen_value = VectorType(model_parameters.n_slice);
            this->betas = VectorType(model_parameters.n_slice);
            this->index = VectorType(model_parameters.n_slice);

        }

    pybind11::tuple get_state()
    {
        return pybind11::make_tuple(
            this->mode_number,
            this->wavenumber,
            this->get_mesh_gradient(),
            this->get_fields_py(),
            this->get_index_py(),
            this->get_betas_py(),
            this->get_eigen_value_py(),
            this->model_parameters
            );
    }

    double get_norm(const size_t &slice, const string &normalization_type) const
    {
        if (normalization_type == "max")
            return this->get_norm_max(slice);
        if (normalization_type == "l2")
            return this->get_norm_l2(slice);
        if (normalization_type == "cmt")
            return this->get_norm_cmt(slice);
    }

    Eigen::VectorXd get_norm_array(const string &normalization_type) const
    {
        if (normalization_type == "max")
            return this->get_norm_max_array();
        if (normalization_type == "l2")
            return this->get_norm_l2_array();
        if (normalization_type == "cmt")
            return this->get_norm_cmt_array();
        if (normalization_type == "scalar_coupling")
            return this->get_norm_scalar_coupling_array();
    }

    double get_norm_cmt(const size_t &slice) const
    {
        // Equation 7.35 from Bures
        double
            itr = model_parameters.itr_list[slice],
            factor = 0.5 * this->model_parameters.dx * itr * this->model_parameters.dy * itr;

        return this->fields.col(slice).cwiseAbs2().sum() * factor;
    }

    double get_norm_max(const size_t &slice) const
    {
        return this->fields.col(slice).cwiseAbs().maxCoeff();
    }

    double get_norm_l2(const size_t &slice) const
    {
        return sqrt(this->fields.col(slice).cwiseAbs2().sum());
    }

    Eigen::VectorXd get_norm_cmt_array() const
    {
        Eigen::VectorXd
            term_0 = 0.5 * this->fields.cwiseAbs2().colwise().sum(),
            term_1 = model_parameters.itr_list.cwiseProduct(model_parameters.itr_list) * this->model_parameters.dx * this->model_parameters.dy,
            norm_array = term_0.cwiseProduct(term_1);

        return norm_array;
    }

    Eigen::VectorXd get_norm_scalar_coupling_array() const
    {
        return this->fields.cwiseAbs2().colwise().sum() * (2 * PI);
    }

    Eigen::VectorXd get_norm_l2_array() const
    {
        double itr = 1;
        double
            dx = this->model_parameters.dx * itr,
            dy = this->model_parameters.dy * itr,
            dA = dx * dy;

        return this->fields.cwiseAbs2().colwise().sum() * dA;
    }

    Eigen::VectorXd get_norm_max_array() const
    {
        return this->fields.colwise().maxCoeff();
    }

    double get_overlap_integral(const SuperMode& other_supermode, const size_t &slice) const
    {
        return this->get_overlap_integral(other_supermode, slice, slice);
    }

    double get_overlap_integral(const SuperMode& other_supermode, const size_t &slice_0, const size_t &slice_1) const
    {
        return this->fields.col(slice_0).transpose() * other_supermode.fields.col(slice_0);
    }

    Eigen::MatrixXd get_overlap_integrals_with_mode(const SuperMode& other_supermode) const
    {
        Eigen::VectorXd
            overlap_integrals = this->fields.cwiseProduct(other_supermode.fields).colwise().sum().cwiseAbs();

        return overlap_integrals;
    }

    double get_trapez_integral(const Eigen::VectorXd &mesh, const double &dx, const double &dy) const
    {
        double
            dA = dx * dy,
            integral = mesh.sum() * dA;

        return integral;
    }

    Eigen::VectorXd get_gradient_overlap_with_mode(const SuperMode& other_supermode) const
    {
        Eigen::MatrixXd
            field_0 = this->fields,
            field_1 = other_supermode.fields,
            overlap = field_0.cwiseProduct(field_1).array().colwise() * mesh_gradient.array();

        Eigen::VectorXd
            dx = this->model_parameters.dx * model_parameters.itr_list,
            dy = this->model_parameters.dy * model_parameters.itr_list,
            dA = dx.cwiseProduct(dy),
            integral = overlap.colwise().sum(),
            gradient_overlap = integral.cwiseAbs() * dA;

        Eigen::VectorXd output(model_parameters.itr_list.size());

        double gradient_overlap_, dx_scaled, dy_scaled, itr;
        for (size_t slice = 0; slice < model_parameters.itr_list.size(); ++slice)
        {
            itr = model_parameters.itr_list[slice];
            dx_scaled = this->model_parameters.dx * itr;
            dy_scaled = this->model_parameters.dy * itr;
            gradient_overlap_ = get_trapez_integral(integral, dx_scaled, dy_scaled);
            output << gradient_overlap;
        }


        return gradient_overlap;
    }

    void normalize_fields()
    {
        for(size_t slice = 0; slice < this->fields.cols(); ++ slice)
            this->normalize_field_slice_l2(slice);

    }


    void arrange_fields()
    {
        for(size_t slice = 0; slice < this->fields.cols() - 1; ++ slice)
        {
            Eigen::VectorXd
                field_0 = this->fields.col(slice),
                field_1 = this->fields.col(slice + 1);


            double overlap = field_0.cwiseProduct(field_1).sum();
            if (overlap < 0)
                this->fields.col(slice + 1) *= -1;
        }
    }


    void normalize_field_slice_l2(const size_t slice)
    {
        double norm = this->get_field_slice_norm_custom(slice);

        this->fields.col(slice) /= sqrt(norm);
    }

    double get_field_slice_norm_custom(const size_t slice) const
    {
        Eigen::VectorXd field = this->fields.col(slice);

        double
            itr = this->model_parameters.itr_list[slice],
            dx = this->model_parameters.dx * itr,
            dy = this->model_parameters.dy * itr,
            dA = dx * dy,
            norm = field.cwiseAbs2().sum() * dA;

        return norm;
    }

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> get_normalized_coupling_with_mode(const SuperMode& other_supermode) const
    {
        Eigen::VectorXd
            integral = this->get_gradient_overlap_with_mode(other_supermode),
            beta_0 = this->betas,
            beta_1 = other_supermode.betas,
            delta_beta = (beta_0 - beta_1).cwiseInverse(),
            beta_prod = (beta_0.cwiseProduct(beta_1)).cwiseSqrt().cwiseInverse(),
            denominator = delta_beta.cwiseProduct(beta_prod);

            std::complex<double> scalar = - 0.5 * (std::complex<double>) J * pow(this->wavenumber, 2);

        return scalar * denominator.cwiseProduct(integral);
    }

    Eigen::VectorXd get_beating_length_with_mode(const SuperMode& other_supermode) const
    {
        Eigen::VectorXd
            beta_0 = this->betas,
            beta_1 = other_supermode.betas;

        return (beta_0 - beta_1).cwiseAbs().cwiseInverse() * (2 * PI);
    }

    Eigen::VectorXd get_adiabatic_with_mode(const SuperMode& other_supermode) const
    {
        Eigen::VectorXd delta_beta = this->betas - other_supermode.betas;

        Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> coupling = this->get_normalized_coupling_with_mode(other_supermode);

        return delta_beta.cwiseProduct(coupling.cwiseInverse()).cwiseAbs();
    }

    pybind11::array_t<double> get_overlap_integrals_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_overlap_integrals_with_mode(supermode),
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> get_gradient_overlap_with_mode_py(const SuperMode& supermode, const string &normalization_type="l2") const
    {
        return templated_eigen_to_ndarray(
            this->get_gradient_overlap_with_mode(supermode),
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<std::complex<double>> get_normalized_coupling_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_normalized_coupling_with_mode(supermode),
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> get_adiabatic_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_adiabatic_with_mode(supermode),
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> get_beating_length_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_beating_length_with_mode(supermode),
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> get_fields_py()
    {
        return templated_eigen_to_ndarray(
            this->fields,
            { model_parameters.n_slice, model_parameters.nx, model_parameters.ny }
        );
    }

    pybind11::array_t<double> get_index_py()
    {
        return templated_eigen_to_ndarray(
            this->index,
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> get_eigen_value_py()
    {
        return templated_eigen_to_ndarray(
            this->eigen_value,
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> get_betas_py()
    {
        return templated_eigen_to_ndarray(
            this->betas,
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> get_itr_list()
    {
        return templated_eigen_to_ndarray(
            this->model_parameters.itr_list,
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> get_mesh_gradient()
    {
        return templated_eigen_to_ndarray(
            this->mesh_gradient,
            { this->model_parameters.nx, this->model_parameters.ny }
        );
    }

    SuperMode get_supermode_from_tuple(pybind11::tuple tuple)
    {
        return SuperMode{
            tuple[0].cast<size_t>(),                             // mode_number
            tuple[1].cast<double>(),                             // k_initial
            tuple[2].cast<pybind11::array_t<double>>(),          // mesh_gradient
            tuple[3].cast<pybind11::array_t<double>>(),          // fields
            tuple[4].cast<pybind11::array_t<double>>(),          // index
            tuple[5].cast<pybind11::array_t<double>>(),          // betas
            tuple[6].cast<pybind11::array_t<double>>(),          // eigen_values
            tuple[7].cast<ModelParameters>()                     // model parameters
        }; // load
    }

};



