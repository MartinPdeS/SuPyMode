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

            return 0.5 * this->fields.col(slice).cwiseAbs2().sum() * this->dx_initial * itr * this->dy_initial * itr;
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

            term_1 = itr_list.cwiseProduct(itr_list) * this->dx_initial * this->dy_initial;

            norm_array = term_0.cwiseProduct(term_1);

            return norm_array;
        }
};

struct SuperMode
{
    size_t mode_number;
    double wavenumber, wavelength, dx_initial, dy_initial;
    pybind11::array_t<double> itr_list_py, mesh_gradient_py;
    Eigen::MatrixXd fields;
    Eigen::VectorXd index, betas, eigen_value, itr_list, mesh_gradient;
    size_t nx, ny, n_slice;

    SuperMode(){}
    SuperMode(
        size_t mode_number,
        double wavenumber,
        double dx_initial,
        double dy_initial,
        pybind11::array_t<double> itr_list_py,
        pybind11::array_t<double> mesh_gradient_py,
        pybind11::array_t<double> fields_py,
        pybind11::array_t<double> index_py,
        pybind11::array_t<double> betas_py,
        pybind11::array_t<double> eigen_value_py
        )
        : mode_number(mode_number), wavenumber(wavenumber), dx_initial(dx_initial), dy_initial(dy_initial)
        {
            this->wavelength = 2. * PI / wavenumber;
            this->nx = mesh_gradient_py.request().shape[0];
            this->ny = mesh_gradient_py.request().shape[1];
            this->n_slice = itr_list_py.request().shape[0];

            double *field_ptr = (double*) fields_py.request().ptr;
            Eigen::Map<Eigen::MatrixXd> mapping_fields(field_ptr, nx * ny, n_slice);
            this->fields = mapping_fields;

            double *index_ptr = (double*) index_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_index(index_ptr, n_slice, 1);
            this->index = mapping_index;

            double *beta_ptr = (double*) betas_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_betas(beta_ptr, n_slice, 1);
            this->betas = mapping_betas;

            double *eigen_value_ptr = (double*) eigen_value_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_eigen_value(eigen_value_ptr, n_slice, 1);
            this->eigen_value = mapping_eigen_value;

            double *itr_list_ptr = (double*) itr_list_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_itr_list(itr_list_ptr, n_slice);
            this->itr_list = mapping_itr_list;

            double *mesh_gradient_ptr = (double*) mesh_gradient_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_mesh_gradient(mesh_gradient_ptr, nx*ny);
            this->mesh_gradient = mapping_mesh_gradient;
        }

    SuperMode(
        const size_t mode_number,
        const double wavenumber,
        const double dx_initial,
        const double dy_initial,
        const Eigen::VectorXd &mesh_gradient,
        const Eigen::VectorXd &itr_list,
        const size_t nx,
        const size_t ny
        )
        : mode_number(mode_number), wavenumber(wavenumber), dx_initial(dx_initial), dy_initial(dy_initial), itr_list(itr_list), mesh_gradient(mesh_gradient)
        {
            this->wavelength = 2. * PI / wavenumber;
            this->nx = nx;
            this->ny = ny;
            this->n_slice = itr_list.size();
            this->fields = MatrixType(nx * ny, n_slice);
            this->eigen_value = VectorType(n_slice);
            this->betas = VectorType(n_slice);
            this->index = VectorType(n_slice);

        }

    pybind11::tuple get_state()
    {
        return pybind11::make_tuple(
            this->mode_number,
            this->wavenumber,
            this->dx_initial,
            this->dy_initial,
            this->get_itr_list(),
            this->get_mesh_gradient(),
            this->get_fields_py(),
            this->get_index_py(),
            this->get_betas_py(),
            this->get_eigen_value_py()
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
            itr = itr_list[slice],
            factor = 0.5 * this->dx_initial * itr * this->dy_initial * itr;

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
            term_1 = itr_list.cwiseProduct(itr_list) * this->dx_initial * this->dy_initial, 
            norm_array = term_0.cwiseProduct(term_1);

        return norm_array;
    }

    Eigen::VectorXd get_norm_scalar_coupling_array() const
    {
        return this->fields.cwiseAbs2().colwise().sum() * (2 * PI);
    }

    Eigen::VectorXd get_norm_l2_array() const
    {
        return this->fields.cwiseAbs2().colwise().sum();
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
        Eigen::VectorXd 
            mode_field_0 = this->get_normalized_field(slice_0, "l2"),
            mode_field_1 = other_supermode.get_normalized_field(slice_1, "l2");

        return mode_field_0.transpose() * mode_field_1;
    }

    Eigen::MatrixXd get_overlap_integrals_with_mode(const SuperMode& other_supermode) const
    {
        Eigen::MatrixXd  
            normalized_field_0 = this->get_normalized_fields("l2"),
            normalized_field_1 = other_supermode.get_normalized_fields("l2");

        Eigen::VectorXd     
            overlap_integrals = normalized_field_0.cwiseProduct(normalized_field_1).colwise().sum().cwiseAbs();

        return overlap_integrals;
    }

    Eigen::VectorXd get_normalized_field(const size_t slice, const string &normalization_type) const
    {
        double norm_field = this->get_norm(slice, normalization_type);

        Eigen::VectorXd field = this->fields.col(slice);

        return field / sqrt(norm_field);
    }

    Eigen::MatrixXd get_normalized_fields(const string &normalization_type) const
    {
        Eigen::VectorXd norm_array = this->get_norm_array(normalization_type);

        Eigen::MatrixXd normalized_fields = this->fields.array().rowwise() / norm_array.cwiseSqrt().transpose().array();

        return normalized_fields;
    }

    Eigen::VectorXd get_gradient_overlap_with_mode(const SuperMode& other_supermode, const string &normalization_type) const
    {
        Eigen::MatrixXd 
            field_0 = this->get_normalized_fields(normalization_type),
            field_1 = other_supermode.get_normalized_fields(normalization_type),                         
            overlap = field_0.cwiseProduct(field_1).array().colwise() * mesh_gradient.array();

        return overlap.colwise().sum().cwiseAbs();
    }

    void normalize_fields()
    {
        Eigen::VectorXd norm_array = this->get_norm_array("l2");

        this->fields = this->fields.array().rowwise() / norm_array.transpose().array().cwiseSqrt();
    }

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> get_normalized_coupling_with_mode(const SuperMode& other_supermode) const
    {
        Eigen::VectorXd 
            integral = this->get_gradient_overlap_with_mode(other_supermode, "scalar_coupling"),
            beta_0 = this->betas,
            beta_1 = other_supermode.betas,
            delta_beta = (beta_0 - beta_1),
            term0 = delta_beta.cwiseInverse(),
            term1 = (beta_0.cwiseProduct(beta_1)).cwiseSqrt().cwiseInverse();

        std::complex<double> scalar = - (std::complex<double>) 0.5 * J * pow(this->wavenumber, 2);

        return scalar * term0.cwiseProduct(term1).cwiseProduct(integral);
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
            { n_slice }
        );
    }

    pybind11::array_t<double> get_gradient_overlap_with_mode_py(const SuperMode& supermode, const string &normalization_type="l2") const
    {
        return templated_eigen_to_ndarray(
            this->get_gradient_overlap_with_mode(supermode, normalization_type),
            { n_slice }
        );
    }

    pybind11::array_t<std::complex<double>> get_normalized_coupling_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_normalized_coupling_with_mode(supermode),
            { n_slice }
        );
    }

    pybind11::array_t<double> get_adiabatic_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_adiabatic_with_mode(supermode),
            { n_slice }
        );
    }

    pybind11::array_t<double> get_beating_length_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_beating_length_with_mode(supermode),
            { n_slice }
        );
    }

    pybind11::array_t<double> get_fields_py()
    {
        return templated_eigen_to_ndarray(
            this->fields,
            { n_slice, nx, ny }
        );
    }

    pybind11::array_t<double> get_index_py()
    {
        return templated_eigen_to_ndarray(
            this->index,
            { n_slice }
        );
    }

    pybind11::array_t<double> get_eigen_value_py()
    {
        return templated_eigen_to_ndarray(
            this->eigen_value,
            { n_slice }
        );
    }

    pybind11::array_t<double> get_betas_py()
    {
        return templated_eigen_to_ndarray(
            this->betas,
            { n_slice }
        );
    }

    pybind11::array_t<double> get_itr_list()
    {
        return templated_eigen_to_ndarray(
            this->itr_list,
            { n_slice }
        );
    }

    pybind11::array_t<double> get_mesh_gradient()
    {
        return templated_eigen_to_ndarray(
            this->mesh_gradient,
            { nx, ny }
        );
    }

    SuperMode get_supermode_from_tuple(pybind11::tuple tuple)
    {
        return SuperMode{
            tuple[0].cast<size_t>(),                             // mode_number
            tuple[1].cast<double>(),                             // k_initial
            tuple[2].cast<double>(),                             // dx
            tuple[3].cast<double>(),                             // dy
            tuple[4].cast<pybind11::array_t<double>>(),          // itr_list
            tuple[5].cast<pybind11::array_t<double>>(),          // mesh_gradient
            tuple[6].cast<pybind11::array_t<double>>(),          // fields
            tuple[7].cast<pybind11::array_t<double>>(),          // index
            tuple[8].cast<pybind11::array_t<double>>(),          // betas
            tuple[9].cast<pybind11::array_t<double>>()           // eigen_values
        }; // load
    }

};



