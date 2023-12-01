#pragma once

#include "model_parameters.cpp"
#include "supermode.h"
#include "definitions.cpp"
#include "utils.cpp"
#include "numpy_interface.cpp"
#include "equations.cpp"


    SuperMode::SuperMode(){}
    SuperMode::SuperMode(
        size_t mode_number,
        pybind11::array_t<double> fields_py,
        pybind11::array_t<double> index_py,
        pybind11::array_t<double> betas_py,
        pybind11::array_t<double> eigen_value_py,
        ModelParameters model_parameters
        )
        : mode_number(mode_number)
        {
            this->model_parameters = model_parameters;

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
        }

    SuperMode::SuperMode(const size_t mode_number, const ModelParameters &model_parameters)
        : mode_number(mode_number)
        {
            this->model_parameters = model_parameters;
            this->fields = MatrixType(this->model_parameters.field_size, model_parameters.n_slice);
            this->eigen_value = VectorType(model_parameters.n_slice);
            this->betas = VectorType(model_parameters.n_slice);
            this->index = VectorType(model_parameters.n_slice);

        }

    pybind11::tuple SuperMode::get_state()
    {
        return pybind11::make_tuple(
            this->mode_number,
            this->get_fields_py(),
            this->get_index_py(),
            this->get_betas_py(),
            this->get_eigen_value_py(),
            this->model_parameters
            );
    }

    double SuperMode::get_norm(const size_t &slice, const std::string &normalization_type) const
    {
        if (normalization_type == "max")
            return this->get_norm_max(slice);
        if (normalization_type == "l2")
            return this->get_norm_l2(slice);
        if (normalization_type == "cmt")
            return this->get_norm_cmt(slice);
    }

    double SuperMode::get_norm_cmt(const size_t &slice) const
    {
        // Equation 7.35 from Bures
        double
            itr = model_parameters.itr_list[slice],
            factor = 0.5 * this->model_parameters.dx * itr * this->model_parameters.dy * itr;

        return this->fields.col(slice).cwiseAbs2().sum() * factor;
    }

    double SuperMode::get_norm_max(const size_t &slice) const
    {
        return this->fields.col(slice).cwiseAbs().maxCoeff();
    }

    double SuperMode::get_norm_l2(const size_t &slice) const
    {
        return sqrt(this->fields.col(slice).cwiseAbs2().sum());
    }

    Eigen::VectorXd SuperMode::get_norm_cmt_array() const
    {
        Eigen::VectorXd
            term_0 = 0.5 * this->fields.cwiseAbs2().colwise().sum(),
            term_1 = model_parameters.itr_list.cwiseProduct(model_parameters.itr_list) * this->model_parameters.dx * this->model_parameters.dy,
            norm_array = term_0.cwiseProduct(term_1);

        return norm_array;
    }

    Eigen::VectorXd SuperMode::get_norm_scalar_coupling_array() const
    {
        return this->fields.cwiseAbs2().colwise().sum() * (2 * PI);
    }

    Eigen::VectorXd SuperMode::get_norm_max_array() const
    {
        return this->fields.colwise().maxCoeff();
    }

    double SuperMode::get_overlap_integral(const SuperMode& other_supermode, const size_t &slice) const
    {
        return this->get_overlap_integral(other_supermode, slice, slice);
    }

    double SuperMode::get_overlap_integral(const SuperMode& other_supermode, const size_t &slice_0, const size_t &slice_1) const
    {
        return this->fields.col(slice_0).transpose() * other_supermode.fields.col(slice_0);
    }

    Eigen::MatrixXd SuperMode::get_overlap_integrals_with_mode(const SuperMode& other_supermode) const
    {
        Eigen::VectorXd
            overlap_integrals = this->fields.cwiseProduct(other_supermode.fields).colwise().sum().cwiseAbs();

        return overlap_integrals;
    }

    void SuperMode::normalize_fields()
    {
        for(size_t slice = 0; slice < this->fields.cols(); ++ slice)
            this->normalize_field_slice_l2(slice);

    }

    void SuperMode::arrange_fields()
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

    void SuperMode::normalize_field_slice_l2(const size_t slice)
    {
        double norm = this->get_field_slice_norm_custom(slice);

        this->fields.col(slice) /= sqrt(norm);
    }

    double SuperMode::get_field_slice_norm_custom(const size_t slice) const
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

    Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> SuperMode::get_normalized_coupling_with_mode(const SuperMode& other_supermode) const
    {
        return get_mode_normalized_coupling(*this, other_supermode, this->model_parameters);
    }

    Eigen::VectorXd SuperMode::get_beating_length_with_mode(const SuperMode& other_supermode) const
    {
        return get_mode_beating_length_with_mode(*this, other_supermode, this->model_parameters);
    }

    Eigen::VectorXd SuperMode::get_adiabatic_with_mode(const SuperMode& other_supermode) const
    {
        return get_mode_adiabatic_with_mode(*this, other_supermode, model_parameters);
    }

    pybind11::array_t<double> SuperMode::get_overlap_integrals_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_overlap_integrals_with_mode(supermode),
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<std::complex<double>> SuperMode::get_normalized_coupling_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_normalized_coupling_with_mode(supermode),
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> SuperMode::get_adiabatic_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_adiabatic_with_mode(supermode),
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> SuperMode::get_beating_length_with_mode_py(const SuperMode& supermode) const
    {
        return templated_eigen_to_ndarray(
            this->get_beating_length_with_mode(supermode),
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> SuperMode::get_fields_py()
    {
        return templated_eigen_to_ndarray(
            this->fields,
            { model_parameters.n_slice, model_parameters.nx, model_parameters.ny }
        );
    }

    pybind11::array_t<double> SuperMode::get_index_py()
    {
        return templated_eigen_to_ndarray(
            this->index,
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> SuperMode::get_eigen_value_py()
    {
        return templated_eigen_to_ndarray(
            this->eigen_value,
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> SuperMode::get_betas_py()
    {
        return templated_eigen_to_ndarray(
            this->betas,
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> SuperMode::get_itr_list()
    {
        return templated_eigen_to_ndarray(
            this->model_parameters.itr_list,
            { model_parameters.n_slice }
        );
    }

    pybind11::array_t<double> SuperMode::get_mesh_gradient()
    {
        return templated_eigen_to_ndarray(
            this->model_parameters.mesh_gradient,
            { this->model_parameters.nx, this->model_parameters.ny }
        );
    }

    SuperMode SuperMode::get_supermode_from_tuple(pybind11::tuple tuple)
    {
        return SuperMode{
            tuple[0].cast<size_t>(),                             // mode_number
            tuple[1].cast<pybind11::array_t<double>>(),          // fields
            tuple[2].cast<pybind11::array_t<double>>(),          // index
            tuple[3].cast<pybind11::array_t<double>>(),          // betas
            tuple[4].cast<pybind11::array_t<double>>(),          // eigen_values
            tuple[5].cast<ModelParameters>()                     // model parameters
        }; // load
    }



