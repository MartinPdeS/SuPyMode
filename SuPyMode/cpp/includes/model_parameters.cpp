#pragma once

#include <definitions.cpp>



class ModelParameters
{
    public:

        double
            wavelength,
            wavenumber,
            dx,
            dy,
            ditr;

        size_t
            nx,
            ny,
            n_slice,
            field_size;

        pybind11::array_t<double>
            itr_list_py,
            mesh_gradient_py;

        Eigen::VectorXd
            mesh_gradient,
            itr_list,
            dx_scaled,
            dy_scaled,
            wavenumber_scaled;

        int
            debug_mode;

        ModelParameters(){}

        ModelParameters(
            double wavelength,
            pybind11::array_t<double> mesh_gradient_py,
            pybind11::array_t<double> itr_list_py,
            double dx,
            double dy
       ):
        wavelength(wavelength),
        itr_list_py(itr_list_py),
        mesh_gradient_py(mesh_gradient_py),
        dx(dx),
        dy(dy)
        {
            this->nx = mesh_gradient_py.request().shape[0];
            this->ny = mesh_gradient_py.request().shape[1];

            this->wavenumber = 2 * PI / this->wavelength;
            this->field_size = nx * ny;
            this->n_slice = itr_list_py.size();


            double *itr_list_ptr = (double*) itr_list_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_itr_list(itr_list_ptr, n_slice);
            this->itr_list = mapping_itr_list;

            double *mesh_gradient_ptr = (double*) mesh_gradient_py.request().ptr;
            Eigen::Map<Eigen::VectorXd> mapping_mesh_gradient(mesh_gradient_ptr, nx * ny);
            this->mesh_gradient = mapping_mesh_gradient;


            this->ditr = abs(itr_list[1] - itr_list[0]);

            this->compute_scaled_parameters();
        }

    void compute_scaled_parameters()
    {
        this->dx_scaled = itr_list * dx;
        this->dy_scaled = itr_list * dy;
        this->wavenumber_scaled = itr_list * this->wavenumber;
    }

    pybind11::tuple get_state()
    {
        return pybind11::make_tuple(
            this->wavelength,
            this->mesh_gradient_py,
            this->itr_list_py,
            this->dx,
            this->dy,
            this->n_slice
            );
    }

    std::ostream &operator<<(std::ostream &os)
    {
        return os << 'dx: '<< this->dx
                  << 'dy: '<< this->dy
                  << 'nx: '<< this->nx
                  << 'ny: '<< this->ny
                  << 'n_slice: '<< this->n_slice;
    }

};