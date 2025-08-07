#include "model_parameters.h"


ModelParameters::ModelParameters(
    const pybind11::array_t<double> &mesh_py,
    const pybind11::array_t<double> &x_vector_py,
    const pybind11::array_t<double> &y_vector_py,
    const pybind11::array_t<double> &itr_list_py,
    const double wavelength,
    const double dx,
    const double dy,
    const int debug_mode)
    :
        wavelength(wavelength),
        dx(dx),
        dy(dy),
        mesh_py(mesh_py),
        x_vector_py(x_vector_py),
        y_vector_py(y_vector_py),
        itr_list_py(itr_list_py),
        debug_mode(debug_mode)
{
    init();
}

void ModelParameters::compute_scaled_parameters() {
    dx_scaled = itr_list * dx;
    dy_scaled = itr_list * dy;
    wavenumber_scaled = itr_list * wavenumber;
}

pybind11::tuple
ModelParameters::get_pickle(const ModelParameters &model_parameters) {
    return pybind11::make_tuple(
        model_parameters.mesh_py,
        model_parameters.x_vector_py,
        model_parameters.y_vector_py,
        model_parameters.itr_list_py,
        model_parameters.wavelength,
        model_parameters.dx,
        model_parameters.dy,
        model_parameters.debug_mode
    );
}

ModelParameters
ModelParameters::build_from_tuple(const pybind11::tuple &tuple) {
    return ModelParameters{
        tuple[0].cast<pybind11::array_t<double>>(),
        tuple[1].cast<pybind11::array_t<double>>(),
        tuple[2].cast<pybind11::array_t<double>>(),
        tuple[3].cast<pybind11::array_t<double>>(),
        tuple[4].cast<double>(),
        tuple[5].cast<double>(),
        tuple[6].cast<double>(),
        tuple[7].cast<int>()
    };
}

// Definition of the friend function
std::ostream& operator<<(std::ostream& os, const ModelParameters& params) {
    os << "dx: " << params.dx
       << " dy: " << params.dy
       << " nx: " << params.nx
       << " ny: " << params.ny
       << " n_slice: " << params.n_slice;
    return os;
}


void ModelParameters::init() {
    nx = mesh_py.request().shape[1];
    ny = mesh_py.request().shape[0];
    n_slice = itr_list_py.request().size;
    wavenumber = 2 * PI / wavelength;

    mesh = NumpyInterface::convert_py_to_eigen<double>(mesh_py, ny, nx);
    x_vector = NumpyInterface::convert_py_to_eigen<double>(x_vector_py, nx);
    y_vector = NumpyInterface::convert_py_to_eigen<double>(y_vector_py, ny);
    itr_list = NumpyInterface::convert_py_to_eigen<double>(itr_list_py, n_slice);
    mesh_gradient = get_rho_gradient_time_rho(mesh.cwiseProduct(mesh), y_vector, x_vector);
    mesh_gradient_py = NumpyInterface::eigen_to_ndarray<double>(mesh_gradient, {ny, nx});

    this->compute_scaled_parameters();
    ditr = std::abs(itr_list[1] - itr_list[0]);
}



// -
