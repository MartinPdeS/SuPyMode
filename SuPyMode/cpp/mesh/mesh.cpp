#include "mesh.h"


std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
compute_gradient_2p(const Eigen::MatrixXd& image, double dx, double dy) {
    int y_size = image.rows();
    int x_size = image.cols();

    // Initialize matrices to store the gradients in x and y directions
    Eigen::MatrixXd gradient_x(y_size, x_size);
    Eigen::MatrixXd gradient_y(y_size, x_size);

    for (int y_index = 0; y_index < y_size; ++y_index) {
        for (int x_index = 0; x_index < x_size; ++x_index) {

            // Compute dxVal
            double dxVal = 0.0;
            if (x_index> 0 && x_index + 1 < x_size) {
                dxVal = (image(y_index, x_index + 1) - image(y_index, x_index - 1)) / (2 * dx);
            } else if (x_index == 0 && x_index + 2 < x_size) {  // Use forward difference at the left boundary
                dxVal = (-3 * image(y_index, x_index) + 4 * image(y_index, x_index + 1) - image(y_index, x_index + 2)) / (2 * dx);
            } else if (x_index == x_size - 1 && x_index - 2 >= 0) {  // Use backward difference at the right boundary
                dxVal = (3 * image(y_index, x_index) - 4 * image(y_index, x_index - 1) + image(y_index, x_index - 2)) / (2 * dx);
            } else if (x_index == 0 && x_index + 1 < x_size) {  // Forward difference when there is only one neighboring point
                dxVal = (image(y_index, x_index + 1) - image(y_index, x_index)) / dx;
            } else if (x_index == x_size - 1 && x_index - 1 >= 0) {  // Backward difference when there is only one neighboring point
                dxVal = (image(y_index, x_index) - image(y_index, x_index - 1)) / dx;
            }

            // Compute dyVal
            double dyVal = 0.0;
            if (y_index> 0 && y_index + 1 < y_size) {
                dyVal = (image(y_index + 1, x_index) - image(y_index - 1, x_index)) / (2 * dy);
            } else if (y_index == 0 && y_index + 2 < y_size) {  // Use forward difference at the top boundary
                dyVal = (-3 * image(y_index, x_index) + 4 * image(y_index + 1, x_index) - image(y_index + 2, x_index)) / (2 * dy);
            } else if (y_index == y_size - 1 && y_index - 2 >= 0) {  // Use backward difference at the bottom boundary
                dyVal = (3 * image(y_index, x_index) - 4 * image(y_index - 1, x_index) + image(y_index - 2, x_index)) / (2 * dy);
            } else if (y_index == 0 && y_index + 1 < y_size) {  // Forward difference when there is only one neighboring point
                dyVal = (image(y_index + 1, x_index) - image(y_index, x_index)) / dy;
            } else if (y_index == y_size - 1 && y_index - 1 >= 0) {  // Backward difference when there is only one neighboring point
                dyVal = (image(y_index, x_index) - image(y_index - 1, x_index)) / dy;
            }

            gradient_x(y_index, x_index) = dxVal;
            gradient_y(y_index, x_index) = dyVal;
        }
    }

    return {gradient_x, gradient_y};
}


// Rows is y-axis -- Cols is x--axis
std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
compute_gradient_5p(const Eigen::MatrixXd& image, const double dx, const double dy) {
    size_t y_size = image.rows();
    size_t x_size = image.cols();

    // Initialize matrices to store the gradients
    Eigen::MatrixXd gradient_x(y_size, x_size);
    Eigen::MatrixXd gradient_y(y_size, x_size);

    for (size_t y_index = 0; y_index < y_size; ++y_index) {
        for (size_t x_index = 0; x_index < x_size; ++x_index) {
            double dxVal = 0.0;
            double dyVal = 0.0;

            // Calculate dxVal
            if (x_index >= 2 && x_index + 2 < x_size) {
                dxVal = (image(y_index, x_index - 2) - 8 * image(y_index, x_index - 1) + 8 * image(y_index, x_index + 1) - image(y_index, x_index + 2)) / (12 * dx);
            } else if (x_index > 0 && x_index + 1 < x_size) {
                dxVal = (image(y_index, x_index + 1) - image(y_index, x_index - 1)) / (2 * dx);
            } else if (x_index == 0 && x_index + 1 < x_size) {
                dxVal = (image(y_index, x_index + 1) - image(y_index, x_index)) / dx;
            } else if (x_index == x_size - 1 && x_index - 1 >= 0) {
                dxVal = (image(y_index, x_index) - image(y_index, x_index - 1)) / dx;
            }

            // Calculate dyVal
            if (y_index >= 2 && y_index + 2 < y_size) {
                dyVal = (image(y_index - 2, x_index) - 8 * image(y_index - 1, x_index) + 8 * image(y_index + 1, x_index) - image(y_index + 2, x_index)) / (12 * dy);
            } else if (y_index > 0 && y_index + 1 < y_size) {
                dyVal = (image(y_index + 1, x_index) - image(y_index - 1, x_index)) / (2 * dy);
            } else if (y_index == 0 && y_index + 1 < y_size) {
                dyVal = (image(y_index + 1, x_index) - image(y_index, x_index)) / dy;
            } else if (y_index == y_size - 1 && y_index - 1 >= 0) {
                dyVal = (image(y_index, x_index) - image(y_index - 1, x_index)) / dy;
            }

            gradient_x(y_index, x_index) = dxVal;
            gradient_y(y_index, x_index) = dyVal;
        }
    }

    return {gradient_x, gradient_y};
}


std::pair<Eigen::MatrixXd, Eigen::MatrixXd>
compute_gradient_7p(const Eigen::MatrixXd& image, double dx, double dy) {
    size_t y_size = image.rows();
    size_t x_size = image.cols();

    // Initialize matrices to store the gradients
    Eigen::MatrixXd gradient_x(y_size, x_size);
    Eigen::MatrixXd gradient_y(y_size, x_size);

    for (size_t y_index = 0; y_index < y_size; ++y_index) {
        for (size_t x_index = 0; x_index < x_size; ++x_index) {
            double dxVal = 0.0;
            if (x_index>= 3 && x_index + 3 < x_size) {
                dxVal = (-image(y_index, x_index - 3) + 9 * image(y_index, x_index - 2) - 45 * image(y_index, x_index - 1) + 45 * image(y_index, x_index + 1) - 9 * image(y_index, x_index + 2) + image(y_index, x_index + 3)) / (60 * dx);
            } else if (x_index> 0 && x_index + 1 < x_size) {
                dxVal = (image(y_index, x_index + 1) - image(y_index, x_index - 1)) / (2 * dx);  // Central difference
            } else if (x_index == 0 && x_index + 2 < x_size) {  // Forward difference at the left boundary
                dxVal = (-3 * image(y_index, x_index) + 4 * image(y_index, x_index + 1) - image(y_index, x_index + 2)) / (2 * dx);
            } else if (x_index == x_size - 1 && x_index - 2 >= 0) {  // Backward difference at the right boundary
                dxVal = (3 * image(y_index, x_index) - 4 * image(y_index, x_index - 1) + image(y_index, x_index - 2)) / (2 * dx);
            } else if (x_index == 0 && x_index + 1 < x_size) {  // Forward difference with one neighboring point
                dxVal = (image(y_index, x_index + 1) - image(y_index, x_index)) / dx;
            } else if (x_index == x_size - 1 && x_index - 1 >= 0) {  // Backward difference with one neighboring point
                dxVal = (image(y_index, x_index) - image(y_index, x_index - 1)) / dx;
            }

            double dyVal = 0.0;
            if (y_index>= 3 && y_index + 3 < y_size) {
                dyVal = (-image(y_index - 3, x_index) + 9 * image(y_index - 2, x_index) - 45 * image(y_index - 1, x_index) + 45 * image(y_index + 1, x_index) - 9 * image(y_index + 2, x_index) + image(y_index + 3, x_index)) / (60 * dy);
            } else if (y_index > 0 && y_index + 1 < y_size) {
                dyVal = (image(y_index + 1, x_index) - image(y_index - 1, x_index)) / (2 * dy);  // Central difference
            } else if (y_index == 0 && y_index + 2 < y_size) {  // Forward difference at the top boundary
                dyVal = (-3 * image(y_index, x_index) + 4 * image(y_index + 1, x_index) - image(y_index + 2, x_index)) / (2 * dy);
            } else if (y_index == y_size - 1 && y_index - 2 >= 0) {  // Backward difference at the bottom boundary
                dyVal = (3 * image(y_index, x_index) - 4 * image(y_index - 1, x_index) + image(y_index - 2, x_index)) / (2 * dy);
            } else if (y_index == 0 && y_index + 1 < y_size) {  // Forward difference with one neighboring point
                dyVal = (image(y_index + 1, x_index) - image(y_index, x_index)) / dy;
            } else if (y_index == y_size - 1 && y_index - 1 >= 0) {  // Backward difference with one neighboring point
                dyVal = (image(y_index, x_index) - image(y_index - 1, x_index)) / dy;
            }

            gradient_x(y_index, x_index) = dxVal;
            gradient_y(y_index, x_index) = dyVal;
        }
    }

    return {gradient_x, gradient_y};
}


Eigen::MatrixXd
get_rho_gradient_time_rho(const Eigen::MatrixXd &mesh, const Eigen::VectorXd &y_vector, const Eigen::VectorXd &x_vector) {
    size_t x_size = x_vector.size();
    size_t y_size = y_vector.size();

    double dx = x_vector(1) - x_vector(0);
    double dy = y_vector(1) - y_vector(0);

    auto gradients = compute_gradient_5p(mesh, dy, dx);
    Eigen::MatrixXd gradient_x = gradients.first;
    Eigen::MatrixXd gradient_y = gradients.second;

    Eigen::MatrixXd rho_gradient(y_size, x_size);

    double angle, cos_angle, sin_angle;

    for (size_t x_index = 0; x_index < x_size; ++x_index) {
        for (size_t y_index = 0; y_index < y_size; ++y_index) {
            double x = x_vector(x_index);
            double y = y_vector(y_index);
            double rho = std::sqrt(pow(x, 2) + pow(y, 2));

            angle = std::atan2(y, x);
            cos_angle = std::cos(angle);
            sin_angle = std::sin(angle);

            rho_gradient(y_index, x_index) = (gradient_x(y_index, x_index) * cos_angle + gradient_y(y_index, x_index) * sin_angle) * rho;
       }
    }

    return rho_gradient;
}


Eigen::MatrixXd
get_rho_gradient(const Eigen::MatrixXd &mesh, const Eigen::VectorXd &y_vector, const Eigen::VectorXd &x_vector) {
    size_t x_size = x_vector.size();
    size_t y_size = y_vector.size();

    double dx = x_vector(1) - x_vector(0);
    double dy = y_vector(1) - y_vector(0);

    auto gradients = compute_gradient_5p(mesh, dy, dx);
    Eigen::MatrixXd gradient_x = gradients.first;
    Eigen::MatrixXd gradient_y = gradients.second;

    Eigen::MatrixXd rho_gradient(y_size, x_size);

    double angle, cos_angle, sin_angle;

    for (size_t x_index = 0; x_index < x_size; ++x_index) {
        for (size_t y_index = 0; y_index < y_size; ++y_index) {
            double x = x_vector(x_index);
            double y = y_vector(y_index);

            angle = std::atan2(y, x);
            cos_angle = std::cos(angle);
            sin_angle = std::sin(angle);
            rho_gradient(y_index, x_index) = (gradient_x(y_index, x_index) * cos_angle + gradient_y(y_index, x_index) * sin_angle);
       }
    }

    return rho_gradient;
}


pybind11::array_t<double>
get_rho_gradient_py(const pybind11::array_t<double> &mesh_py, const pybind11::array_t<double> &x_vector_py, const pybind11::array_t<double> &y_vector_py) {
    size_t x_size = x_vector_py.request().size;
    size_t y_size = y_vector_py.request().size;

    Eigen::VectorXd x_vector = numy_interface::convert_py_to_eigen(x_vector_py, x_size);
    Eigen::VectorXd y_vector = numy_interface::convert_py_to_eigen(y_vector_py, y_size);
    Eigen::MatrixXd mesh = numy_interface::convert_py_to_eigen(mesh_py, y_size, x_size);

    Eigen::MatrixXd rho_gradient = get_rho_gradient(mesh, y_vector, x_vector);

    return numy_interface::eigen_to_ndarray<double>(rho_gradient, {y_size, x_size});
}


pybind11::array_t<double>
get_rho_gradient_time_rho_py(const pybind11::array_t<double> &mesh_py, const pybind11::array_t<double> &x_vector_py, const pybind11::array_t<double> &y_vector_py) {
    size_t x_size = x_vector_py.request().size;
    size_t y_size = y_vector_py.request().size;

    Eigen::VectorXd x_vector = numy_interface::convert_py_to_eigen(x_vector_py, x_size);
    Eigen::VectorXd y_vector = numy_interface::convert_py_to_eigen(y_vector_py, y_size);
    Eigen::MatrixXd mesh = numy_interface::convert_py_to_eigen(mesh_py, y_size, x_size);

    Eigen::MatrixXd rho_gradient = get_rho_gradient_time_rho(mesh, y_vector, x_vector);

    return numy_interface::eigen_to_ndarray<double>(rho_gradient, {y_size, x_size});
}


