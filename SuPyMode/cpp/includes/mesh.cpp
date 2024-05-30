#include <utility>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include "numpy_interface.cpp"
#include <unsupported/Eigen/MatrixFunctions>


std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_gradient_2p(const Eigen::MatrixXd& image, double dx, double dy) {
    int rows = image.rows();
    int cols = image.cols();

    // Initialize matrices to store the gradients in x and y directions
    Eigen::MatrixXd gradient_x(rows, cols);
    Eigen::MatrixXd gradient_y(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Compute dxVal
            double dxVal = 0.0;
            if (j > 0 && j + 1 < cols) {
                dxVal = (image(i, j + 1) - image(i, j - 1)) / (2 * dx);
            } else if (j == 0 && j + 2 < cols) {  // Use forward difference at the left boundary
                dxVal = (-3 * image(i, j) + 4 * image(i, j + 1) - image(i, j + 2)) / (2 * dx);
            } else if (j == cols - 1 && j - 2 >= 0) {  // Use backward difference at the right boundary
                dxVal = (3 * image(i, j) - 4 * image(i, j - 1) + image(i, j - 2)) / (2 * dx);
            } else if (j == 0 && j + 1 < cols) {  // Forward difference when there is only one neighboring point
                dxVal = (image(i, j + 1) - image(i, j)) / dx;
            } else if (j == cols - 1 && j - 1 >= 0) {  // Backward difference when there is only one neighboring point
                dxVal = (image(i, j) - image(i, j - 1)) / dx;
            }

            // Compute dyVal
            double dyVal = 0.0;
            if (i > 0 && i + 1 < rows) {
                dyVal = (image(i + 1, j) - image(i - 1, j)) / (2 * dy);
            } else if (i == 0 && i + 2 < rows) {  // Use forward difference at the top boundary
                dyVal = (-3 * image(i, j) + 4 * image(i + 1, j) - image(i + 2, j)) / (2 * dy);
            } else if (i == rows - 1 && i - 2 >= 0) {  // Use backward difference at the bottom boundary
                dyVal = (3 * image(i, j) - 4 * image(i - 1, j) + image(i - 2, j)) / (2 * dy);
            } else if (i == 0 && i + 1 < rows) {  // Forward difference when there is only one neighboring point
                dyVal = (image(i + 1, j) - image(i, j)) / dy;
            } else if (i == rows - 1 && i - 1 >= 0) {  // Backward difference when there is only one neighboring point
                dyVal = (image(i, j) - image(i - 1, j)) / dy;
            }

            gradient_x(i, j) = dxVal;
            gradient_y(i, j) = dyVal;
        }
    }

    return {gradient_x, gradient_y};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_gradient_5p(const Eigen::MatrixXd& image, const double dx, const double dy) {
    int rows = image.rows();
    int cols = image.cols();

    // Initialize matrices to store the gradients
    Eigen::MatrixXd gradient_x(rows, cols);
    Eigen::MatrixXd gradient_y(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double dxVal = 0.0;
            double dyVal = 0.0;

            // Calculate dxVal
            if (j >= 2 && j + 2 < cols) {
                dxVal = (image(i, j - 2) - 8 * image(i, j - 1) + 8 * image(i, j + 1) - image(i, j + 2)) / (12 * dx);
            } else if (j > 0 && j + 1 < cols) {
                dxVal = (image(i, j + 1) - image(i, j - 1)) / (2 * dx);
            } else if (j == 0 && j + 1 < cols) {
                dxVal = (image(i, j + 1) - image(i, j)) / dx;
            } else if (j == cols - 1 && j - 1 >= 0) {
                dxVal = (image(i, j) - image(i, j - 1)) / dx;
            }

            // Calculate dyVal
            if (i >= 2 && i + 2 < rows) {
                dyVal = (image(i - 2, j) - 8 * image(i - 1, j) + 8 * image(i + 1, j) - image(i + 2, j)) / (12 * dy);
            } else if (i > 0 && i + 1 < rows) {
                dyVal = (image(i + 1, j) - image(i - 1, j)) / (2 * dy);
            } else if (i == 0 && i + 1 < rows) {
                dyVal = (image(i + 1, j) - image(i, j)) / dy;
            } else if (i == rows - 1 && i - 1 >= 0) {
                dyVal = (image(i, j) - image(i - 1, j)) / dy;
            }

            gradient_x(i, j) = dxVal;
            gradient_y(i, j) = dyVal;
        }
    }

    return {gradient_x, gradient_y};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_gradient_7p(const Eigen::MatrixXd& image, double dx, double dy) {
    int rows = image.rows();
    int cols = image.cols();

    // Initialize matrices to store the gradients
    Eigen::MatrixXd gradient_x(rows, cols);
    Eigen::MatrixXd gradient_y(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double dxVal = 0.0;
            if (j >= 3 && j + 3 < cols) {
                dxVal = (-image(i, j - 3) + 9 * image(i, j - 2) - 45 * image(i, j - 1) + 45 * image(i, j + 1) - 9 * image(i, j + 2) + image(i, j + 3)) / (60 * dx);
            } else if (j > 0 && j + 1 < cols) {
                dxVal = (image(i, j + 1) - image(i, j - 1)) / (2 * dx);  // Central difference
            } else if (j == 0 && j + 2 < cols) {  // Forward difference at the left boundary
                dxVal = (-3 * image(i, j) + 4 * image(i, j + 1) - image(i, j + 2)) / (2 * dx);
            } else if (j == cols - 1 && j - 2 >= 0) {  // Backward difference at the right boundary
                dxVal = (3 * image(i, j) - 4 * image(i, j - 1) + image(i, j - 2)) / (2 * dx);
            } else if (j == 0 && j + 1 < cols) {  // Forward difference with one neighboring point
                dxVal = (image(i, j + 1) - image(i, j)) / dx;
            } else if (j == cols - 1 && j - 1 >= 0) {  // Backward difference with one neighboring point
                dxVal = (image(i, j) - image(i, j - 1)) / dx;
            }

            double dyVal = 0.0;
            if (i >= 3 && i + 3 < rows) {
                dyVal = (-image(i - 3, j) + 9 * image(i - 2, j) - 45 * image(i - 1, j) + 45 * image(i + 1, j) - 9 * image(i + 2, j) + image(i + 3, j)) / (60 * dy);
            } else if (i > 0 && i + 1 < rows) {
                dyVal = (image(i + 1, j) - image(i - 1, j)) / (2 * dy);  // Central difference
            } else if (i == 0 && i + 2 < rows) {  // Forward difference at the top boundary
                dyVal = (-3 * image(i, j) + 4 * image(i + 1, j) - image(i + 2, j)) / (2 * dy);
            } else if (i == rows - 1 && i - 2 >= 0) {  // Backward difference at the bottom boundary
                dyVal = (3 * image(i, j) - 4 * image(i - 1, j) + image(i - 2, j)) / (2 * dy);
            } else if (i == 0 && i + 1 < rows) {  // Forward difference with one neighboring point
                dyVal = (image(i + 1, j) - image(i, j)) / dy;
            } else if (i == rows - 1 && i - 1 >= 0) {  // Backward difference with one neighboring point
                dyVal = (image(i, j) - image(i - 1, j)) / dy;
            }

            gradient_x(i, j) = dxVal;
            gradient_y(i, j) = dyVal;
        }
    }

    return {gradient_x, gradient_y};
}


Eigen::MatrixXd get_rho_gradient_time_rho(const Eigen::MatrixXd &mesh, const Eigen::VectorXd &x_vector, const Eigen::VectorXd &y_vector) {
    size_t x_size = x_vector.size();
    size_t y_size = y_vector.size();

    double dx = x_vector(1) - x_vector(0);
    double dy = y_vector(1) - y_vector(0);

    auto gradients = compute_gradient_5p(mesh, dy, dx);
    Eigen::MatrixXd gradient_x = gradients.first;
    Eigen::MatrixXd gradient_y = gradients.second;

    Eigen::MatrixXd rho_gradient(x_size, y_size);

    double angle, cos_angle, sin_angle;

    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            double x = x_vector(i);
            double y = y_vector(j);
            double rho = std::sqrt(pow(x, 2) + pow(y, 2));

            angle = std::atan2(x, y);
            cos_angle = std::cos(angle);
            sin_angle = std::sin(angle);
            rho_gradient(i, j) = (gradient_x(i, j) * cos_angle + gradient_y(i, j) * sin_angle) * rho;
       }
    }

    return rho_gradient;
}

Eigen::MatrixXd get_rho_gradient(const Eigen::MatrixXd &mesh, const Eigen::VectorXd &x_vector, const Eigen::VectorXd &y_vector) {
    size_t x_size = x_vector.size();
    size_t y_size = y_vector.size();

    double dx = x_vector(1) - x_vector(0);
    double dy = y_vector(1) - y_vector(0);

    auto gradients = compute_gradient_5p(mesh, dy, dx);
    Eigen::MatrixXd gradient_x = gradients.first;
    Eigen::MatrixXd gradient_y = gradients.second;

    Eigen::MatrixXd rho_gradient(x_size, y_size);

    double angle, cos_angle, sin_angle;

    for (int i = 0; i < x_size; ++i) {
        for (int j = 0; j < y_size; ++j) {
            double x = x_vector(i);
            double y = y_vector(j);

            angle = std::atan2(x, y);
            cos_angle = std::cos(angle);
            sin_angle = std::sin(angle);
            rho_gradient(i, j) = (gradient_x(i, j) * cos_angle + gradient_y(i, j) * sin_angle);
       }
    }

    return rho_gradient;
}



pybind11::array_t<double> get_rho_gradient_py(const pybind11::array_t<double> &mesh_py, const pybind11::array_t<double> &x_vector_py, const pybind11::array_t<double> &y_vector_py) {
    size_t x_size = x_vector_py.request().size;
    size_t y_size = y_vector_py.request().size;

    Eigen::VectorXd x_vector = convert_py_to_eigen(x_vector_py, x_size);
    Eigen::VectorXd y_vector = convert_py_to_eigen(y_vector_py, y_size);
    Eigen::MatrixXd mesh = convert_py_to_eigen(mesh_py, x_size, y_size);

    Eigen::MatrixXd rho_gradient = get_rho_gradient(mesh, x_vector, y_vector);

    return eigen_to_ndarray<double>(rho_gradient, {y_size, x_size});
}

pybind11::array_t<double> get_rho_gradient_time_rho_py(const pybind11::array_t<double> &mesh_py, const pybind11::array_t<double> &x_vector_py, const pybind11::array_t<double> &y_vector_py) {
    size_t x_size = x_vector_py.request().size;
    size_t y_size = y_vector_py.request().size;

    Eigen::VectorXd x_vector = convert_py_to_eigen(x_vector_py, x_size);
    Eigen::VectorXd y_vector = convert_py_to_eigen(y_vector_py, y_size);
    Eigen::MatrixXd mesh = convert_py_to_eigen(mesh_py, x_size, y_size);

    Eigen::MatrixXd rho_gradient = get_rho_gradient_time_rho(mesh, x_vector, y_vector);

    return eigen_to_ndarray<double>(rho_gradient, {y_size, x_size});
}


PYBIND11_MODULE(Example, module)
{
    module.def("get_rho_gradient_5p", &get_rho_gradient_py, pybind11::arg("mesh"), pybind11::arg("x_vector"), pybind11::arg("y_vector"));
    module.def("get_rho_gradient_time_rho_5p", &get_rho_gradient_time_rho_py, pybind11::arg("mesh"), pybind11::arg("x_vector"), pybind11::arg("y_vector"));
}
