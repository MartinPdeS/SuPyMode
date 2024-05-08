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
            // Use central differences when possible
            double dxVal = 0.0;
            if (j > 0 && j + 1 < cols) {
                dxVal = (image(i, j + 1) - image(i, j - 1)) / (2 * dx);
            } else if (j == 0) {  // Use forward difference at the left boundary
                dxVal = (-3 * image(i, j) + 4 * image(i, j + 1) - image(i, j + 2)) / (2 * dx);
            } else if (j == cols - 1) {  // Use backward difference at the right boundary
                dxVal = (3 * image(i, j) - 4 * image(i, j - 1) + image(i, j - 2)) / (2 * dx);
            }

            double dyVal = 0.0;
            if (i > 0 && i + 1 < rows) {
                dyVal = (image(i + 1, j) - image(i - 1, j)) / (2 * dy);
            } else if (i == 0) {  // Use forward difference at the top boundary
                dyVal = (-3 * image(i, j) + 4 * image(i + 1, j) - image(i + 2, j)) / (2 * dy);
            } else if (i == rows - 1) {  // Use backward difference at the bottom boundary
                dyVal = (3 * image(i, j) - 4 * image(i - 1, j) + image(i - 2, j)) / (2 * dy);
            }

            gradient_x(i, j) = dxVal;
            gradient_y(i, j) = dyVal;
        }
    }

    return {gradient_x, gradient_y};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> compute_gradient_5p(const Eigen::MatrixXd& image, double dx, double dy) {
    int rows = image.rows();
    int cols = image.cols();

    // Initialize matrices to store the gradients
    Eigen::MatrixXd gradient_x(rows, cols);
    Eigen::MatrixXd gradient_y(rows, cols);

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double dxVal = 0.0;
            if (j >= 2 && j + 2 < cols) {
                dxVal = (image(i, j - 2) - 8 * image(i, j - 1) + 8 * image(i, j + 1) - image(i, j + 2)) / (12 * dx);
            } else {
                // Apply lower-order schemes at the boundaries
                if (j < 2 || j + 2 >= cols) {
                    dxVal = (j > 0 && j + 1 < cols) ? (image(i, j + 1) - image(i, j - 1)) / (2 * dx) : 0; // Fallback to central difference where possible
                }
            }

            double dyVal = 0.0;
            if (i >= 2 && i + 2 < rows) {
                dyVal = (image(i - 2, j) - 8 * image(i - 1, j) + 8 * image(i + 1, j) - image(i + 2, j)) / (12 * dy);
            } else {
                // Apply lower-order schemes at the boundaries
                if (i < 2 || i + 2 >= rows) {
                    dyVal = (i > 0 && i + 1 < rows) ? (image(i + 1, j) - image(i - 1, j)) / (2 * dy) : 0; // Fallback to central difference where possible
                }
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
            } else {
                // Apply lower-order schemes at the boundaries
                if (j < 3 || j + 3 >= cols) {
                    dxVal = (j > 0 && j + 1 < cols) ? (image(i, j + 1) - image(i, j - 1)) / (2 * dx) : 0; // Fallback to central difference where possible
                }
            }

            double dyVal = 0.0;
            if (i >= 3 && i + 3 < rows) {
                dyVal = (-image(i - 3, j) + 9 * image(i - 2, j) - 45 * image(i - 1, j) + 45 * image(i + 1, j) - 9 * image(i + 2, j) + image(i + 3, j)) / (60 * dy);
            } else {
                // Apply lower-order schemes at the boundaries
                if (i < 3 || i + 3 >= rows) {
                    dyVal = (i > 0 && i + 1 < rows) ? (image(i + 1, j) - image(i - 1, j)) / (2 * dy) : 0; // Fallback to central difference where possible
                }
            }

            gradient_x(i, j) = dxVal;
            gradient_y(i, j) = dyVal;
        }
    }

    return {gradient_x, gradient_y};
}


pybind11::array_t<double> get_example(pybind11::array_t<double> &mesh_py, pybind11::array_t<double> &x_space_py, pybind11::array_t<double> &y_space_py) {

    Eigen::VectorXd x_space = convert_py_to_eigen<double>(x_space_py, x_space_py.request().size);
    Eigen::VectorXd y_space = convert_py_to_eigen<double>(y_space_py, y_space_py.request().size);
    Eigen::MatrixXd image = convert_py_to_eigen<double>(mesh_py);

    size_t x_size = x_space.size();
    size_t y_size = y_space.size();
    double dx = x_space(1) - x_space(0);
    double dy = y_space(1) - y_space(0);

    auto gradients = compute_gradient_5p(image, dx, dy);
    Eigen::MatrixXd gradient_x = gradients.first;
    Eigen::MatrixXd gradient_y = gradients.second;

    Eigen::MatrixXd rho_gradient(x_size, y_size);

    double angle, cos_angle, sin_angle;

    for (int x_idx = 0; x_idx < x_size; ++x_idx) {
        for (int y_idx = 0; y_idx < y_size; ++y_idx) {
            double x = x_space(x_idx);
            double y = y_space(y_idx);

            angle = std::atan2(x, y);
            cos_angle = std::cos(angle);
            sin_angle = std::sin(angle);
            rho_gradient(x_idx, y_idx) = gradient_x(x_idx, y_idx) * cos_angle + gradient_y(x_idx, y_idx) * sin_angle;
       }
    }

    return eigen_to_ndarray<double>(rho_gradient, {x_size, y_size});
}


PYBIND11_MODULE(Example, module)
{
    module.def("get_rho_gradient_5p", &get_example, pybind11::arg("mesh"), pybind11::arg("x_space"), pybind11::arg("y_space"));
}
