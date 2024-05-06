#include <iostream>
#include <Eigen/Dense>
#include <pybind11/pybind11.h>
#include "numpy_interface.cpp"
#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;

// Function to compute the gradients in x and y directions
std::pair<MatrixXd, MatrixXd> compute_gradient_2p(const MatrixXd& image, double dx = 1, double dy = 1) {
    int rows = image.rows();
    int cols = image.cols();

    // Compute gradients in x and y directions
    MatrixXd gradientX(rows, cols);
    MatrixXd gradientY(rows, cols);

    // Compute central differences for gradients in x and y directions
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double dxVal = (j + 1 < cols ? (image(i, j + 1) - image(i, j)) / dx : 0.0);
            double dyVal = (i + 1 < rows ? (image(i + 1, j) - image(i, j)) / dy : 0.0);

            gradientX(i, j) = dxVal;
            gradientY(i, j) = dyVal;
        }
    }

    return {gradientX, gradientY};
}

std::pair<MatrixXd, MatrixXd> compute_gradient_5p(const MatrixXd &image, double dx = 1, double dy = 1) {
    int rows = image.rows();
    int cols = image.cols();

    // Compute gradients in x and y directions
    MatrixXd gradientX(rows, cols);
    MatrixXd gradientY(rows, cols);

    // Compute gradients using the five-point stencil for second-order accuracy
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            // Central differences for x gradient
            double dxVal = (j >= 2 && j + 2 < cols ? (image(i, j - 2) - 8 * image(i, j - 1) + 8 * image(i, j + 1) - image(i, j + 2)) / (12 * dx) : 0.0);
            // Central differences for y gradient
            double dyVal = (i >= 2 && i + 2 < rows ? (image(i - 2, j) - 8 * image(i - 1, j) + 8 * image(i + 1, j) - image(i + 2, j)) / (12 * dy) : 0.0);

            gradientX(i, j) = dxVal;
            gradientY(i, j) = dyVal;
        }
    }

    return {gradientX, gradientY};
}



pybind11::array_t<double> get_example() {
    // Example circular image represented by an Eigen matrix
    size_t size = 100;
    MatrixXd image(size, size);
    image.setZero();

    // Draw a circle in the image
    double radius = size / 4.0;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double dx = j - size / 2.0;
            double dy = i - size / 2.0;
            if (dx * dx + dy * dy <= radius * radius) {
                image(i, j) = 1.0; // inside the circle
            }
        }
    }

    // Compute gradients
    auto gradients = computeGradient(image);
    MatrixXd gradientX = gradients.first;
    MatrixXd gradientY = gradients.second;

    MatrixXd rho_gradient(size, size);

    double angle, cos_angle, sin_angle;

    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            angle = std::atan2(i - 50, j - 50);
            cos_angle = std::cos(angle);
            sin_angle = std::sin(angle);
            rho_gradient(i, j) = gradientX(i, j) * cos_angle + gradientY(i, j) * sin_angle;
       }
    }

    return eigen_to_ndarray<double>(rho_gradient, {size, size});
}


PYBIND11_MODULE(Example, module)
{
    module.def("get_example", &get_example, "Test method");
}
