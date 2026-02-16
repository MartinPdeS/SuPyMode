#pragma once

#include <vector>
#include <stdexcept>
#include <cstddef>

// -----------------------------
// Small numeric helpers
// -----------------------------
inline std::vector<double> linspace(double start, double stop, std::size_t number_of_points) {
    if (number_of_points == 0) return {};
    if (number_of_points == 1) return {start};

    std::vector<double> values(number_of_points);
    const double step = (stop - start) / static_cast<double>(number_of_points - 1);

    for (std::size_t i = 0; i < number_of_points; ++i) {
        values[i] = start + step * static_cast<double>(i);
    }
    return values;
}

// 1D gradient with edge_order=2 like numpy.gradient(..., edge_order=2) for uniform spacing in index,
// but we allow nonuniform x by dividing by dx at each point.
inline std::vector<double> gradient_1d(const std::vector<double>& y, const std::vector<double>& x) {
    const std::size_t n = y.size();
    if (x.size() != n) throw std::invalid_argument("gradient: x and y must have same size.");
    if (n < 3) throw std::invalid_argument("gradient: need at least 3 points for edge_order=2.");

    std::vector<double> dy_dx(n);

    // Interior: central difference
    for (std::size_t i = 1; i + 1 < n; ++i) {
        const double dx = x[i + 1] - x[i - 1];
        if (dx == 0.0) throw std::runtime_error("gradient: repeated x value.");
        dy_dx[i] = (y[i + 1] - y[i - 1]) / dx;
    }

    // Edges: second order one sided
    // i = 0
    {
        const double dx0 = x[1] - x[0];
        const double dx1 = x[2] - x[1];
        if (dx0 == 0.0 || dx1 == 0.0) throw std::runtime_error("gradient: repeated x value.");
        // Use a local quadratic fit style formula on nonuniform grid is more complex.
        // Here we assume approximately uniform spacing as in typical linspace usage.
        // For linspace this matches numpy edge_order=2 behavior closely.
        const double dx = x[1] - x[0];
        dy_dx[0] = (-3.0 * y[0] + 4.0 * y[1] - 1.0 * y[2]) / (2.0 * dx);
    }
    // i = n-1
    {
        const double dx = x[n - 1] - x[n - 2];
        if (dx == 0.0) throw std::runtime_error("gradient: repeated x value.");
        dy_dx[n - 1] = (3.0 * y[n - 1] - 4.0 * y[n - 2] + 1.0 * y[n - 3]) / (2.0 * dx);
    }

    return dy_dx;
}
