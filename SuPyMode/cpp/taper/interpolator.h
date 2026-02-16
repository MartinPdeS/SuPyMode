// interpolator_1d.h
#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

class Interpolator1D {
public:
    Interpolator1D() = default;

    Interpolator1D(
        std::vector<double> x_values,
        std::vector<double> y_values,
        bool bounds_error,
        double fill_value
    )
        : x_(std::move(x_values)),
          y_(std::move(y_values)),
          bounds_error_(bounds_error),
          fill_value_(fill_value)
    {
        if (x_.size() != y_.size()) throw std::invalid_argument("Interpolator1D: x and y size mismatch.");
        if (x_.size() < 2) throw std::invalid_argument("Interpolator1D: need at least 2 points.");
        if (!std::is_sorted(x_.begin(), x_.end())) throw std::invalid_argument("Interpolator1D: x must be sorted ascending.");
    }

    double operator()(double x_query) const {
        const double x_min = x_.front();
        const double x_max = x_.back();

        if (x_query < x_min || x_query > x_max) {
            if (bounds_error_) throw std::out_of_range("Interpolator1D: query out of bounds.");
            return fill_value_;
        }

        if (x_query == x_min) return y_.front();
        if (x_query == x_max) return y_.back();

        auto it = std::upper_bound(x_.begin(), x_.end(), x_query);
        const std::size_t right_index = static_cast<std::size_t>(std::distance(x_.begin(), it));
        const std::size_t left_index = right_index - 1;

        const double x0 = x_[left_index];
        const double x1 = x_[right_index];
        const double y0 = y_[left_index];
        const double y1 = y_[right_index];

        const double dx = x1 - x0;
        if (dx == 0.0) return y0;

        const double t = (x_query - x0) / dx;
        return y0 + t * (y1 - y0);
    }

    std::vector<double> operator()(const std::vector<double>& x_queries) const {
        std::vector<double> y_queries(x_queries.size());
        for (std::size_t i = 0; i < x_queries.size(); ++i) y_queries[i] = (*this)(x_queries[i]);
        return y_queries;
    }

    bool bounds_error() const { return bounds_error_; }
    double fill_value() const { return fill_value_; }
    double x_min() const { return x_.front(); }
    double x_max() const { return x_.back(); }

private:
    std::vector<double> x_;
    std::vector<double> y_;
    bool bounds_error_ = false;
    double fill_value_ = 0.0;
};


// Linear interpolation similar to scipy interp1d(kind="linear") with
// bounds_error configurable and fill_value outside bounds.
class LinearInterpolator {
public:
    LinearInterpolator() = default;

    LinearInterpolator(
        std::vector<double> x_values,
        std::vector<double> y_values,
        bool bounds_error,
        double fill_value)
        : x_(std::move(x_values)),
          y_(std::move(y_values)),
          bounds_error_(bounds_error),
          fill_value_(fill_value)
    {
        if (x_.size() != y_.size()) throw std::invalid_argument("Interpolator: x and y size mismatch.");
        if (x_.size() < 2) throw std::invalid_argument("Interpolator: need at least 2 points.");
        if (!std::is_sorted(x_.begin(), x_.end())) {
            throw std::invalid_argument("Interpolator: x must be sorted ascending.");
        }
    }

    double operator()(double x_query) const {
        const double x_min = x_.front();
        const double x_max = x_.back();

        if (x_query < x_min || x_query > x_max) {
            if (bounds_error_) throw std::out_of_range("Interpolator: query out of bounds.");
            return fill_value_;
        }

        // Handle exact endpoints quickly
        if (x_query == x_min) return y_.front();
        if (x_query == x_max) return y_.back();

        // Find rightmost index such that x_[index] <= x_query
        auto it = std::upper_bound(x_.begin(), x_.end(), x_query);
        const std::size_t right_index = static_cast<std::size_t>(std::distance(x_.begin(), it));
        const std::size_t left_index = right_index - 1;

        const double x0 = x_[left_index];
        const double x1 = x_[right_index];
        const double y0 = y_[left_index];
        const double y1 = y_[right_index];

        const double dx = x1 - x0;
        if (dx == 0.0) return y0;

        const double t = (x_query - x0) / dx;
        return y0 + t * (y1 - y0);
    }

    std::vector<double> operator()(const std::vector<double>& x_queries) const {
        std::vector<double> y_queries(x_queries.size());
        for (std::size_t i = 0; i < x_queries.size(); ++i) y_queries[i] = (*this)(x_queries[i]);
        return y_queries;
    }

private:
    std::vector<double> x_;
    std::vector<double> y_;
    bool bounds_error_ = false;
    double fill_value_ = 0.0;
};
