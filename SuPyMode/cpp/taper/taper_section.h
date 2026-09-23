// alpha_profile.h
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>
#include <utils/math.h>


// -----------------------------
// TaperSection
// -----------------------------
class TaperSection {
public:
    TaperSection() = default;

    TaperSection(
        std::vector<double> z_array,
        std::vector<double> radius_array,
        double heating_length_initial = std::numeric_limits<double>::quiet_NaN(),
        double heating_length_final = std::numeric_limits<double>::quiet_NaN())
        : z_array_(std::move(z_array)),
          radius_array_(std::move(radius_array)),
          heating_length_initial_(heating_length_initial),
          heating_length_final_(heating_length_final)
    {
        if (z_array_.size() != radius_array_.size()) {
            throw std::invalid_argument("TaperSection: z_array and radius_array size mismatch.");
        }
        if (z_array_.empty()) throw std::invalid_argument("TaperSection: arrays must not be empty.");
        if (!std::is_sorted(z_array_.begin(), z_array_.end())) {
            throw std::invalid_argument("TaperSection: z_array must be sorted ascending.");
        }
    }

    double z_initial() const { return z_array_.front(); }
    double z_final() const { return z_array_.back(); }
    double radius_initial() const { return radius_array_.front(); }
    double radius_final() const { return radius_array_.back(); }
    bool is_constant() const { return radius_array_.front() == radius_array_.back(); }

    double heating_length_initial() const { return heating_length_initial_; }
    double heating_length_final() const { return heating_length_final_; }

    const std::vector<double>& z_array() const { return z_array_; }
    const std::vector<double>& radius_array() const { return radius_array_; }

    // Similar to scipy interp1d(bounds_error=False, fill_value=0)
    LinearInterpolator interpolation(bool bounds_error = false, double fill_value = 0.0) const {
        return LinearInterpolator(z_array_, radius_array_, bounds_error, fill_value);
    }

private:
    std::vector<double> z_array_;
    std::vector<double> radius_array_;
    double heating_length_initial_ = std::numeric_limits<double>::quiet_NaN();
    double heating_length_final_ = std::numeric_limits<double>::quiet_NaN();
};
