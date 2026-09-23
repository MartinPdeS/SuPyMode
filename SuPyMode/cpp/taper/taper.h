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
#include "interpolator.h"
#include "taper_section.h"

class AlphaProfile {
public:
    AlphaProfile(
        double initial_radius = 1.0,
        std::size_t n_point = 200,
        bool symmetric = false,
        std::string label = "profile",
        bool add_end_of_taper_section = true)
        : initial_radius_(initial_radius),
          n_point_(n_point),
          symmetric_(symmetric),
          label_(std::move(label)),
          add_end_of_taper_section_(add_end_of_taper_section)
    {}

    // Section access
    const TaperSection& first_section() const {
        if (section_list_.empty()) throw std::runtime_error("first_section: no sections.");
        return section_list_.front();
    }

    const TaperSection& last_section() const {
        if (section_list_.empty()) throw std::runtime_error("last_section: no sections.");
        return section_list_.back();
    }

    // State getters (valid after initialize)
    const std::vector<double>& distance() const {
        ensure_initialized_("distance");
        return distance_;
    }

    const std::vector<double>& radius() const {
        ensure_initialized_("radius");
        return radius_;
    }

    const std::vector<double>& itr_list() const {
        ensure_initialized_("itr_list");
        return itr_list_;
    }

    const std::vector<double>& adiabatic() const {
        ensure_initialized_("adiabatic");
        return adiabatic_;
    }

    const std::vector<double>& taper_angle() const {
        ensure_initialized_("taper_angle");
        return taper_angle_;
    }

    double initial_radius() const { return initial_radius_; }
    std::size_t n_point() const { return n_point_; }
    bool symmetric() const { return symmetric_; }
    const std::string& label() const { return label_; }

    double last_z() const {
        if (section_list_.empty()) return 0.0;
        return last_section().z_final();
    }

    double total_length() const {
        if (section_list_.empty()) return 0.0;
        return last_section().z_final();
    }

    double last_radius() const {
        if (section_list_.empty()) return initial_radius_;
        return last_section().radius_final();
    }

    double smallest_itr() const {
        ensure_initialized_("smallest_itr");
        if (itr_list_.empty()) throw std::runtime_error("smallest_itr: empty itr_list.");
        return *std::min_element(itr_list_.begin(), itr_list_.end());
    }

    // Sections construction
    TaperSection get_constant_custom_section(double length, double radius, double start_z = 0.0, size_t n_point = 100) const;

    void add_constant_segment(double length, std::size_t n_point = 100);

    Interpolator1D get_itr_vs_distance_interpolation_object() const {
        ensure_initialized_("get_itr_vs_distance_interpolation");
        return Interpolator1D(distance_, itr_list_, false, 0.0);
    }

    Interpolator1D get_distance_vs_itr_interpolation_object() const {
        ensure_initialized_("get_distance_vs_itr_interpolation");
        return Interpolator1D(itr_to_distance_x_, itr_to_distance_y_, true, 0.0);
    }

    void add_end_of_taper_segment(std::size_t n_point = 100);

    // Core segment formula (returns radius array, final_radius, final_heating_length)
    std::tuple<std::vector<double>, double, double> compute_radius_from_segment(
        double alpha,
        double initial_heating_length,
        double stretching_length,
        double initial_radius,
        const std::vector<double>& distance) const;

    void assert_conditions(double alpha, double stretching_length, double initial_heating_length) const {
        if (!(initial_heating_length > 0.0)) {
            throw std::invalid_argument("Initial heating length must be positive.");
        }
        if (alpha == 0.0) {
            throw std::invalid_argument("Alpha must not be zero to avoid division by zero.");
        }
        if (alpha < 0.0 && stretching_length >= initial_heating_length / std::abs(alpha)) {
            throw std::invalid_argument("Stretching length for negative alpha exceeds the viable limit.");
        }
    }

    void add_taper_custom_segment(
        double alpha,
        double initial_heating_length,
        double initial_radius,
        double stretching_length,
        double start_z = 0.0,
        std::size_t n_point = 100);


    void add_taper_segment(double alpha, double initial_heating_length, double stretching_length, size_t n_point = 100);

    // Interpolation over assembled sections (z query)
    std::vector<double> compute_radius_from_segment_from_interpolation(const std::vector<double>& z) const;

    std::vector<double> compute_adiabatic(const std::vector<double>& distance, const std::vector<double>& radius) const;

    std::vector<double> compute_taper_angle(const std::vector<double>& distance, const std::vector<double>& radius) const;

    // Evaluation helpers similar to your Python methods
    std::vector<double> evaluate_adiabatic_factor(const std::vector<double>& itr) const;

    double evaluate_itr_at_distance(double z) const;

    std::vector<double> evaluate_distance_vs_itr(const std::vector<double>& itr_query) const;

    LinearInterpolator get_itr_vs_distance_interpolation() const;

    void initialize();



private:
    void ensure_initialized_(const char* field_name) const;

private:
    double initial_radius_ = 1.0;
    std::size_t n_point_ = 200;
    bool symmetric_ = false;
    std::string label_ = "profile";
    bool add_end_of_taper_section_ = true;

    std::vector<TaperSection> section_list_;

    bool initialized_ = false;
    std::vector<double> distance_;
    std::vector<double> radius_;
    std::vector<double> itr_list_;
    std::vector<double> adiabatic_;
    std::vector<double> taper_angle_;
    std::vector<double> itr_to_distance_x_;
    std::vector<double> itr_to_distance_y_;
    std::vector<double> itr_to_adiabatic_x_;
    std::vector<double> itr_to_adiabatic_y_;

};
