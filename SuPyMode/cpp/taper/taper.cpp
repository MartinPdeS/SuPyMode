#include "taper.h"


TaperSection AlphaProfile::get_constant_custom_section(
    double length,
    double radius,
    double start_z,
    std::size_t n_point) const
{
    auto z_array = linspace(start_z, start_z + length, n_point);
    std::vector<double> radius_array(n_point, radius);
    return TaperSection(std::move(z_array), std::move(radius_array));
}

void AlphaProfile::add_taper_segment(double alpha, double initial_heating_length, double stretching_length, size_t n_point) {
    add_taper_custom_segment(
        alpha,
        initial_heating_length,
        last_radius(),
        stretching_length,
        last_z(),
        n_point);
}

void AlphaProfile::add_end_of_taper_segment(std::size_t n_point) {
    if (section_list_.empty()) return;
    if (last_section().is_constant()) return;

    const double heating_length_final = last_section().heating_length_final();
    if (!std::isfinite(heating_length_final)) {
        throw std::runtime_error("add_end_of_taper_segment: last section has no heating_length_final.");
    }

    const double length = heating_length_final / 2.0;
    section_list_.push_back(get_constant_custom_section(length, last_radius(), last_z(), n_point));
    initialized_ = false;
}

// Interpolation over assembled sections (z query)
std::vector<double> AlphaProfile::compute_radius_from_segment_from_interpolation(const std::vector<double>& z) const {
    std::vector<double> radius_value(z.size(), 0.0);

    for (const auto& section : section_list_) {
        const auto interpolator = section.interpolation(false, 0.0);
        for (std::size_t i = 0; i < z.size(); ++i) {
            const double evaluation = interpolator(z[i]);
            if (evaluation != 0.0) radius_value[i] = evaluation;
        }
    }

    return radius_value;
}

std::vector<double> AlphaProfile::compute_adiabatic(const std::vector<double>& distance, const std::vector<double>& radius) const {
    if (distance.size() != radius.size()) throw std::invalid_argument("compute_adiabatic: size mismatch.");
    if (distance.size() < 3) throw std::invalid_argument("compute_adiabatic: need at least 3 points.");

    std::vector<double> log_radius(radius.size());
    for (std::size_t i = 0; i < radius.size(); ++i) {
        if (!(radius[i] > 0.0)) throw std::runtime_error("compute_adiabatic: non physical radius.");
        log_radius[i] = std::log(radius[i]);
    }

    const auto d_z = gradient_1d(distance, distance); // derivative of x wrt x yields ~1, but we want spacing.
    // The above is not right. We want dz along index. For linspace, dz is constant, but with general x:
    // Use y=distance and x=index is awkward. We will compute dz directly from distance:
    std::vector<double> dz(distance.size());
    {
        // dz here should represent derivative of z wrt index step ~1, but numpy returns spacing per point.
        // Easiest: compute gradient of z with x=index (uniform) then use it, but that equals local dz per step.
        std::vector<double> index(distance.size());
        for (std::size_t i = 0; i < distance.size(); ++i) index[i] = static_cast<double>(i);
        dz = gradient_1d(distance, index);
    }

    std::vector<double> ditr;
    {
        std::vector<double> index(radius.size());
        for (std::size_t i = 0; i < radius.size(); ++i) index[i] = static_cast<double>(i);
        ditr = gradient_1d(log_radius, index);
    }

    std::vector<double> adiabatic_value(distance.size());
    for (std::size_t i = 0; i < distance.size(); ++i) {
        adiabatic_value[i] = std::abs(ditr[i] / dz[i]);
    }
    return adiabatic_value;
}

std::tuple<std::vector<double>, double, double> AlphaProfile::compute_radius_from_segment(
    double alpha,
    double initial_heating_length,
    double stretching_length,
    double initial_radius,
    const std::vector<double>& distance) const
{
    assert_conditions(alpha, stretching_length, initial_heating_length);

    std::vector<double> radius(distance.size());

    const double term2 = (1.0 - alpha) * initial_heating_length;
    const double term3 = -1.0 / (2.0 * alpha);

    for (std::size_t i = 0; i < distance.size(); ++i) {
        const double term0 = 2.0 * alpha * distance[i];
        const double base = 1.0 + term0 / term2;
        radius[i] = initial_radius * std::pow(base, term3);
        if (!(radius[i] > 0.0)) {
            throw std::runtime_error("compute_radius_from_segment: non physical radius encountered.");
        }
    }

    const double final_radius =
        initial_radius * std::pow(1.0 + alpha * stretching_length / initial_heating_length, -1.0 / (2.0 * alpha));
    const double final_heating_length = initial_heating_length + alpha * stretching_length;

    return {radius, final_radius, final_heating_length};
}


void AlphaProfile::add_taper_custom_segment(
    double alpha,
    double initial_heating_length,
    double initial_radius,
    double stretching_length,
    double start_z,
    size_t n_point)
{
    if (alpha == 0.0) alpha = 0.01; // match Python safeguard

    const double z_0 = (1.0 - alpha) * stretching_length / 2.0;
    auto distance_local = linspace(0.0, z_0, n_point);

    if (distance_local.empty() || distance_local.front() != 0.0) {
        throw std::runtime_error("add_taper_custom_segment: distance must start at 0.");
    }

    auto [radius_local, final_radius, final_heating_length] = compute_radius_from_segment(
        alpha, initial_heating_length, stretching_length, initial_radius, distance_local);

    // Shift by start_z
    for (auto& z_value : distance_local) z_value += start_z;

    section_list_.push_back(
        TaperSection(std::move(distance_local), std::move(radius_local), initial_heating_length, final_heating_length));
    initialized_ = false;
}

void AlphaProfile::add_constant_segment(double length, size_t n_point) {
    const double start_z = last_z();
    const double radius_value = last_radius();
    section_list_.push_back(get_constant_custom_section(length, radius_value, start_z, n_point));
    initialized_ = false;
}

std::vector<double> AlphaProfile::compute_taper_angle(const std::vector<double>& distance, const std::vector<double>& radius) const {
    if (distance.size() != radius.size()) throw std::invalid_argument("compute_taper_angle: size mismatch.");
    if (distance.size() < 3) throw std::invalid_argument("compute_taper_angle: need at least 3 points.");

    std::vector<double> index(distance.size());
    for (std::size_t i = 0; i < distance.size(); ++i) index[i] = static_cast<double>(i);

    const auto d_z = gradient_1d(distance, index);
    const auto d_rho = gradient_1d(radius, index);

    std::vector<double> angle(distance.size());
    for (std::size_t i = 0; i < distance.size(); ++i) angle[i] = std::abs(d_rho[i] / d_z[i]);
    return angle;
}

std::vector<double> AlphaProfile::evaluate_adiabatic_factor(const std::vector<double>& itr) const {
    ensure_initialized_("evaluate_adiabatic_factor");

    if (itr_to_adiabatic_x_.empty() || itr_to_adiabatic_y_.empty()) {
        throw std::runtime_error("evaluate_adiabatic_factor: mapping not built. Call initialize().");
    }

    const double nan_value = std::numeric_limits<double>::quiet_NaN();
    LinearInterpolator interpolator(itr_to_adiabatic_x_, itr_to_adiabatic_y_, false, nan_value);
    return interpolator(itr);
}

double AlphaProfile::evaluate_itr_at_distance(double z) const {
    ensure_initialized_("evaluate_itr_at_distance");

    LinearInterpolator interpolator(distance_, itr_list_, false, 0.0);
    return interpolator(z);
}

std::vector<double> AlphaProfile::evaluate_distance_vs_itr(const std::vector<double>& itr_query) const {
    ensure_initialized_("evaluate_distance_vs_itr");

    if (itr_to_distance_x_.empty() || itr_to_distance_y_.empty()) {
        throw std::runtime_error("evaluate_distance_vs_itr: inverse mapping was not built. Call initialize().");
    }

    // bounds_error=true matches your Python intent.
    LinearInterpolator interpolator(itr_to_distance_x_, itr_to_distance_y_, true, 0.0);
    return interpolator(itr_query);
}

LinearInterpolator AlphaProfile::get_itr_vs_distance_interpolation() const {
    ensure_initialized_("get_itr_vs_distance_interpolation");
    return LinearInterpolator(distance_, itr_list_, false, 0.0);
}


void AlphaProfile::initialize() {
    if (add_end_of_taper_section_) {
        add_end_of_taper_segment();
    }

    const double z_end = last_z();
    if (!(z_end > 0.0)) {
        throw std::runtime_error("initialize: profile has no length. Add sections first.");
    }

    // Assemble on the half profile grid (0 -> z_end)
    std::vector<double> distance_local = linspace(0.0, z_end, n_point_);
    std::vector<double> radius_local = compute_radius_from_segment_from_interpolation(distance_local);

    std::vector<double> itr_local(radius_local.size());
    for (std::size_t i = 0; i < radius_local.size(); ++i) {
        itr_local[i] = radius_local[i] / initial_radius_;
    }

    std::vector<double> adiabatic_local = compute_adiabatic(distance_local, radius_local);
    std::vector<double> taper_angle_local = compute_taper_angle(distance_local, radius_local);

    // Write assembled arrays (full profile if symmetric)
    if (symmetric_) {
        const std::size_t n = distance_local.size();

        distance_.assign(2 * n - 1, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            distance_[i] = distance_local[i];
        }
        for (std::size_t i = 1; i < n; ++i) {
            distance_[n - 1 + i] =
                distance_local.back() + distance_local.back() - distance_local[n - 1 - i];
        }

        itr_list_.clear();
        radius_.clear();
        adiabatic_.clear();
        taper_angle_.clear();

        itr_list_.reserve(2 * n - 1);
        radius_.reserve(2 * n - 1);
        adiabatic_.reserve(2 * n - 1);
        taper_angle_.reserve(2 * n - 1);

        itr_list_.insert(itr_list_.end(), itr_local.begin(), itr_local.end());
        radius_.insert(radius_.end(), radius_local.begin(), radius_local.end());
        adiabatic_.insert(adiabatic_.end(), adiabatic_local.begin(), adiabatic_local.end());
        taper_angle_.insert(taper_angle_.end(), taper_angle_local.begin(), taper_angle_local.end());

        for (std::size_t i = 1; i < n; ++i) {
            itr_list_.push_back(itr_local[n - 1 - i]);
            radius_.push_back(radius_local[n - 1 - i]);
            adiabatic_.push_back(adiabatic_local[n - 1 - i]);
            taper_angle_.push_back(taper_angle_local[n - 1 - i]);
        }
    } else {
        distance_ = std::move(distance_local);
        radius_ = std::move(radius_local);
        itr_list_ = std::move(itr_local);
        adiabatic_ = std::move(adiabatic_local);
        taper_angle_ = std::move(taper_angle_local);
    }

    // Build monotonic helper mappings for any interpolation that uses ITR as x.
    // We choose the branch from start -> waist, because it is single-valued.
    if (itr_list_.size() < 2) {
        throw std::runtime_error("initialize: itr_list too small.");
    }

    const auto waist_iterator = std::min_element(itr_list_.begin(), itr_list_.end());
    const std::size_t waist_index =
        static_cast<std::size_t>(std::distance(itr_list_.begin(), waist_iterator));

    if (waist_index < 1) {
        throw std::runtime_error("initialize: waist_index invalid for monotonic mapping.");
    }

    auto build_monotonic_mapping = [](
        const std::vector<double>& x_source,
        const std::vector<double>& y_source,
        std::size_t end_index_inclusive,
        std::vector<double>& x_out,
        std::vector<double>& y_out)
    {
        x_out.assign(x_source.begin(), x_source.begin() + end_index_inclusive + 1);
        y_out.assign(y_source.begin(), y_source.begin() + end_index_inclusive + 1);

        std::reverse(x_out.begin(), x_out.end());
        std::reverse(y_out.begin(), y_out.end());

        std::vector<double> x_filtered;
        std::vector<double> y_filtered;
        x_filtered.reserve(x_out.size());
        y_filtered.reserve(y_out.size());

        for (std::size_t i = 0; i < x_out.size(); ++i) {
            if (x_filtered.empty() || x_out[i] > x_filtered.back()) {
                x_filtered.push_back(x_out[i]);
                y_filtered.push_back(y_out[i]);
            }
        }

        x_out = std::move(x_filtered);
        y_out = std::move(y_filtered);

        if (x_out.size() < 2) {
            throw std::runtime_error("initialize: monotonic mapping has fewer than 2 points.");
        }
    };

    build_monotonic_mapping(itr_list_, distance_, waist_index, itr_to_distance_x_, itr_to_distance_y_);
    build_monotonic_mapping(itr_list_, adiabatic_, waist_index, itr_to_adiabatic_x_, itr_to_adiabatic_y_);

    initialized_ = true;
}

void AlphaProfile::ensure_initialized_(const char* field_name) const {
    if (!initialized_) {
        throw std::runtime_error(std::string("Profile has not been initialized yet. Call initialize() before accessing ") + field_name + ".");
    }
}

// -
