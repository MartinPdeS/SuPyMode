#pragma once

#include <string>
#include <vector>
#include <tuple>
#include <array>
#include <utility>
#include <stdexcept>

// Uses your existing types
// enum class BoundaryCondition { Zero, Symmetric, AntiSymmetric };
// class Boundaries { public: BoundaryCondition left, right, top, bottom; ... };

class ModeLabeler {
public:
    ModeLabeler() = default;

    std::string get(const Boundaries& boundaries, const size_t mode_number) const
    {
        const BoundaryCondition x_parity = get_parity(boundaries.left, boundaries.right);
        const BoundaryCondition y_parity = get_parity(boundaries.top, boundaries.bottom);

        const std::vector<ModeTuple> filtered_modes = get_filtered_mode_list(x_parity, y_parity);

        if (mode_number >= filtered_modes.size()) {
            throw std::out_of_range(
                "Mode number " + std::to_string(mode_number) +
                " exceeds available modes. Max allowed: " +
                (filtered_modes.empty() ? std::string("none") : std::to_string(filtered_modes.size() - 1))
            );
        }

        const auto& [azimuthal, radial, sub_label] = filtered_modes[mode_number];
        return "LP" + std::to_string(azimuthal) + std::to_string(radial) + sub_label;
    }

private:
    using ModeTuple = std::tuple<int, int, std::string>;

    struct ModeEntry {
        int azimuthal;
        int radial;
        std::string sub_label;
        BoundaryCondition x_parity;
        BoundaryCondition y_parity;
    };

private:
    static BoundaryCondition get_parity(
        const BoundaryCondition boundary_1,
        const BoundaryCondition boundary_2
    ) {
        if (boundary_1 == BoundaryCondition::Symmetric || boundary_2 == BoundaryCondition::Symmetric)
            return BoundaryCondition::Symmetric;

        if (boundary_1 == BoundaryCondition::AntiSymmetric || boundary_2 == BoundaryCondition::AntiSymmetric)
            return BoundaryCondition::AntiSymmetric;

        return BoundaryCondition::Zero;
    }

    static const std::vector<std::pair<int, int>>& get_azimuthal_radial_pairs()
    {
        static const std::vector<std::pair<int, int>> pairs = {
            {0, 1}, {1, 1}, {2, 1}, {0, 2}, {3, 1}, {1, 2}, {4, 1},
            {2, 2}, {0, 3}, {5, 1}, {3, 2}, {1, 3}, {6, 1},
        };
        return pairs;
    }

    static const std::vector<ModeEntry>& get_mode_catalog()
    {
        static const std::vector<ModeEntry> catalog = []() {
            std::vector<ModeEntry> entries;
            entries.reserve(64);

            for (const auto& [azimuthal, radial] : get_azimuthal_radial_pairs())
            {
                std::vector<std::array<BoundaryCondition, 2>> parities;
                std::vector<std::string> sublabels;

                if (azimuthal == 0) {
                    parities = { {BoundaryCondition::Symmetric, BoundaryCondition::Symmetric} };
                    sublabels = { "" };
                }
                else if (azimuthal % 2 == 1) {
                    parities = {
                        {BoundaryCondition::AntiSymmetric, BoundaryCondition::Symmetric},
                        {BoundaryCondition::Symmetric, BoundaryCondition::AntiSymmetric},
                    };
                    sublabels = { "_a", "_b" };
                }
                else {
                    parities = {
                        {BoundaryCondition::Symmetric, BoundaryCondition::Symmetric},
                        {BoundaryCondition::AntiSymmetric, BoundaryCondition::AntiSymmetric},
                    };
                    sublabels = { "_a", "_b" };
                }

                for (size_t i = 0; i < parities.size(); ++i)
                {
                    ModeEntry entry;
                    entry.azimuthal = azimuthal;
                    entry.radial = radial;
                    entry.sub_label = sublabels[i];
                    entry.x_parity = parities[i][0];
                    entry.y_parity = parities[i][1];
                    entries.push_back(std::move(entry));
                }
            }

            return entries;
        }();

        return catalog;
    }

    static std::vector<ModeTuple> get_filtered_mode_list(
        const BoundaryCondition x_parity,
        const BoundaryCondition y_parity
    ) {
        const auto& catalog = get_mode_catalog();

        std::vector<ModeTuple> filtered_modes;
        filtered_modes.reserve(catalog.size());

        for (const ModeEntry& entry : catalog)
        {
            const bool x_ok = (x_parity == BoundaryCondition::Zero) || (x_parity == entry.x_parity);
            const bool y_ok = (y_parity == BoundaryCondition::Zero) || (y_parity == entry.y_parity);

            if (x_ok && y_ok)
                filtered_modes.emplace_back(entry.azimuthal, entry.radial, entry.sub_label);
        }

        return filtered_modes;
    }
};
