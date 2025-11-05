#pragma once

#include "../core/spherical_types.hpp"
#include "sdt_constants.hpp"
#include <functional>

namespace hsml::physics {

/**
 * SDT Field - Spatial displacement field in pure spherical space
 * NO FORCES - Only spatial displacement
 */
template<typename T = double>
class SDTField {
public:
    using FieldFunction = std::function<sdt::SphericalCoord<T>(const sdt::SphericalCoord<T>&)>;

private:
    sdt::SphericalCoord<T> center_;
    T strength_;
    T range_;
    FieldFunction field_fn_;

public:
    explicit SDTField(
        const sdt::SphericalCoord<T>& center = {},
        T strength = T(1),
        T range = T(10)
    ) noexcept
        : center_(center)
        , strength_(strength)
        , range_(range)
    {
        // Default radial field
        field_fn_ = [this](const sdt::SphericalCoord<T>& pos) -> sdt::SphericalCoord<T> {
            T distance = center_.spherical_distance(pos);
            if (distance < SDTConstants::EPSILON_DISTANCE) {
                return {T(0), T(0), T(0)};
            }

            // Radial displacement
            T radial_strength = strength_ * std::exp(-distance / range_);

            // Direction in spherical space
            T delta_r = pos.r - center_.r;
            T delta_theta = pos.theta - center_.theta;
            T delta_phi = pos.phi - center_.phi;

            T magnitude = std::sqrt(delta_r * delta_r + delta_theta * delta_theta + delta_phi * delta_phi);
            if (magnitude < T(0.001)) {
                return {T(0), T(0), T(0)};
            }

            return {
                radial_strength * delta_r / magnitude,
                radial_strength * delta_theta / magnitude,
                radial_strength * delta_phi / magnitude
            };
        };
    }

    // Calculate displacement at a point
    [[nodiscard]] sdt::SphericalCoord<T> calculate_displacement_at(
        const sdt::SphericalCoord<T>& position) const {
        return field_fn_(position);
    }

    // Update field dynamics
    void update(T delta_time) {
        // Fields can evolve over time
        // Default implementation: static field
        (void)delta_time;
    }

    // Setters
    void set_center(const sdt::SphericalCoord<T>& center) noexcept { center_ = center; }
    void set_strength(T strength) noexcept { strength_ = strength; }
    void set_range(T range) noexcept { range_ = std::max(T(0.001), range); }

    void set_field_function(FieldFunction fn) noexcept {
        field_fn_ = std::move(fn);
    }

    // Getters
    [[nodiscard]] const sdt::SphericalCoord<T>& center() const noexcept { return center_; }
    [[nodiscard]] T strength() const noexcept { return strength_; }
    [[nodiscard]] T range() const noexcept { return range_; }
};

// Type aliases
using SDTFieldf = SDTField<float>;
using SDTFieldd = SDTField<double>;

} // namespace hsml::physics
