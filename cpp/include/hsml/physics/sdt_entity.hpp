#pragma once

#include "../core/spherical_types.hpp"
#include "sdt_constants.hpp"
#include <string>
#include <memory>

namespace hsml::physics {

/**
 * SDT Entity - A physical object in spherical space
 * PURE SPHERICAL - No Cartesian coordinates ever
 */
template<typename T = double>
class SDTEntity {
private:
    sdt::SphericalCoord<T> position_;
    sdt::State21D<T> state_;
    sdt::MatterState matter_state_;
    T radius_;
    T mass_;
    T density_;
    std::string id_;

public:
    explicit SDTEntity(
        const sdt::SphericalCoord<T>& pos = {},
        T radius = T(1),
        T mass = T(1),
        sdt::MatterState matter = sdt::MatterState::SOLID,
        std::string id = ""
    ) noexcept
        : position_(pos)
        , state_{}
        , matter_state_(matter)
        , radius_(radius)
        , mass_(mass)
        , density_(mass / (T(4)/T(3) * M_PI * radius * radius * radius))
        , id_(std::move(id))
    {}

    // Position and movement
    [[nodiscard]] const sdt::SphericalCoord<T>& position() const noexcept { return position_; }
    [[nodiscard]] sdt::SphericalCoord<T>& position() noexcept { return position_; }

    void set_position(const sdt::SphericalCoord<T>& pos) noexcept {
        position_ = pos;
    }

    void displace(const sdt::SphericalCoord<T>& displacement) noexcept {
        position_.r += displacement.r;
        position_.theta = sdt::SphericalCoord<T>::safe_angle(position_.theta + displacement.theta);
        position_.phi = sdt::SphericalCoord<T>::safe_angle(position_.phi + displacement.phi);
    }

    // 21D State
    [[nodiscard]] const sdt::State21D<T>& state() const noexcept { return state_; }
    [[nodiscard]] sdt::State21D<T>& state() noexcept { return state_; }

    void set_state(const sdt::State21D<T>& s) noexcept { state_ = s; }

    // Matter state
    [[nodiscard]] sdt::MatterState matter_state() const noexcept { return matter_state_; }
    void set_matter_state(sdt::MatterState state) noexcept { matter_state_ = state; }

    // Physical properties
    [[nodiscard]] T radius() const noexcept { return radius_; }
    [[nodiscard]] T mass() const noexcept { return mass_; }
    [[nodiscard]] T density() const noexcept { return density_; }

    void set_radius(T r) noexcept {
        radius_ = std::max(T(0.001), r);
        update_density();
    }

    void set_mass(T m) noexcept {
        mass_ = std::max(T(0.001), m);
        update_density();
    }

    // Identification
    [[nodiscard]] const std::string& id() const noexcept { return id_; }
    void set_id(std::string new_id) noexcept { id_ = std::move(new_id); }

    // Calculate energy from 21D state
    [[nodiscard]] T total_energy() const noexcept {
        return state_.energy_potential + state_.energy_kinetic + state_.energy_binding +
               state_.energy_resonant + state_.energy_displacement + state_.energy_quantum;
    }

    // Calculate flux magnitude
    [[nodiscard]] T flux_magnitude() const noexcept {
        return state_.flux_density;
    }

private:
    void update_density() noexcept {
        const T volume = (T(4)/T(3)) * M_PI * radius_ * radius_ * radius_;
        density_ = mass_ / std::max(volume, T(1e-10));
    }
};

// Type aliases
using SDTEntityf = SDTEntity<float>;
using SDTEntityd = SDTEntity<double>;

} // namespace hsml::physics
