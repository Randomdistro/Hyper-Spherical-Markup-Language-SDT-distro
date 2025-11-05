#pragma once

#include <cmath>
#include <type_traits>
#include <concepts>
#include <numbers>
#include <algorithm>
#include <vector>
#include <array>

namespace hsml::core {

// Unified constexpr spherical coordinates
// Combines the best features from both .h and .hpp versions

template<typename T>
concept arithmetic_type = std::is_arithmetic_v<T>;

template<typename T>
concept floating_point_type = std::is_floating_point_v<T>;

// Modern vector3 for internal calculations (not for rendering!)
template<floating_point_type T>
struct vector3 {
    T x, y, z;

    constexpr vector3(T x_val = T{}, T y_val = T{}, T z_val = T{}) noexcept
        : x(x_val), y(y_val), z(z_val) {}

    // Vector operations
    constexpr vector3 operator+(const vector3& other) const noexcept {
        return vector3(x + other.x, y + other.y, z + other.z);
    }

    constexpr vector3 operator*(T scalar) const noexcept {
        return vector3(x * scalar, y * scalar, z * scalar);
    }

    constexpr T magnitude() const noexcept {
        return std::sqrt(x * x + y * y + z * z);
    }

    constexpr vector3 normalized() const noexcept {
        T mag = magnitude();
        if (mag < std::numeric_limits<T>::epsilon()) {
            return vector3(T{}, T{}, T{});
        }
        return *this * (T{1} / mag);
    }
};

// Main spherical coordinates class
template<floating_point_type T>
class spherical_coords {
private:
    T radius_;
    T theta_;  // polar angle
    T phi_;    // azimuthal angle

public:
    // Constructors
    constexpr spherical_coords(T r = T{}, T theta = T{}, T phi = T{}) noexcept
        : radius_(r), theta_(theta), phi_(phi) {}

    // Getters
    constexpr T radius() const noexcept { return radius_; }
    constexpr T theta() const noexcept { return theta_; }
    constexpr T phi() const noexcept { return phi_; }

    // Setters with validation
    constexpr void set_radius(T r) noexcept { radius_ = std::max(T{}, r); }
    constexpr void set_theta(T t) noexcept { theta_ = std::clamp(t, T{}, std::numbers::pi_v<T>); }
    constexpr void set_phi(T p) noexcept { phi_ = std::fmod(p, std::numbers::pi_v<T> * T{2}); }

    // Conversions (marked for potential removal in pure spherical future)
    constexpr vector3<T> to_cartesian() const noexcept {
        T sin_theta = std::sin(theta_);
        T cos_theta = std::cos(theta_);
        T sin_phi = std::sin(phi_);
        T cos_phi = std::cos(phi_);

        return vector3<T>(
            radius_ * sin_theta * cos_phi,
            radius_ * sin_theta * sin_phi,
            radius_ * cos_theta
        );
    }

    // Pure spherical operations
    constexpr T angular_distance(const spherical_coords& other) const noexcept {
        T cos_angle = std::sin(theta_) * std::sin(other.theta_) *
                     std::cos(phi_ - other.phi_) +
                     std::cos(theta_) * std::cos(other.theta_);

        return std::acos(std::clamp(cos_angle, T{-1}, T{1}));
    }

    constexpr spherical_coords slerp(const spherical_coords& other, T t) const noexcept {
        if (t <= T{}) return *this;
        if (t >= T{1}) return other;

        T angular_sep = angular_distance(other);

        if (std::abs(angular_sep) < std::numeric_limits<T>::epsilon()) {
            return interpolate_linear(other, t);
        }

        T weight1 = std::sin((T{1} - t) * angular_sep) / std::sin(angular_sep);
        T weight2 = std::sin(t * angular_sep) / std::sin(angular_sep);

        T new_radius = radius_ * weight1 + other.radius_ * weight2;
        T new_theta = theta_ * weight1 + other.theta_ * weight2;
        T new_phi = phi_ * weight1 + other.phi_ * weight2;

        return spherical_coords(new_radius, new_theta, new_phi);
    }

    constexpr spherical_coords interpolate_linear(const spherical_coords& other, T t) const noexcept {
        if (t <= T{}) return *this;
        if (t >= T{1}) return other;

        T new_radius = radius_ + t * (other.radius_ - radius_);
        T new_theta = theta_ + t * (other.theta_ - theta_);
        T new_phi = phi_ + t * (other.phi_ - phi_);

        return spherical_coords(new_radius, new_theta, new_phi);
    }

    // Validation
    constexpr bool is_valid() const noexcept {
        return std::isfinite(radius_) && std::isfinite(theta_) && std::isfinite(phi_) &&
               radius_ >= T{} && theta_ >= T{} && theta_ <= std::numbers::pi_v<T>;
    }

    // Comparison
    constexpr bool approximately_equal(const spherical_coords& other,
                                     T epsilon = std::numeric_limits<T>::epsilon() * T{10}) const noexcept {
        return std::abs(radius_ - other.radius_) < epsilon &&
               std::abs(theta_ - other.theta_) < epsilon &&
               std::abs(phi_ - other.phi_) < epsilon;
    }

    constexpr bool operator==(const spherical_coords& other) const noexcept {
        return approximately_equal(other);
    }

    constexpr bool operator!=(const spherical_coords& other) const noexcept {
        return !approximately_equal(other);
    }
};

// Type aliases
using spherical_coords_f32 = spherical_coords<float>;
using spherical_coords_f64 = spherical_coords<double>;

// Utility functions
template<floating_point_type T>
constexpr T solid_angle(T radius) noexcept {
    T cos_angle = std::cos(radius);
    return T{2} * std::numbers::pi_v<T> * (T{1} - cos_angle);
}

template<floating_point_type T>
constexpr std::vector<spherical_coords<T>> generate_sphere_points(size_t count) {
    std::vector<spherical_coords<T>> points;
    points.reserve(count);

    T golden_ratio = (T{1} + std::sqrt(T{5})) / T{2};

    for (size_t i = 0; i < count; ++i) {
        T y = T{1} - (T{2} * i) / (count - 1);
        T radius = std::sqrt(T{1} - y * y);
        T theta = std::acos(y);
        T phi = (T{2} * std::numbers::pi_v<T> * i) / golden_ratio;

        points.emplace_back(T{1}, theta, phi);
    }

    return points;
}

} // namespace hsml::core
