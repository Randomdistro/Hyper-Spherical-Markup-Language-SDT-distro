#pragma once

#include <array>
#include <cmath>
#include <numbers>
#include <type_traits>
#include <vector>
#include <algorithm>
#include <concepts>

namespace hsml::core {

// Unified constexpr solid angle implementation
// Combines the best features from both .h and .hpp versions

template<typename T>
concept arithmetic_type = std::is_arithmetic_v<T>;

template<typename T>
concept floating_point_type = std::is_floating_point_v<T>;

// Modern steradian value structure
template<floating_point_type T>
struct steradian_value {
    T theta;
    T solid_angle;
    T differential;

    constexpr steradian_value(T t = T{}, T sa = T{}, T diff = T{}) noexcept
        : theta(t), solid_angle(sa), differential(diff) {}
};

// Pixel coordinate structure
struct pixel_coordinate {
    size_t x;
    size_t y;

    constexpr pixel_coordinate(size_t x_val = 0, size_t y_val = 0) noexcept
        : x(x_val), y(y_val) {}
};

// Solid angle computation utilities
template<floating_point_type T>
class solid_angle_computer {
public:
    // Compute solid angle of a spherical cap
    static constexpr T spherical_cap(T radius) noexcept {
        T cos_radius = std::cos(radius);
        return T{2} * std::numbers::pi_v<T> * (T{1} - cos_radius);
    }

    // Compute solid angle between two polar angles
    static constexpr T between_angles(T theta1, T theta2) noexcept {
        T cos1 = std::cos(theta1);
        T cos2 = std::cos(theta2);
        return T{2} * std::numbers::pi_v<T> * (cos1 - cos2);
    }

    // Compute differential solid angle element
    static constexpr T differential(T theta, T d_theta, T d_phi) noexcept {
        T sin_theta = std::sin(theta);
        return sin_theta * d_theta * d_phi;
    }

    // Check if a point is within a solid angle region
    static constexpr bool is_within_region(T point_theta, T point_phi,
                                         T center_theta, T center_phi, T radius) noexcept {
        T angular_distance = angular_distance(point_theta, point_phi, center_theta, center_phi);
        return angular_distance <= radius;
    }

    // Calculate angular distance between two points on sphere
    static constexpr T angular_distance(T theta1, T phi1, T theta2, T phi2) noexcept {
        T cos_angle = std::sin(theta1) * std::sin(theta2) * std::cos(phi1 - phi2) +
                     std::cos(theta1) * std::cos(theta2);

        return std::acos(std::clamp(cos_angle, T{-1}, T{1}));
    }
};

// Compile-time solid angle table generation
template<floating_point_type T, size_t Resolution>
class solid_angle_table {
private:
    std::array<steradian_value<T>, Resolution> table_;

public:
    constexpr solid_angle_table() : table_{} {
        // Generate lookup table at compile time
        for (size_t i = 0; i < Resolution; ++i) {
            T theta = (T{i} / T{Resolution - 1}) * std::numbers::pi_v<T>;
            T solid_angle_val = solid_angle_computer<T>::spherical_cap(theta);
            T differential_val = solid_angle_computer<T>::differential(
                theta, std::numbers::pi_v<T> / T{Resolution}, T{2} * std::numbers::pi_v<T} / T{Resolution}
            );

            table_[i] = steradian_value<T>(theta, solid_angle_val, differential_val);
        }
    }

    // Lookup solid angle value
    constexpr T lookup(T theta) const noexcept {
        if (theta <= T{} || theta >= std::numbers::pi_v<T>) {
            return T{};
        }

        size_t index = static_cast<size_t>((theta / std::numbers::pi_v<T>) * (Resolution - 1));
        index = std::clamp(index, size_t{0}, Resolution - 1);

        return table_[index].solid_angle;
    }

    // Get table size
    static constexpr size_t size() noexcept { return Resolution; }

    // Access raw table data
    constexpr const auto& data() const noexcept { return table_; }
};

// Pre-computed tables for common resolutions
using solid_angle_table_f32_64 = solid_angle_table<float, 64>;
using solid_angle_table_f64_128 = solid_angle_table<double, 128>;

// Pixel-based solid angle calculations
template<floating_point_type T>
class pixel_solid_angle_calculator {
private:
    T pixel_width_;
    T pixel_height_;
    T focal_length_;

public:
    constexpr pixel_solid_angle_calculator(T pixel_width, T pixel_height, T focal_length) noexcept
        : pixel_width_(pixel_width), pixel_height_(pixel_height), focal_length_(focal_length) {}

    // Calculate solid angle subtended by a single pixel
    constexpr T pixel_solid_angle(const pixel_coordinate& pixel) const noexcept {
        // Simplified calculation - would need more complex optics for real camera
        T half_width = pixel_width_ / T{2};
        T half_height = pixel_height_ / T{2};

        T angle_x = std::atan(half_width / focal_length_);
        T angle_y = std::atan(half_height / focal_length_);

        return T{4} * angle_x * angle_y; // Approximation
    }

    // Calculate solid angle for a region of pixels
    constexpr T region_solid_angle(const pixel_coordinate& start,
                                 const pixel_coordinate& end) const noexcept {
        T total_solid_angle = T{};

        for (size_t y = start.y; y <= end.y; ++y) {
            for (size_t x = start.x; x <= end.x; ++x) {
                total_solid_angle += pixel_solid_angle(pixel_coordinate(x, y));
            }
        }

        return total_solid_angle;
    }
};

// Utility functions for solid angle operations
template<floating_point_type T>
constexpr std::vector<pixel_coordinate> generate_pixel_grid(size_t width, size_t height) {
    std::vector<pixel_coordinate> pixels;
    pixels.reserve(width * height);

    for (size_t y = 0; y < height; ++y) {
        for (size_t x = 0; x < width; ++x) {
            pixels.emplace_back(x, y);
        }
    }

    return pixels;
}

template<floating_point_type T>
constexpr T total_sphere_solid_angle() noexcept {
    return T{4} * std::numbers::pi_v<T>;
}

template<floating_point_type T>
constexpr T hemisphere_solid_angle() noexcept {
    return T{2} * std::numbers::pi_v<T>;
}

// Validation functions
template<floating_point_type T>
constexpr bool is_valid_solid_angle(T solid_angle) noexcept {
    return std::isfinite(solid_angle) && solid_angle >= T{} &&
           solid_angle <= total_sphere_solid_angle<T>();
}

template<floating_point_type T>
constexpr bool is_valid_steradian(T steradian) noexcept {
    return is_valid_solid_angle(steradian);
}

} // namespace hsml::core
