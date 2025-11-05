#pragma once

#include <array>
#include <cmath>
#include <type_traits>
#include <cstdint>
#include <immintrin.h>
#include <algorithm>
#include <utility>


namespace hsml::sdt {

// Type constraint for numeric types
template<typename T>
struct is_spherical_numeric {
    static constexpr bool value = std::is_floating_point_v<T> || std::is_integral_v<T>;
};

// Pure spherical coordinate - the only truth
template<typename T = double>
struct alignas(32) SphericalCoord {
    static_assert(is_spherical_numeric<T>::value, "T must be a floating point or integral type");
    T r;      // radius (can be any positive value, including 0.001 - small existence!)
    T theta;  // 1-360 degrees (no zero position!)
    T phi;    // 1-360 degrees (no zero position!)
    T _pad;   // padding for SIMD alignment

    constexpr SphericalCoord(T radius = T(1), T th = T(181), T ph = T(181)) noexcept
        : r(radius > T(0) ? radius : T(1e-12)), theta(safe_angle(th)), phi(safe_angle(ph)), _pad(T(0)) {}
    
    // Ensure angles remain in safe 1-360 range (NO ZERO!)
    // The circle has positions 1-360, where 361 wraps to 1
    static constexpr T safe_angle(T angle) noexcept {
        // Normalize to 1-360 range
        // If angle is 361 or higher, wrap down
        while (angle > T(360)) angle -= T(360);
        // If angle is 0 or lower, wrap up
        while (angle < T(1)) angle += T(360);

        return angle;
    }
    
    // Distance in pure spherical space
    [[nodiscard]] T spherical_distance(const SphericalCoord& other) const noexcept {
        // Convert 1-360 degrees to radians
        // Position 1 = 1 degree, position 360 = 360 degrees
        T theta1_rad = theta * M_PI / T(180);
        T theta2_rad = other.theta * M_PI / T(180);
        T phi1_rad = phi * M_PI / T(180);
        T phi2_rad = other.phi * M_PI / T(180);

        // Spherical law of cosines for great circle distance
        T cos_angle = std::sin(theta1_rad) * std::sin(theta2_rad) * std::cos(phi1_rad - phi2_rad) +
                      std::cos(theta1_rad) * std::cos(theta2_rad);

        // 3D Euclidean distance between points on sphere
        T dist_sq = r * r + other.r * other.r - T(2) * r * other.r * cos_angle;
        // Distance can be tiny (0.0001) but not negative - that's physically impossible!
        return std::sqrt(std::max(dist_sq, T(0)));
    }
};

// 21-dimensional state vector
template<typename T = double>
struct alignas(64) State21D {
    static constexpr size_t DIMENSIONS = 21;
    
    // Dimensions organized by levels (NO ZERO - starting at 1!)
    enum Level : uint8_t {
        UNITY_POINT = 1,  // Level 1: Unity Point (dimension 1)
        LINE = 2,         // Level 2: Line (dimensions 2-3)
        PLANE = 3,        // Level 3: Plane (dimensions 4-6)
        SPHERE = 4,       // Level 4: Sphere (dimensions 7-10)
        FLUX = 5,         // Level 5: Flux (dimensions 11-15)
        ENERGY = 6        // Level 6: Energy (dimensions 16-21)
    };
    
    union {
        std::array<T, DIMENSIONS> dims;
        struct {
            // Level 1: Unity Point (1D) - The singular point (no zero, just unity)
            T unity_point;

            // Level 2: Line (2D)
            T line_x, line_y;

            // Level 3: Plane (3D)
            T plane_x, plane_y, plane_z;

            // Level 4: Sphere (4D)
            T sphere_r, sphere_theta, sphere_phi, sphere_t;

            // Level 5: Flux (5D)
            T flux_density, flux_direction, flux_rotation, flux_resonance, flux_coherence;

            // Level 6: Energy (6D)
            T energy_potential, energy_kinetic, energy_binding,
              energy_resonant, energy_displacement, energy_quantum;
        };
    };
    
    constexpr State21D() noexcept : dims{} {
        // Initialize to unity sphere state
        // Note: Small values (0.001) are valid existence at tiny scales!
        unity_point = T(1);
        sphere_r = T(1);
        sphere_theta = T(181);
        sphere_phi = T(181);
    }
    
    // Get dimension range for a level (1-indexed, no zero!)
    static constexpr std::pair<size_t, size_t> level_range(Level level) noexcept {
        constexpr std::array<std::pair<size_t, size_t>, 7> ranges = {{
            {1, 1},    // Padding (unused - levels start at 1)
            {1, 1},    // UNITY_POINT: dimension 1
            {2, 3},    // LINE: dimensions 2-3
            {4, 6},    // PLANE: dimensions 4-6
            {7, 10},   // SPHERE: dimensions 7-10
            {11, 15},  // FLUX: dimensions 11-15
            {16, 21}   // ENERGY: dimensions 16-21
        }};
        return ranges[static_cast<size_t>(level)];
    }
    
    // SIMD-optimized norm calculation
    [[nodiscard]] T norm() const noexcept {
        if constexpr (std::is_same_v<T, float>) {
            // Manual accumulation for now (SIMD optimization later)
            T sum_of_squares = T(0);
            for (size_t i = 0; i < DIMENSIONS; ++i) {
                sum_of_squares += dims[i] * dims[i];
            }
            return std::sqrt(sum_of_squares);
        } else {
            T sum_of_squares = T(0);
            for (const auto& d : dims) {
                sum_of_squares += d * d;
            }
            return std::sqrt(sum_of_squares);
        }
    }
};

// Matter states with spherical properties (NO ZERO - starting at 1!)
enum class MatterState : uint8_t {
    SOLID = 1,    // Coherent matter
    LIQUID = 2,   // Flowing matter
    GAS = 3,      // Dispersed matter
    PLASMA = 4,   // Ionized matter
    BARYONIC = 5, // Particulate matter
    QUANTUM = 6,  // Beyond classical states
    VORTEX = 7,   // Effect of matter
    VOID = 8      // Absence of matter
};

// Steradian representation for viewport
template<typename T = double>
struct Steradian {
    T solid_angle;  // In steradians (can be tiny like 0.0001, but not negative!)
    SphericalCoord<T> center;

    constexpr Steradian(T angle = T(4) * M_PI, SphericalCoord<T> c = {}) noexcept
        : solid_angle(std::max(angle, T(0))), center(c) {}
    
    // Check if a point is within this steradian
    [[nodiscard]] bool contains(const SphericalCoord<T>& point) const noexcept {
        T distance = center.spherical_distance(point);
        T max_angle = std::sqrt(solid_angle / M_PI);
        return distance <= center.r * max_angle;
    }
};

// Four-corner steradian viewport
template<typename T = double>
struct SteradianViewport {
    std::array<SphericalCoord<T>, 4> corners;
    T interpolation_quality = T(1);  // 1-361, where 1 is perfect spherical interpolation
    
    // Calculate the bubble volume enclosed by the four corners
    [[nodiscard]] T bubble_volume() const noexcept {
        // Using spherical tetrahedron volume formula
        // V = (1/3) * r³ * Ω where Ω is solid angle
        T avg_r = (corners[0].r + corners[1].r + corners[2].r + corners[3].r) / T(4);
        T solid_angle = calculate_solid_angle();
        return (T(1)/T(3)) * avg_r * avg_r * avg_r * solid_angle;
    }
    
private:
    [[nodiscard]] T calculate_solid_angle() const noexcept {
        // Girard's theorem for spherical excess
        // TODO: Implement proper spherical polygon area calculation
        return T(4) * M_PI;  // Full sphere for now
    }
};

// Resonance frequency in Hz
using Resonance = double;

// Spatial displacement vector (not a force!)
template<typename T>
struct Displacement : SphericalCoord<T> {
    using SphericalCoord<T>::SphericalCoord;  // Inherit constructors
};

} // namespace hsml::sdt