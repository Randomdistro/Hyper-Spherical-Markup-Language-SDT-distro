#pragma once

#include <array>
#include <cmath>
#include <concepts>
#include <immintrin.h>

namespace hsml::sdt {

// Concept for numeric types that can represent spherical coordinates
template<typename T>
concept SphericalNumeric = std::floating_point<T> || std::integral<T>;

// Pure spherical coordinate - the only truth
template<SphericalNumeric T = double>
struct alignas(32) SphericalCoord {
    T r;      // radius
    T theta;  // 1-361 degrees (zero-division safe)
    T phi;    // 1-361 degrees (zero-division safe)
    T _pad;   // padding for SIMD alignment
    
    constexpr SphericalCoord(T radius = T(1), T th = T(181), T ph = T(181)) noexcept
        : r(radius), theta(safe_angle(th)), phi(safe_angle(ph)), _pad(T(0)) {}
    
    // Ensure angles remain in safe 1-361 range
    static constexpr T safe_angle(T angle) noexcept {
        if (angle < T(1)) return T(1);
        if (angle > T(361)) return T(361);
        return angle;
    }
    
    // Distance in pure spherical space
    [[nodiscard]] T spherical_distance(const SphericalCoord& other) const noexcept {
        // Using law of cosines in spherical space
        T theta1_rad = theta * M_PI / T(180);
        T theta2_rad = other.theta * M_PI / T(180);
        T phi1_rad = phi * M_PI / T(180);
        T phi2_rad = other.phi * M_PI / T(180);
        
        T cos_angle = std::sin(theta1_rad) * std::sin(theta2_rad) * std::cos(phi1_rad - phi2_rad) +
                      std::cos(theta1_rad) * std::cos(theta2_rad);
        
        return std::sqrt(r * r + other.r * other.r - T(2) * r * other.r * cos_angle);
    }
};

// 21-dimensional state vector
template<typename T = double>
struct alignas(64) State21D {
    static constexpr size_t DIMENSIONS = 21;
    
    // Dimensions organized by levels
    enum Level : uint8_t {
        ZERO_POINT = 0,  // Dimensions 0 (1D total)
        LINE = 1,        // Dimensions 1-2 (3D total)
        PLANE = 2,       // Dimensions 3-5 (6D total)
        SPHERE = 3,      // Dimensions 6-9 (10D total)
        FLUX = 4,        // Dimensions 10-14 (15D total)
        ENERGY = 5       // Dimensions 15-20 (21D total)
    };
    
    union {
        std::array<T, DIMENSIONS> dims;
        struct {
            // Level 0: Zero Point (1D)
            T zero_point;
            
            // Level 1: Line (2D)
            T line_x, line_y;
            
            // Level 2: Plane (3D)
            T plane_x, plane_y, plane_z;
            
            // Level 3: Sphere (4D)
            T sphere_r, sphere_theta, sphere_phi, sphere_t;
            
            // Level 4: Flux (5D)
            T flux_density, flux_direction, flux_rotation, flux_resonance, flux_coherence;
            
            // Level 5: Energy (6D)
            T energy_potential, energy_kinetic, energy_binding, 
              energy_resonant, energy_displacement, energy_quantum;
        };
    };
    
    constexpr State21D() noexcept : dims{} {
        // Initialize to unity sphere state
        zero_point = T(1);
        sphere_r = T(1);
        sphere_theta = T(181);
        sphere_phi = T(181);
    }
    
    // Get dimension range for a level
    static constexpr std::pair<size_t, size_t> level_range(Level level) noexcept {
        constexpr std::array<std::pair<size_t, size_t>, 6> ranges = {{
            {0, 1}, {1, 3}, {3, 6}, {6, 10}, {10, 15}, {15, 21}
        }};
        return ranges[static_cast<size_t>(level)];
    }
    
    // SIMD-optimized norm calculation
    [[nodiscard]] T norm() const noexcept {
        if constexpr (std::same_as<T, float>) {
            __m256 sum = _mm256_setzero_ps();
            for (size_t i = 0; i < DIMENSIONS; i += 8) {
                __m256 v = _mm256_load_ps(&dims[i]);
                sum = _mm256_fmadd_ps(v, v, sum);
            }
            float result[8];
            _mm256_store_ps(result, sum);
            T total = 0;
            for (int i = 0; i < 8; ++i) total += result[i];
            return std::sqrt(total);
        } else {
            T sum = 0;
            for (const auto& d : dims) sum += d * d;
            return std::sqrt(sum);
        }
    }
};

// Matter states with spherical properties
enum class MatterState : uint8_t {
    SOLID = 0,
    LIQUID = 1,
    GAS = 2,
    PLASMA = 3,
    QUANTUM = 4,  // Beyond classical states
    VOID = 5      // Absence of matter
};

// Steradian representation for viewport
template<typename T = double>
struct Steradian {
    T solid_angle;  // In steradians (0 to 4π)
    SphericalCoord<T> center;
    
    constexpr Steradian(T angle = T(4) * M_PI, SphericalCoord<T> c = {}) noexcept
        : solid_angle(angle), center(c) {}
    
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
    T interpolation_quality = T(1);  // 0-1, where 1 is perfect spherical interpolation
    
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
template<typename T = double>
using Displacement = SphericalCoord<T>;

} // namespace hsml::sdt