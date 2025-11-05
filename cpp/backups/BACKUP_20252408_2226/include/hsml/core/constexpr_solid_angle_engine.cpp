#pragma once

#include <array>
#include <concepts>
#include <type_traits>
#include <cmath>
#include <numbers>
#include "spherical_coords.h"

namespace hsml {
namespace core {

// Mathematical constants
template<typename T>
constexpr T π = std::numbers::pi_v<T>;

template<typename T>
constexpr T TWO_π = T{2} * π<T>;

// Steradian value structure
template<typename T>
struct steradian_value {
    T theta;
    T solid_angle;
    T differential;
    
    constexpr bool operator==(const steradian_value& other) const noexcept {
        return theta == other.theta && 
               solid_angle == other.solid_angle && 
               differential == other.differential;
    }
};

// Pixel coordinate structure
struct pixel_coordinate {
    int x, y;
    constexpr bool operator==(const pixel_coordinate& other) const noexcept {
        return x == other.x && y == other.y;
    }
};

// Display geometry structure
struct display_geometry {
    int width, height;
    double field_of_view;
    double aspect_ratio;
    
    constexpr bool operator==(const display_geometry& other) const noexcept {
        return width == other.width && 
               height == other.height && 
               field_of_view == other.field_of_view && 
               aspect_ratio == other.aspect_ratio;
    }
};

// Steradian calculations with compile-time optimization
template<typename T>
class constexpr_solid_angle_engine {
    static_assert(std::is_floating_point_v<T>);
    
public:
    // Precompute lookup table at compile time
    template<size_t N>
    static consteval auto precompute_lookup_table() noexcept {
        std::array<steradian_value<T>, N> table{};
        
        for (size_t i = 0; i < N; ++i) {
            T theta = static_cast<T>(i) * (π<T> / static_cast<T>(N - 1));
            table[i] = steradian_value<T>{
                .theta = theta,
                .solid_angle = T{2} * π<T> * (T{1} - std::cos(theta)),
                .differential = T{2} * π<T> * std::sin(theta)
            };
        }
        
        return table;
    }
    
    // Static lookup table - computed at compile time
    static constexpr auto lookup_table = precompute_lookup_table<1024>();
    
    // Calculate pixel steradian with compile-time/runtime dispatch
    [[nodiscard]] static constexpr auto 
    calculate_pixel_steradian(pixel_coordinate pixel, display_geometry geometry) noexcept -> T {
        if (std::is_constant_evaluated()) {
            return compile_time_calculation(pixel, geometry);
        } else {
            return runtime_interpolation(pixel, geometry);
        }
    }
    
    // Calculate solid angle for a spherical cap
    [[nodiscard]] static constexpr auto 
    calculate_spherical_cap_solid_angle(T cap_height, T sphere_radius) noexcept -> T {
        if (cap_height <= T{0}) return T{0};
        if (cap_height >= T{2} * sphere_radius) return T{4} * π<T> * sphere_radius * sphere_radius;
        
        return T{2} * π<T> * sphere_radius * cap_height;
    }
    
    // Calculate solid angle for a cone
    [[nodiscard]] static constexpr auto 
    calculate_cone_solid_angle(T half_angle_radians) noexcept -> T {
        return T{2} * π<T> * (T{1} - std::cos(half_angle_radians));
    }
    
    // Calculate solid angle for a spherical rectangle
    [[nodiscard]] static constexpr auto 
    calculate_spherical_rectangle_solid_angle(T theta_min, T theta_max, T phi_min, T phi_max) noexcept -> T {
        T delta_phi = phi_max - phi_min;
        T cos_theta_min = std::cos(theta_min);
        T cos_theta_max = std::cos(theta_max);
        
        return delta_phi * (cos_theta_min - cos_theta_max);
    }
    
    // Calculate solid angle between two directions
    [[nodiscard]] static constexpr auto 
    calculate_solid_angle_between_directions(const SphericalCoords& dir1, const SphericalCoords& dir2) noexcept -> T {
        T cos_angle = std::cos(dir1.theta()) * std::cos(dir2.theta()) + 
                     std::sin(dir1.theta()) * std::sin(dir2.theta()) * 
                     std::cos(dir1.phi() - dir2.phi());
        
        cos_angle = std::clamp(cos_angle, T{-1}, T{1});
        return T{2} * π<T> * (T{1} - cos_angle);
    }
    
    // Monte Carlo solid angle estimation
    template<typename F>
    [[nodiscard]] static constexpr auto 
    monte_carlo_solid_angle_estimation(F&& is_inside_region, size_t num_samples) noexcept -> T {
        if consteval {
            // Compile-time estimation with limited samples
            size_t inside_count = 0;
            for (size_t i = 0; i < std::min(num_samples, size_t{100}); ++i) {
                T theta = static_cast<T>(i) * π<T> / static_cast<T>(std::min(num_samples, size_t{100}));
                T phi = static_cast<T>(i) * TWO_π<T> / static_cast<T>(std::min(num_samples, size_t{100}));
                
                SphericalCoords sample(1.0, theta, phi);
                if (is_inside_region(sample)) {
                    ++inside_count;
                }
            }
            
            return T{4} * π<T> * static_cast<T>(inside_count) / static_cast<T>(std::min(num_samples, size_t{100}));
        } else {
            // Runtime estimation with full samples
            size_t inside_count = 0;
            for (size_t i = 0; i < num_samples; ++i) {
                T theta = static_cast<T>(i) * π<T> / static_cast<T>(num_samples);
                T phi = static_cast<T>(i) * TWO_π<T> / static_cast<T>(num_samples);
                
                SphericalCoords sample(1.0, theta, phi);
                if (is_inside_region(sample)) {
                    ++inside_count;
                }
            }
            
            return T{4} * π<T> * static_cast<T>(inside_count) / static_cast<T>(num_samples);
        }
    }
    
    // Adaptive grid for solid angle calculations
    struct adaptive_grid_cell {
        T theta_min, theta_max;
        T phi_min, phi_max;
        T solid_angle;
        size_t subdivision_level;
        bool needs_refinement;
        
        constexpr adaptive_grid_cell(T t_min, T t_max, T p_min, T p_max, size_t level = 0) noexcept
            : theta_min(t_min), theta_max(t_max), phi_min(p_min), phi_max(p_max)
            , solid_angle(calculate_spherical_rectangle_solid_angle(t_min, t_max, p_min, p_max))
            , subdivision_level(level)
            , needs_refinement(false) {}
    };
    
    // Create adaptive grid
    [[nodiscard]] static constexpr auto 
    create_adaptive_grid(size_t initial_theta_divisions, size_t initial_phi_divisions) noexcept 
        -> std::vector<adaptive_grid_cell> {
        
        std::vector<adaptive_grid_cell> grid;
        grid.reserve(initial_theta_divisions * initial_phi_divisions);
        
        T theta_step = π<T> / static_cast<T>(initial_theta_divisions);
        T phi_step = TWO_π<T> / static_cast<T>(initial_phi_divisions);
        
        for (size_t i = 0; i < initial_theta_divisions; ++i) {
            T theta_min = static_cast<T>(i) * theta_step;
            T theta_max = static_cast<T>(i + 1) * theta_step;
            
            for (size_t j = 0; j < initial_phi_divisions; ++j) {
                T phi_min = -π<T> + static_cast<T>(j) * phi_step;
                T phi_max = -π<T> + static_cast<T>(j + 1) * phi_step;
                
                grid.emplace_back(theta_min, theta_max, phi_min, phi_max);
            }
        }
        
        return grid;
    }
    
private:
    // Compile-time calculation
    [[nodiscard]] static constexpr auto 
    compile_time_calculation(pixel_coordinate pixel, display_geometry geometry) noexcept -> T {
        // Simplified compile-time calculation
        T normalized_x = static_cast<T>(pixel.x) / static_cast<T>(geometry.width);
        T normalized_y = static_cast<T>(pixel.y) / static_cast<T>(geometry.height);
        
        T theta = (normalized_y - T{0.5}) * geometry.field_of_view;
        T phi = (normalized_x - T{0.5}) * geometry.field_of_view * geometry.aspect_ratio;
        
        return calculate_spherical_rectangle_solid_angle(
            theta - T{0.5} * geometry.field_of_view / static_cast<T>(geometry.height),
            theta + T{0.5} * geometry.field_of_view / static_cast<T>(geometry.height),
            phi - T{0.5} * geometry.field_of_view * geometry.aspect_ratio / static_cast<T>(geometry.width),
            phi + T{0.5} * geometry.field_of_view * geometry.aspect_ratio / static_cast<T>(geometry.width)
        );
    }
    
    // Runtime interpolation using lookup table
    [[nodiscard]] static constexpr auto 
    runtime_interpolation(pixel_coordinate pixel, display_geometry geometry) noexcept -> T {
        T normalized_x = static_cast<T>(pixel.x) / static_cast<T>(geometry.width);
        T normalized_y = static_cast<T>(pixel.y) / static_cast<T>(geometry.height);
        
        T theta = (normalized_y - T{0.5}) * geometry.field_of_view;
        T phi = (normalized_x - T{0.5}) * geometry.field_of_view * geometry.aspect_ratio;
        
        // Use lookup table for interpolation
        size_t index = static_cast<size_t>(std::abs(theta) * static_cast<T>(lookup_table.size() - 1) / π<T>);
        index = std::min(index, lookup_table.size() - 1);
        
        return lookup_table[index].solid_angle * 
               geometry.field_of_view * geometry.aspect_ratio / 
               (static_cast<T>(geometry.width) * static_cast<T>(geometry.height));
    }
};

// Compile-time tests
namespace compile_time_tests {
    consteval bool test_lookup_table_creation() {
        constexpr auto table = constexpr_solid_angle_engine<double>::lookup_table;
        return table.size() == 1024 && 
               table[0].theta == 0.0 && 
               table[1023].theta == std::numbers::pi_v<double>;
    }
    
    consteval bool test_spherical_cap_calculation() {
        constexpr auto solid_angle = constexpr_solid_angle_engine<double>::calculate_spherical_cap_solid_angle(1.0, 1.0);
        return solid_angle > 0.0 && solid_angle < 4.0 * std::numbers::pi_v<double>;
    }
    
    consteval bool test_cone_calculation() {
        constexpr auto solid_angle = constexpr_solid_angle_engine<double>::calculate_cone_solid_angle(std::numbers::pi_v<double> / 4.0);
        return solid_angle > 0.0 && solid_angle < 4.0 * std::numbers::pi_v<double>;
    }
    
    consteval bool test_pixel_steradian_calculation() {
        constexpr pixel_coordinate pixel{100, 100};
        constexpr display_geometry geometry{800, 600, 1.0, 4.0/3.0};
        constexpr auto solid_angle = constexpr_solid_angle_engine<double>::calculate_pixel_steradian(pixel, geometry);
        return solid_angle > 0.0;
    }
    
    static_assert(test_lookup_table_creation());
    static_assert(test_spherical_cap_calculation());
    static_assert(test_cone_calculation());
    static_assert(test_pixel_steradian_calculation());
}

} // namespace core
} // namespace hsml
