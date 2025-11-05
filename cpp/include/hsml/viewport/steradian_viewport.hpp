#pragma once

#include "../core/spherical_types.hpp"
#include <array>
#include <functional>
#include <memory>

namespace hsml::viewport {

template<typename T = double>
class SteradianViewport {
private:
    static constexpr size_t CORNER_COUNT = 4;
    
    std::array<sdt::SphericalCoord<T>, CORNER_COUNT> corners_;
    sdt::SphericalCoord<T> user_position_;
    T interpolation_quality_ = T(1.0);
    T bubble_radius_ = T(1000000.0);  // Default 1000km radius (in meters)
    
    // Interpolation cache for performance
    struct InterpolationCache {
        bool valid = false;
        std::array<T, CORNER_COUNT> weights{};
        sdt::SphericalCoord<T> result{};
    };
    mutable InterpolationCache cache_;
    
public:
    explicit SteradianViewport(T bubble_radius = T(1000000.0)) noexcept
        : bubble_radius_(bubble_radius) {
        // Use the bubble radius directly for corner positions
        T corner_radius = bubble_radius;
        
        // Initialize corners to default steradian coverage
        corners_[0] = {corner_radius, T(90), T(0)};    // Front-left
        corners_[1] = {corner_radius, T(90), T(90)};   // Front-right
        corners_[2] = {corner_radius, T(270), T(0)};   // Back-left
        corners_[3] = {corner_radius, T(270), T(90)};  // Back-right
        
        // Initialize default user position at center with typical screen viewing distance (400mm)
        user_position_ = {T(0.4), T(180), T(180)};  // r=0.4m, theta=180°, phi=180°
    }
    
    // Set the four corner positions defining the viewport
    void set_corners(const std::array<sdt::SphericalCoord<T>, 4>& corners) noexcept {
        corners_ = corners;
        invalidate_cache();
    }
    
    // Update user position for viewport following
    void set_user_position(const sdt::SphericalCoord<T>& position) noexcept {
        if (user_position_.spherical_distance(position) > T(0.001)) {
            user_position_ = position;
            invalidate_cache();
        }
    }
    
    // Set interpolation quality (0 = linear, 1 = perfect spherical)
    void set_interpolation_quality(T quality) noexcept {
        interpolation_quality_ = std::clamp(quality, T(0), T(1));
        invalidate_cache();
    }
    
    // Calculate current viewport position using four-corner interpolation
    [[nodiscard]] sdt::SphericalCoord<T> calculate_viewport_position() const noexcept {
        if (cache_.valid) {
            return cache_.result;
        }
        
        auto weights = calculate_corner_weights();
        auto result = interpolate_spherical_coordinates(weights);
        
        // Cache the result
        cache_.weights = weights;
        cache_.result = result;
        cache_.valid = true;
        
        return result;
    }
    
    // Calculate the volumetric bubble enclosed by the viewport
    [[nodiscard]] T calculate_bubble_volume() const noexcept {
        // Using spherical simplex (tetrahedron) volume calculation
        auto center = calculate_viewport_position();
        T total_volume = T(0);
        
        // Calculate volume of each tetrahedral segment
        for (size_t i = 0; i < CORNER_COUNT; ++i) {
            size_t j = (i + 1) % CORNER_COUNT;
            total_volume += tetrahedron_volume(center, corners_[i], corners_[j], user_position_);
        }
        
        return total_volume;
    }
    
    // Check if a point is within the viewport bubble
    [[nodiscard]] bool contains_point(const sdt::SphericalCoord<T>& point) const noexcept {
        auto center = calculate_viewport_position();
        T distance = center.spherical_distance(point);
        return distance <= bubble_radius_;
    }
    
    // Get solid angle coverage of the viewport
    [[nodiscard]] T solid_angle_coverage() const noexcept {
        // Calculate solid angle subtended by the four corners
        T total_angle = T(0);
        
        for (size_t i = 0; i < CORNER_COUNT; ++i) {
            size_t j = (i + 1) % CORNER_COUNT;
            size_t k = (i + 2) % CORNER_COUNT;
            
            // Spherical triangle angle using spherical trigonometry
            total_angle += spherical_triangle_angle(corners_[i], corners_[j], corners_[k]);
        }
        
        return total_angle;
    }
    
    // Update viewport for smooth following
    void update_viewport_following(T delta_time) noexcept {
        // Smooth interpolation towards user position
        auto target = calculate_optimal_viewport_position();
        auto current = calculate_viewport_position();
        
        T lerp_factor = T(1) - std::exp(-delta_time * T(5)); // Smooth exponential decay
        
        // Interpolate each corner
        for (auto& corner : corners_) {
            corner = spherical_lerp(corner, target, lerp_factor);
        }
        
        invalidate_cache();
    }
    
    // Getters
    [[nodiscard]] const auto& corners() const noexcept { return corners_; }
    [[nodiscard]] const auto& user_position() const noexcept { return user_position_; }
    [[nodiscard]] T interpolation_quality() const noexcept { return interpolation_quality_; }
    [[nodiscard]] T bubble_radius() const noexcept { return bubble_radius_; }
    
private:
    void invalidate_cache() const noexcept {
        cache_.valid = false;
    }
    
    // Calculate weights for each corner based on user position
    [[nodiscard]] std::array<T, CORNER_COUNT> calculate_corner_weights() const noexcept {
        std::array<T, CORNER_COUNT> weights{};
        T total_weight = T(0);
        
        for (size_t i = 0; i < CORNER_COUNT; ++i) {
            T distance = user_position_.spherical_distance(corners_[i]);
            
            // Avoid division by zero and handle large distances
            if (distance < T(0.001)) {
                weights.fill(T(0));
                weights[i] = T(1);
                return weights;
            }
            
            // Use clamped distance to avoid numerical instability with very large values
            T clamped_distance = std::min(distance, T(1000.0));
            
            // Inverse distance weighting with quality adjustment
            T weight = T(1) / (clamped_distance * clamped_distance);
            if (interpolation_quality_ < T(1)) {
                // Blend with linear interpolation
                T linear_weight = T(1) / clamped_distance;
                weight = interpolation_quality_ * weight + (T(1) - interpolation_quality_) * linear_weight;
            }
            
            weights[i] = weight;
            total_weight += weight;
        }
        
        // Normalize weights
        if (total_weight > T(0)) {
            for (auto& weight : weights) {
                weight /= total_weight;
            }
        }
        
        return weights;
    }
    
    // Spherical interpolation of coordinates using weights
    [[nodiscard]] sdt::SphericalCoord<T> interpolate_spherical_coordinates(
        const std::array<T, CORNER_COUNT>& weights) const noexcept {
        
        sdt::SphericalCoord<T> result{T(0), T(0), T(0)};
        
        // Weighted average in spherical space
        for (size_t i = 0; i < CORNER_COUNT; ++i) {
            result.r += weights[i] * corners_[i].r;
            result.theta += weights[i] * corners_[i].theta;
            result.phi += weights[i] * corners_[i].phi;
        }
        
        // Ensure angles remain in safe range
        result.theta = sdt::SphericalCoord<T>::safe_angle(result.theta);
        result.phi = sdt::SphericalCoord<T>::safe_angle(result.phi);
        
        return result;
    }
    
    // Calculate optimal viewport position for following user
    [[nodiscard]] sdt::SphericalCoord<T> calculate_optimal_viewport_position() const noexcept {
        // Position viewport ahead of user movement
        // This is a simplified version - could be enhanced with velocity prediction
        return user_position_;
    }
    
    // Spherical linear interpolation
    [[nodiscard]] sdt::SphericalCoord<T> spherical_lerp(
        const sdt::SphericalCoord<T>& a,
        const sdt::SphericalCoord<T>& b,
        T t) const noexcept {
        
        return {
            a.r + t * (b.r - a.r),
            sdt::SphericalCoord<T>::safe_angle(a.theta + t * (b.theta - a.theta)),
            sdt::SphericalCoord<T>::safe_angle(a.phi + t * (b.phi - a.phi))
        };
    }
    
    // Calculate volume of tetrahedron in spherical space
    [[nodiscard]] T tetrahedron_volume(
        const sdt::SphericalCoord<T>& a,
        const sdt::SphericalCoord<T>& b,
        const sdt::SphericalCoord<T>& c,
        const sdt::SphericalCoord<T>& d) const noexcept {
        
        // Using spherical determinant formula
        // This is a simplified calculation - proper spherical geometry would be more complex
        T avg_r = (a.r + b.r + c.r + d.r) / T(4);
        return avg_r * avg_r * avg_r / T(6); // Approximate
    }
    
    // Calculate angle of spherical triangle
    [[nodiscard]] T spherical_triangle_angle(
        const sdt::SphericalCoord<T>& a,
        const sdt::SphericalCoord<T>& b,
        const sdt::SphericalCoord<T>& c) const noexcept {
        
        // Convert to unit sphere for calculation
        T theta_a = a.theta * M_PI / T(180);
        T theta_b = b.theta * M_PI / T(180);
        T theta_c = c.theta * M_PI / T(180);
        T phi_a = a.phi * M_PI / T(180);
        T phi_b = b.phi * M_PI / T(180);
        T phi_c = c.phi * M_PI / T(180);
        
        // Using spherical law of cosines
        T cos_a = std::cos(theta_a);
        T cos_b = std::cos(theta_b);
        T cos_c = std::cos(theta_c);
        
        T sin_a = std::sin(theta_a);
        T sin_b = std::sin(theta_b);
        T sin_c = std::sin(theta_c);
        
        // Approximate calculation for now
        return std::acos(cos_a * cos_b + sin_a * sin_b * std::cos(phi_a - phi_b));
    }
};

// Type aliases for common use cases
using SteradianViewportf = SteradianViewport<float>;
using SteradianViewportd = SteradianViewport<double>;

} // namespace hsml::viewport