#pragma once

#include "mathematical_foundation.h"
#include "optimized_vector3.h"
#include <ostream>
#include <string>
#include <unordered_map>
#include <memory>
#include <chrono>

#ifdef HSML_SSE2_AVAILABLE
#include <emmintrin.h>
#endif

namespace hsml {
namespace core {

// ============================================================================
// UNIFIED SPHERICAL COORDINATES - Single, adaptive implementation
// ============================================================================

template<math::PrecisionLevel Level = math::PrecisionLevel::Standard>
class UnifiedSphericalCoords {
private:
    using traits = math::precision_traits<Level>;
    using safe_math = math::SafeMath<Level>;
    using vector_type = OptimizedVector3<Level>;
    
    // Core spherical coordinates with validated ranges
    math::ValidatedValue<double, Level> radius_;
    math::ValidatedValue<double, Level> theta_;  // Polar angle [0, π]
    math::ValidatedValue<double, Level> phi_;    // Azimuthal angle [-π, π]
    
    // Performance optimization: cached Cartesian conversion
    mutable std::unique_ptr<vector_type> cached_cartesian_;
    mutable bool cartesian_cache_valid_ = false;
    mutable std::chrono::steady_clock::time_point cache_timestamp_;
    
    // Thread-local distance calculation cache
    struct DistanceCache {
        std::unordered_map<std::string, double> cache;
        size_t hits = 0;
        size_t misses = 0;
        static constexpr size_t max_size = 10000;
    };
    
    static thread_local DistanceCache distance_cache_;

public:
    // ========================================================================
    // CONSTRUCTION AND INITIALIZATION
    // ========================================================================
    
    constexpr UnifiedSphericalCoords() noexcept
        : radius_(traits::origin_marker)
        , theta_(traits::origin_marker)
        , phi_(traits::origin_marker) {}
    
    UnifiedSphericalCoords(double radius, double theta, double phi)
        : radius_(std::max(0.0, radius))
        , theta_(theta, 0.0, math::constants::pi)
        , phi_(phi, -math::constants::pi, math::constants::pi) {
        normalize_angles();
    }
    
    // Copy constructor with cache invalidation
    UnifiedSphericalCoords(const UnifiedSphericalCoords& other)
        : radius_(other.radius_)
        , theta_(other.theta_)
        , phi_(other.phi_)
        , cartesian_cache_valid_(false) {}
    
    // Assignment operator with cache invalidation
    UnifiedSphericalCoords& operator=(const UnifiedSphericalCoords& other) {
        if (this != &other) {
            radius_ = other.radius_;
            theta_ = other.theta_;
            phi_ = other.phi_;
            invalidate_cache();
        }
        return *this;
    }
    
    // ========================================================================
    // ELEMENT ACCESS AND MODIFICATION
    // ========================================================================
    
    constexpr double radius() const noexcept { return radius_.get(); }
    constexpr double theta() const noexcept { return theta_.get(); }
    constexpr double phi() const noexcept { return phi_.get(); }
    
    void set_radius(double radius) {
        radius_ = std::max(0.0, radius);
        invalidate_cache();
    }
    
    void set_theta(double theta) {
        theta_ = math::ValidatedValue<double, Level>(theta, 0.0, math::constants::pi);
        normalize_angles();
        invalidate_cache();
    }
    
    void set_phi(double phi) {
        phi_ = math::ValidatedValue<double, Level>(phi, -math::constants::pi, math::constants::pi);
        normalize_angles();
        invalidate_cache();
    }
    
    // ========================================================================
    // COORDINATE SYSTEM CONVERSIONS
    // ========================================================================
    
    vector_type to_cartesian() const {
        // Check cache validity (with timestamp for debugging)
        if (cartesian_cache_valid_ && cached_cartesian_) {
            return *cached_cartesian_;
        }
        
        // Compute Cartesian coordinates with trigonometric precision
        const double r = radius_.get();
        const double sin_theta = std::sin(theta_.get());
        const double cos_theta = std::cos(theta_.get());
        const double sin_phi = std::sin(phi_.get());
        const double cos_phi = std::cos(phi_.get());
        
        vector_type result(
            r * sin_theta * cos_phi,  // x
            r * sin_theta * sin_phi,  // y
            r * cos_theta             // z
        );
        
        // Update cache
        cached_cartesian_ = std::make_unique<vector_type>(result);
        cartesian_cache_valid_ = true;
        cache_timestamp_ = std::chrono::steady_clock::now();
        
        return result;
    }
    
    static UnifiedSphericalCoords from_cartesian(const vector_type& cartesian) {
        const double radius = cartesian.magnitude();
        
        // Handle near-zero radius case
        if (safe_math::is_effectively_zero(radius)) {
            return UnifiedSphericalCoords(traits::origin_marker, traits::origin_marker, traits::origin_marker);
        }
        
        // Safe trigonometric calculations
        const double theta = safe_math::safe_acos(cartesian.z() / radius);
        const double phi = std::atan2(cartesian.y(), cartesian.x());
        
        return UnifiedSphericalCoords(radius, theta, phi);
    }
    
    // ========================================================================
    // DISTANCE CALCULATIONS - Unified, adaptive implementation
    // ========================================================================
    
    double distance_to(const UnifiedSphericalCoords& other) const {
        // Check cache first for performance
        const std::string cache_key = generate_cache_key(other);
        auto& cache = distance_cache_;
        
        auto cache_it = cache.cache.find(cache_key);
        if (cache_it != cache.cache.end()) {
            ++cache.hits;
            return cache_it->second;
        }
        
        // Select optimal calculation method based on data characteristics
        double distance;
        if (should_use_simd_calculation(other)) {
            distance = calculate_distance_simd(other);
        } else if (should_use_direct_spherical(other)) {
            distance = calculate_distance_spherical(other);
        } else {
            distance = calculate_distance_cartesian(other);
        }
        
        // Update cache with size management
        if (cache.cache.size() >= cache.max_size) {
            cache.cache.clear(); // Simple eviction strategy
        }
        cache.cache[cache_key] = distance;
        ++cache.misses;
        
        return distance;
    }
    
    double angular_distance(const UnifiedSphericalCoords& other) const {
        const double cos_angle = 
            std::sin(theta_.get()) * std::sin(other.theta_.get()) * std::cos(phi_.get() - other.phi_.get()) +
            std::cos(theta_.get()) * std::cos(other.theta_.get());
        
        return safe_math::safe_acos(std::clamp(cos_angle, -1.0, 1.0));
    }
    
    // ========================================================================
    // INTERPOLATION AND GEOMETRIC OPERATIONS
    // ========================================================================
    
    UnifiedSphericalCoords interpolate(const UnifiedSphericalCoords& other, double t) const {
        t = std::clamp(t, 0.0, 1.0);
        
        if (t <= 0.0) return *this;
        if (t >= 1.0) return other;
        
        // Perform interpolation in Cartesian space for accuracy
        vector_type cart1 = to_cartesian();
        vector_type cart2 = other.to_cartesian();
        vector_type interpolated = vector_type::lerp(cart1, cart2, t);
        
        return from_cartesian(interpolated);
    }
    
    UnifiedSphericalCoords slerp(const UnifiedSphericalCoords& other, double t) const {
        t = std::clamp(t, 0.0, 1.0);
        
        vector_type v1 = to_cartesian().normalized();
        vector_type v2 = other.to_cartesian().normalized();
        vector_type interpolated = vector_type::slerp(v1, v2, t);
        
        // Interpolate radius separately
        double interpolated_radius = radius_.get() * (1.0 - t) + other.radius_.get() * t;
        
        UnifiedSphericalCoords result = from_cartesian(interpolated * interpolated_radius);
        result.set_radius(interpolated_radius);
        
        return result;
    }
    
    UnifiedSphericalCoords rotate_about_axis(const vector_type& axis, double angle) const {
        vector_type cart = to_cartesian();
        vector_type normalized_axis = axis.normalized();
        
        // Rodrigues' rotation formula
        const double cos_angle = std::cos(angle);
        const double sin_angle = std::sin(angle);
        const double one_minus_cos = 1.0 - cos_angle;
        
        const double x = cart.x(), y = cart.y(), z = cart.z();
        const double ax = normalized_axis.x(), ay = normalized_axis.y(), az = normalized_axis.z();
        
        vector_type rotated(
            x * (cos_angle + ax * ax * one_minus_cos) +
            y * (ax * ay * one_minus_cos - az * sin_angle) +
            z * (ax * az * one_minus_cos + ay * sin_angle),
            
            x * (ay * ax * one_minus_cos + az * sin_angle) +
            y * (cos_angle + ay * ay * one_minus_cos) +
            z * (ay * az * one_minus_cos - ax * sin_angle),
            
            x * (az * ax * one_minus_cos - ay * sin_angle) +
            y * (az * ay * one_minus_cos + ax * sin_angle) +
            z * (cos_angle + az * az * one_minus_cos)
        );
        
        return from_cartesian(rotated);
    }
    
    // ========================================================================
    // VALIDATION AND UTILITY METHODS
    // ========================================================================
    
    bool is_valid() const noexcept {
        return std::isfinite(radius_.get()) && std::isfinite(theta_.get()) && std::isfinite(phi_.get()) && 
               radius_.get() >= 0.0 && theta_.get() >= 0.0 && theta_.get() <= math::constants::pi;
    }
    
    bool approximately_equal(const UnifiedSphericalCoords& other) const noexcept {
        return safe_math::are_equal(radius_.get(), other.radius_.get()) &&
               safe_math::are_equal(theta_.get(), other.theta_.get()) &&
               safe_math::are_equal(phi_.get(), other.phi_.get());
    }
    
    double solid_angle_to(const UnifiedSphericalCoords& other) const {
        vector_type v1 = to_cartesian().normalized();
        vector_type v2 = other.to_cartesian().normalized();
        
        double cos_angle = std::clamp(v1.dot(v2), -1.0, 1.0);
        return 2.0 * math::constants::pi * (1.0 - cos_angle);
    }
    
    // ========================================================================
    // CACHE MANAGEMENT AND STATISTICS
    // ========================================================================
    
    static void clear_distance_cache() {
        distance_cache_.cache.clear();
        distance_cache_.hits = 0;
        distance_cache_.misses = 0;
    }
    
    static std::pair<size_t, size_t> get_cache_stats() {
        return {distance_cache_.hits, distance_cache_.misses};
    }
    
    void invalidate_cache() const {
        cartesian_cache_valid_ = false;
        cached_cartesian_.reset();
    }

private:
    // ========================================================================
    // INTERNAL CALCULATION METHODS
    // ========================================================================
    
    void normalize_angles() {
        double phi_val = phi_.get();
        while (phi_val > math::constants::pi) phi_val -= math::constants::two_pi;
        while (phi_val < -math::constants::pi) phi_val += math::constants::two_pi;
        phi_ = phi_val;
        
        double theta_val = std::clamp(theta_.get(), 0.0, math::constants::pi);
        theta_ = theta_val;
    }
    
    std::string generate_cache_key(const UnifiedSphericalCoords& other) const {
        // Use precision-appropriate formatting for cache keys
        constexpr int precision = static_cast<int>(-std::log10(traits::epsilon));
        
        char buffer[256];
        std::snprintf(buffer, sizeof(buffer), 
                     "dist_%.*f_%.*f_%.*f_%.*f_%.*f_%.*f",
                     precision, radius_.get(), precision, theta_.get(), precision, phi_.get(),
                     precision, other.radius_.get(), precision, other.theta_.get(), precision, other.phi_.get());
        return std::string(buffer);
    }
    
    bool should_use_simd_calculation(const UnifiedSphericalCoords& other) const {
#ifdef HSML_SSE2_AVAILABLE
        // Use SIMD for moderate to large radius values where precision loss is acceptable
        return (radius_.get() > 1.0 && other.radius_.get() > 1.0) && 
               (Level >= math::PrecisionLevel::Engineering);
#else
        return false;
#endif
    }
    
    bool should_use_direct_spherical(const UnifiedSphericalCoords& other) const {
        // Use direct spherical calculation when radii are similar and angles are well-conditioned
        const double radius_ratio = radius_.get() / std::max(other.radius_.get(), traits::origin_marker);
        return (radius_ratio > 0.1 && radius_ratio < 10.0) &&
               !safe_math::is_effectively_zero(std::sin(theta_.get())) &&
               !safe_math::is_effectively_zero(std::sin(other.theta_.get()));
    }
    
    double calculate_distance_simd(const UnifiedSphericalCoords& other) const {
#ifdef HSML_SSE2_AVAILABLE
        // SIMD-optimized trigonometric calculations
        const double cos_theta1 = std::cos(theta_.get());
        const double sin_theta1 = std::sin(theta_.get());
        const double cos_theta2 = std::cos(other.theta_.get());
        const double sin_theta2 = std::sin(other.theta_.get());
        const double cos_phi_diff = std::cos(other.phi_.get() - phi_.get());
        
        // Vectorized dot product calculation
        __m128d cos_vals = _mm_set_pd(cos_theta2, cos_theta1);
        __m128d sin_vals = _mm_set_pd(sin_theta2, sin_theta1);
        __m128d cos_prod = _mm_mul_pd(cos_vals, _mm_shuffle_pd(cos_vals, cos_vals, 1));
        __m128d sin_prod = _mm_mul_pd(sin_vals, _mm_shuffle_pd(sin_vals, sin_vals, 1));
        
        double cos_angular = _mm_cvtsd_f64(cos_prod) + _mm_cvtsd_f64(sin_prod) * cos_phi_diff;
        cos_angular = std::clamp(cos_angular, -1.0, 1.0);
        
        // 3D distance using law of cosines
        const double r1 = std::max(radius_.get(), traits::origin_marker);
        const double r2 = std::max(other.radius_.get(), traits::origin_marker);
        
        __m128d r_vals = _mm_set_pd(r2, r1);
        __m128d r_squared = _mm_mul_pd(r_vals, r_vals);
        double r_sum_squared = _mm_cvtsd_f64(r_squared) + 
                              _mm_cvtsd_f64(_mm_unpackhi_pd(r_squared, r_squared));
        
        double distance_squared = r_sum_squared - 2.0 * r1 * r2 * cos_angular;
        return safe_math::safe_sqrt(distance_squared);
#else
        return calculate_distance_spherical(other);
#endif
    }
    
    double calculate_distance_spherical(const UnifiedSphericalCoords& other) const {
        // Direct spherical distance calculation using law of cosines
        const double cos_theta1 = std::cos(theta_.get());
        const double sin_theta1 = std::sin(theta_.get());
        const double cos_theta2 = std::cos(other.theta_.get());
        const double sin_theta2 = std::sin(other.theta_.get());
        const double cos_phi_diff = std::cos(other.phi_.get() - phi_.get());
        
        const double cos_angular = cos_theta1 * cos_theta2 + sin_theta1 * sin_theta2 * cos_phi_diff;
        const double cos_clamped = std::clamp(cos_angular, -1.0, 1.0);
        
        const double r1 = std::max(radius_.get(), traits::origin_marker);
        const double r2 = std::max(other.radius_.get(), traits::origin_marker);
        
        const double distance_squared = r1 * r1 + r2 * r2 - 2.0 * r1 * r2 * cos_clamped;
        return safe_math::safe_sqrt(distance_squared);
    }
    
    double calculate_distance_cartesian(const UnifiedSphericalCoords& other) const {
        // Fallback Cartesian distance for numerical stability
        return (to_cartesian() - other.to_cartesian()).magnitude();
    }
};

// ============================================================================
// STATIC MEMBER DEFINITION
// ============================================================================

template<math::PrecisionLevel Level>
thread_local typename UnifiedSphericalCoords<Level>::DistanceCache 
    UnifiedSphericalCoords<Level>::distance_cache_;

// ============================================================================
// FREE FUNCTIONS AND OPERATORS
// ============================================================================

template<math::PrecisionLevel Level>
inline std::ostream& operator<<(std::ostream& os, const UnifiedSphericalCoords<Level>& coords) {
    return os << "SphericalCoords(r=" << coords.radius() 
              << ", θ=" << coords.theta() 
              << ", φ=" << coords.phi() << ")";
}

// ============================================================================
// TYPE ALIASES FOR DIFFERENT PRECISION LEVELS
// ============================================================================

using QuantumSphericalCoords = UnifiedSphericalCoords<math::PrecisionLevel::Quantum>;
using ScientificSphericalCoords = UnifiedSphericalCoords<math::PrecisionLevel::Scientific>;
using StandardSphericalCoords = UnifiedSphericalCoords<math::PrecisionLevel::Standard>;
using EngineeringSphericalCoords = UnifiedSphericalCoords<math::PrecisionLevel::Engineering>;
using GraphicsSphericalCoords = UnifiedSphericalCoords<math::PrecisionLevel::Graphics>;
using DisplaySphericalCoords = UnifiedSphericalCoords<math::PrecisionLevel::Display>;

// Default precision for backward compatibility
using SphericalCoords = StandardSphericalCoords;

} // namespace core
} // namespace hsml