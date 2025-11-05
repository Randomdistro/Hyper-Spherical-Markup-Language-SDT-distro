#pragma once

#include "../spherical_purity_guard.h"  // 🛡️ ANTI-CARTESIAN PROTECTION
#include "vector3.h"
#include "precision_constants.h"
#include <cmath>
#include <ostream>
#include <algorithm>
#include <utility>
#include <unordered_map>
#include <string>
#include <sstream>
#include <iomanip>

// [MPD Code Monkey]: C++14 compatibility helper
namespace hsml_clamp_ns {
    template<typename T>
    constexpr const T& clamp_compat(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }
}

namespace hsml {
namespace core {

#if defined(HSML_STRICT_SPHERICAL)
#  define HSML_SPHERICAL_DEPRECATED [[deprecated("Cartesian bridge disabled under HSML_STRICT_SPHERICAL")]]
#else
#  define HSML_SPHERICAL_DEPRECATED
#endif

class SphericalCoords {
public:
    // [ASPIE ARCHITECT]: Use centralized mathematical constants for consistency
    using MathConstants = precision::MathematicalConstants;
    using Precision = precision::PrecisionLevels;
    using SafeMath = precision::SafeMath;
    
    static constexpr double PI = MathConstants::PI;
    static constexpr double TWO_PI = MathConstants::TWO_PI;
    static constexpr double HALF_PI = MathConstants::HALF_PI;
    
    // O = Origin (0,0,0) - SDT has no Cartesian so zero is safe
    SphericalCoords() : radius_(0.0), theta_(0.0), phi_(0.0) {}
    SphericalCoords(double radius, double theta, double phi) 
        : radius_(radius), theta_(theta), phi_(phi) {
        normalize_angles();
    }
    
    double radius() const { return radius_; }
    double theta() const { return theta_; }
    double phi() const { return phi_; }
    
    void set_radius(double radius) { radius_ = std::max(0.0, radius); }
    void set_theta(double theta) { theta_ = theta; normalize_angles(); }
    void set_phi(double phi) { phi_ = phi; normalize_angles(); }
    
    // ❌ ERADICATED: Cartesian conversion method completely removed
    // Use pure spherical operations for all calculations!
    
    // ❌ ERADICATED: Cartesian conversion method completely removed
    // Start with spherical coordinates for all operations!
    
    // [ASPIE ARCHITECT]: Consolidated distance calculations - ONLY 2 methods for architectural purity
    
    // Primary optimized distance calculation with caching and SIMD acceleration
    double spherical_distance(const SphericalCoords& other) const;
    
    
    
    // [ASPIE ARCHITECT]: Performance comparison method for validation only
    struct DistanceComparison {
        double optimized_result;
        double fallback_result;
        double performance_ratio;
        bool results_match;
    };
    DistanceComparison compare_distance_methods(const SphericalCoords& other) const;
    
private:
    // [ASPIE ARCHITECT]: Internal calculation methods
    double calculate_distance_scalar(const SphericalCoords& other) const;
    double calculate_distance_simd(const SphericalCoords& other) const;
    
public:
    
    double angular_distance(const SphericalCoords& other) const {
        double cos_angle = std::sin(theta_) * std::sin(other.theta_) * std::cos(phi_ - other.phi_) +
                          std::cos(theta_) * std::cos(other.theta_);
        using hsml_clamp_ns::clamp_compat;
        return std::acos(clamp_compat(cos_angle, -1.0, 1.0));
    }
    
    // Pure spherical interpolation - no Cartesian conversion needed
    SphericalCoords interpolate_spherical(const SphericalCoords& other, double t) const {
        if (t <= 0.0) return *this;
        if (t >= 1.0) return other;

        // Spherical linear interpolation (slerp) for orientation
        double theta_diff = other.theta_ - theta_;
        double phi_diff = other.phi_ - phi_;

        double new_theta = theta_ + theta_diff * t;
        double new_phi = phi_ + phi_diff * t;

        // Linear interpolation for radius
        double new_radius = radius_ + (other.radius_ - radius_) * t;

        return SphericalCoords(new_radius, new_theta, new_phi);
    }
    
    bool is_valid() const {
        return std::isfinite(radius_) && std::isfinite(theta_) && std::isfinite(phi_) && 
               radius_ >= 0.0 && theta_ >= 0.0 && theta_ <= PI;
    }
    
    // Additional utility methods
    SphericalCoords slerp(const SphericalCoords& other, double t) const;
    SphericalCoords rotate_about_axis(const Vector3& axis, double angle) const;
    std::pair<SphericalCoords, SphericalCoords> get_tangent_vectors() const;
    bool is_approximately_equal(const SphericalCoords& other, double epsilon = Precision::STANDARD_PRECISE) const;
    double solid_angle_to(const SphericalCoords& other) const;
    
    // Cache management methods (transposed from TypeScript)
    static void clear_distance_cache();
    static std::pair<size_t, size_t> get_cache_stats(); // returns {hits, misses}
    
 
    
private:
    // O = Origin (0,0,0) - true mathematical origin
    double radius_{0.0};   // True zero radius
    double theta_{0.0};    // True zero polar angle  
    double phi_{0.0};      // True zero azimuthal angle
    
    // Static cache for distance calculations (transposed from TypeScript)
    static std::unordered_map<std::string, double> distance_cache_;
    static size_t cache_hits_;
    static size_t cache_misses_;
    
    void normalize_angles() {
        // Standard spherical coordinate normalization
        // φ: [0, 2π], θ: [0, π]
        while (phi_ > TWO_PI) phi_ -= TWO_PI;
        while (phi_ < 0.0) phi_ += TWO_PI;
        
        using hsml_clamp_ns::clamp_compat;
        theta_ = clamp_compat(theta_, 0.0, PI);
    }
    
    // Helper method to generate cache key
    std::string generate_cache_key(const SphericalCoords& other) const;
};

inline std::ostream& operator<<(std::ostream& os, const SphericalCoords& coords) {
    return os << "SphericalCoords(r=" << coords.radius() 
              << ", θ=" << coords.theta() 
              << ", φ=" << coords.phi() << ")";
}

} // namespace core
} // namespace hsml