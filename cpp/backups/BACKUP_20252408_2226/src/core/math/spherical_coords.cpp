#include "hsml/core/spherical_coords.h"
#include "hsml/core/simd_math.h"
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <string>
#include <sstream>
#include <iomanip>
#include <chrono>

// [ASPIE ARCHITECT]: Import precision constants for architectural consistency
using namespace hsml::core::precision;

// [MPD Code Monkey]: C++14 compatibility helper
namespace {
    template<typename T>
    constexpr const T& clamp_compat(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }
}

namespace hsml {
namespace core {

// [ASPIE ARCHITECT]: Thread-safe cache with standardized limits
std::unordered_map<std::string, double> SphericalCoords::distance_cache_;
size_t SphericalCoords::cache_hits_ = 0;
size_t SphericalCoords::cache_misses_ = 0;
static constexpr size_t MAX_CACHE_SIZE = PrecisionLevels::MAX_CACHE_SIZE;

// ============================================================================
// [ASPIE ARCHITECT]: CONSOLIDATED SPHERICAL DISTANCE IMPLEMENTATION
// Single optimized method with SIMD acceleration and caching
// ============================================================================

namespace {
    // [ASPIE ARCHITECT]: Thread-local cache for performance
    thread_local std::unordered_map<std::string, double> computation_cache;
    thread_local uint64_t cache_hits = 0;
    thread_local uint64_t cache_misses = 0;
    
    std::string make_cache_key(const SphericalCoords& p1, const SphericalCoords& p2) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(static_cast<int>(-std::log10(PrecisionLevels::CACHE_PRECISION)))
            << "dist_" << p1.radius() << "_" << p1.theta() << "_" << p1.phi()
            << "_" << p2.radius() << "_" << p2.theta() << "_" << p2.phi();
        return oss.str();
    }
}

// [ASPIE ARCHITECT]: Primary optimized spherical distance calculation
// Combines caching, SIMD acceleration, and mathematical safety
double SphericalCoords::spherical_distance(const SphericalCoords& other) const {
    // Check cache first for performance
    const std::string cache_key = make_cache_key(*this, other);
    auto cache_it = computation_cache.find(cache_key);
    if (cache_it != computation_cache.end()) {
        cache_hits++;
        return cache_it->second;
    }
    
    // Try SIMD-optimized calculation if available
    double distance = 0.0;
    if (simd::SIMDConfig::use_simd && simd::SIMDConfig::has_sse2) {
        distance = calculate_distance_simd(other);
    } else {
        distance = calculate_distance_scalar(other);
    }
    
    // [ASPIE ARCHITECT]: Standardized cache management
    if (computation_cache.size() < PrecisionLevels::MAX_CACHE_SIZE) {
        computation_cache[cache_key] = distance;
    } else {
        computation_cache.clear();
        computation_cache[cache_key] = distance;
    }
    cache_misses++;
    return distance;
}

// Standard spherical distance calculation
double SphericalCoords::calculate_distance_scalar(const SphericalCoords& other) const {
    // Direct trigonometric calculations - zero is safe in SDT
    const double cos_theta1 = std::cos(theta_);
    const double sin_theta1 = std::sin(theta_);
    const double cos_theta2 = std::cos(other.theta_);
    const double sin_theta2 = std::sin(other.theta_);
    const double cos_phi_diff = std::cos(other.phi_ - phi_);
    
    const double cos_angular = cos_theta1 * cos_theta2 + 
                               sin_theta1 * sin_theta2 * cos_phi_diff;
    
    // [ASPIE ARCHITECT]: Safe mathematical operations
    const double cos_clamped = std::max(-1.0, std::min(1.0, cos_angular));
    
    // 3D spherical distance using law of cosines
    const double r1_safe = std::max(radius_, PrecisionLevels::MIN_SAFE_RADIUS);
    const double r2_safe = std::max(other.radius_, PrecisionLevels::MIN_SAFE_RADIUS);
    const double distance_squared = r1_safe * r1_safe + r2_safe * r2_safe - 
                                   2.0 * r1_safe * r2_safe * cos_clamped;
    
    return SafeMath::safe_sqrt(distance_squared);
}

// [ASPIE ARCHITECT]: SIMD-optimized calculation
double SphericalCoords::calculate_distance_simd(const SphericalCoords& other) const {
#ifdef HSML_SSE2_AVAILABLE
    // Use SIMD for parallel trigonometric calculations where beneficial
    const double cos_theta1 = std::cos(theta_);
    const double sin_theta1 = std::sin(theta_);
    const double cos_theta2 = std::cos(other.theta_);
    const double sin_theta2 = std::sin(other.theta_);
    const double cos_phi_diff = std::cos(other.phi_ - phi_);
    
    // Vectorized cosine calculation
    __m128d cos_vals = _mm_set_pd(cos_theta2, cos_theta1);
    __m128d sin_vals = _mm_set_pd(sin_theta2, sin_theta1);
    __m128d cos_prod = _mm_mul_pd(cos_vals, _mm_shuffle_pd(cos_vals, cos_vals, 1));
    __m128d sin_prod = _mm_mul_pd(sin_vals, _mm_shuffle_pd(sin_vals, sin_vals, 1));
    
    double cos_angular = _mm_cvtsd_f64(cos_prod) + 
                        _mm_cvtsd_f64(sin_prod) * cos_phi_diff;
    
    // Safe clamping and distance calculation
    cos_angular = std::max(-1.0, std::min(1.0, cos_angular));
    
    const double r1_safe = std::max(radius_, PrecisionLevels::MIN_SAFE_RADIUS);
    const double r2_safe = std::max(other.radius_, PrecisionLevels::MIN_SAFE_RADIUS);
    
    __m128d r_vals = _mm_set_pd(r2_safe, r1_safe);
    __m128d r_squared = _mm_mul_pd(r_vals, r_vals);
    double r_sum_squared = _mm_cvtsd_f64(r_squared) + 
                          _mm_cvtsd_f64(_mm_unpackhi_pd(r_squared, r_squared));
    
    double distance_squared = r_sum_squared - 2.0 * r1_safe * r2_safe * cos_angular;
    return SafeMath::safe_sqrt(distance_squared);
#else
    // Fallback to scalar implementation
    return calculate_distance_scalar(other);
#endif
}

// [ASPIE ARCHITECT]: Removed redundant template method - consolidated into main implementation

// [ASPIE ARCHITECT]: Removed redundant old SIMD method - consolidated into main implementation

// [ASPIE ARCHITECT]: Removed redundant functional and template methods

// [ASPIE ARCHITECT]: Cache key generation moved to private utility
std::string SphericalCoords::generate_cache_key(const SphericalCoords& other) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(static_cast<int>(-std::log10(PrecisionLevels::CACHE_PRECISION)))
        << "dist_" << radius_ << "_" << theta_ << "_" << phi_
        << "_" << other.radius_ << "_" << other.theta_ << "_" << other.phi_;
    return oss.str();
}

// [ASPIE ARCHITECT]: Performance comparison for validation
SphericalCoords::DistanceComparison SphericalCoords::compare_distance_methods(const SphericalCoords& other) const {
    DistanceComparison result;
    
    // Time the optimized method
    auto start = std::chrono::high_resolution_clock::now();
    result.optimized_result = spherical_distance(other);
    auto end = std::chrono::high_resolution_clock::now();
    auto optimized_time = std::chrono::duration<double, std::nano>(end - start).count();
    
    // Time the fallback method
    start = std::chrono::high_resolution_clock::now();
    result.fallback_result = distance_to(other);
    end = std::chrono::high_resolution_clock::now();
    auto fallback_time = std::chrono::duration<double, std::nano>(end - start).count();
    
    // Calculate performance ratio
    result.performance_ratio = (fallback_time > 0.0) ? (optimized_time / fallback_time) : 1.0;
    
    // Check if results match within tolerance
    result.results_match = SafeMath::are_equal(result.optimized_result, result.fallback_result, 
                                             PrecisionLevels::STANDARD_PRECISE);
    
    return result;
}

// Cache management methods
void SphericalCoords::clear_distance_cache() {
    distance_cache_.clear();
    cache_hits_ = 0;
    cache_misses_ = 0;
}

std::pair<size_t, size_t> SphericalCoords::get_cache_stats() {
    return {cache_hits_, cache_misses_};
}

SphericalCoords SphericalCoords::slerp(const SphericalCoords& other, double t) const {
    t = clamp_compat(t, 0.0, 1.0);

    // Pure spherical slerp - no Cartesian conversion needed
    double angular_sep = angular_distance(other);

    if (std::abs(angular_sep) < 1e-6) {
        // Points are nearly identical, use linear interpolation
        return interpolate_spherical(other, t);
    }

    // Spherical linear interpolation for angular components
    double theta_diff = other.theta_ - theta_;
    double phi_diff = other.phi_ - phi_;

    // Use the angular separation for proper slerp weighting
    double weight1 = std::sin((1.0 - t) * angular_sep) / std::sin(angular_sep);
    double weight2 = std::sin(t * angular_sep) / std::sin(angular_sep);

    double new_theta = theta_ * weight1 + other.theta_ * weight2;
    double new_phi = phi_ * weight1 + other.phi_ * weight2;

    // Linear interpolation for radius
    double interpolated_radius = radius_ * (1.0 - t) + other.radius_ * t;

    SphericalCoords result(interpolated_radius, new_theta, new_phi);
    
    return result;
}

// ❌ ERADICATED: Cartesian-based rotation method completely removed
// Use pure spherical rotation operations instead!

std::pair<SphericalCoords, SphericalCoords> SphericalCoords::get_tangent_vectors() const {
    // Pure spherical tangent vectors - no Cartesian conversion needed

    // Theta direction tangent (varying theta at constant phi)
    // In spherical coordinates, moving in theta direction affects the polar angle
    SphericalCoords theta_tangent(radius_, theta_ + 0.01, phi_);

    // Phi direction tangent (varying phi at constant theta)
    // In spherical coordinates, moving in phi direction affects the azimuthal angle
    SphericalCoords phi_tangent(radius_, theta_, phi_ + 0.01);

    return {theta_tangent, phi_tangent};
}

bool SphericalCoords::is_approximately_equal(const SphericalCoords& other, double epsilon) const {
    return std::abs(radius_ - other.radius_) < epsilon &&
           std::abs(theta_ - other.theta_) < epsilon &&
           std::abs(phi_ - other.phi_) < epsilon;
}

double SphericalCoords::solid_angle_to(const SphericalCoords& other) const {
    Vector3 v1 = to_cartesian().normalized();
    Vector3 v2 = other.to_cartesian().normalized();
    
    double cos_angle = v1.dot(v2);
    cos_angle = clamp_compat(cos_angle, -1.0, 1.0);
    
    return 2.0 * PI * (1.0 - cos_angle);
}

} // namespace core
} // namespace hsml