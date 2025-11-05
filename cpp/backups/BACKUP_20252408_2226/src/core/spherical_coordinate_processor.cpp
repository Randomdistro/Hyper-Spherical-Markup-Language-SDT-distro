/**
 * HSML Spherical Coordinate Processor - C++20 Implementation
 * Real-time spherical mathematics engine with SIMD optimization
 * Multiple programming paradigms for maximum performance
 */

#include "hsml/core/spherical_coordinate_processor.h"
#include "hsml/core/simd_math.h"
#include <immintrin.h>
#include <algorithm>
#include <execution>
#include <numeric>
#include <chrono>
#include <thread>
#include <functional>

// [MPD Code Monkey - Security Paranoid]: C++14 compatibility helper because SOMEONE used C++17 std::clamp
namespace {
    template<typename T>
    constexpr const T& clamp_compat(const T& v, const T& lo, const T& hi) {
        return (v < lo) ? lo : (hi < v) ? hi : v;
    }
}

namespace hsml::core {

// Static member definitions
std::unique_ptr<SphericalCoordinateProcessor> SphericalCoordinateProcessor::instance_;
std::shared_mutex SphericalCoordinateProcessor::instance_mutex_;

// SIMD-accelerated distance calculations for multiple points
std::vector<double> SphericalCoordinateProcessor::spherical_distances_simd(
    const SphericalCoords& reference, 
    std::span<const SphericalCoords> points) const {
    
    std::vector<double> distances;
    distances.reserve(points.size());
    
    if (!avx2_available_ || points.size() < 4) {
        // Fallback to scalar calculation
        for (const auto& point : points) {
            distances.push_back(spherical_distance<double>(reference, point));
        }
        return distances;
    }
    
    // Load reference point into SIMD registers
    const __m256d ref_theta = _mm256_set1_pd(reference.theta());
    const __m256d ref_phi = _mm256_set1_pd(reference.phi());
    const __m256d ref_sin_theta = _mm256_set1_pd(std::sin(reference.theta()));
    const __m256d ref_cos_theta = _mm256_set1_pd(std::cos(reference.theta()));
    
    // Process 4 points at a time using AVX2
    const size_t simd_iterations = points.size() / 4;
    
    for (size_t i = 0; i < simd_iterations; ++i) {
        const size_t base_idx = i * 4;
        
        // Load 4 points' theta values
        const __m256d theta_vals = _mm256_set_pd(
            points[base_idx].theta(),
            points[base_idx + 1].theta(),
            points[base_idx + 2].theta(),
            points[base_idx + 3].theta()
        );
        
        // Load 4 points' phi values
        const __m256d phi_vals = _mm256_set_pd(
            points[base_idx].phi(),
            points[base_idx + 1].phi(),
            points[base_idx + 2].phi(),
            points[base_idx + 3].phi()
        );
        
        // Calculate sin and cos for theta values
        alignas(32) double theta_array[4];
        _mm256_store_pd(theta_array, theta_vals);
        
        const __m256d sin_theta_vals = _mm256_set_pd(
            std::sin(theta_array[0]), std::sin(theta_array[1]),
            std::sin(theta_array[2]), std::sin(theta_array[3])
        );
        
        const __m256d cos_theta_vals = _mm256_set_pd(
            std::cos(theta_array[0]), std::cos(theta_array[1]),
            std::cos(theta_array[2]), std::cos(theta_array[3])
        );
        
        // Calculate phi differences
        const __m256d phi_diff = _mm256_sub_pd(phi_vals, ref_phi);
        
        // Calculate cos(phi_diff)
        alignas(32) double phi_diff_array[4];
        _mm256_store_pd(phi_diff_array, phi_diff);
        
        const __m256d cos_phi_diff = _mm256_set_pd(
            std::cos(phi_diff_array[0]), std::cos(phi_diff_array[1]),
            std::cos(phi_diff_array[2]), std::cos(phi_diff_array[3])
        );
        
        // Spherical law of cosines: cos(d) = sin(θ1)sin(θ2)cos(φ1-φ2) + cos(θ1)cos(θ2)
        const __m256d term1 = _mm256_mul_pd(
            _mm256_mul_pd(ref_sin_theta, sin_theta_vals), 
            cos_phi_diff
        );
        const __m256d term2 = _mm256_mul_pd(ref_cos_theta, cos_theta_vals);
        const __m256d cos_distance = _mm256_add_pd(term1, term2);
        
        // Convert to distances using acos
        alignas(32) double cos_dist_array[4];
        _mm256_store_pd(cos_dist_array, cos_distance);
        
        for (int j = 0; j < 4; ++j) {
            const double clamped = clamp_compat(cos_dist_array[j], -1.0, 1.0);
            distances.push_back(std::acos(clamped));
        }
    }
    
    // Handle remaining points
    for (size_t i = simd_iterations * 4; i < points.size(); ++i) {
        distances.push_back(spherical_distance<double>(reference, points[i]));
    }
    
    // Update performance metrics
    total_calculations_.fetch_add(points.size(), std::memory_order_relaxed);
    
    return distances;
}

// Gradient calculation in spherical coordinates
SphericalCoords SphericalCoordinateProcessor::gradient_spherical(
    std::function<double(const SphericalCoords&)> scalar_field,
    const SphericalCoords& point, double epsilon) const {
    
    // Check cache first
    const std::string cache_key = std::to_string(point.r()) + "," + 
                                 std::to_string(point.theta()) + "," + 
                                 std::to_string(point.phi()) + "," + 
                                 std::to_string(epsilon);
    
    if (auto cached = gradient_cache_.get(cache_key)) {
        return *cached;
    }
    
    const auto start_time = std::chrono::high_resolution_clock::now();
    
    // Numerical gradient calculation in spherical coordinates
    const double f_center = scalar_field(point);
    
    // ∂f/∂r
    const SphericalCoords r_plus{point.r() + epsilon, point.theta(), point.phi()};
    const double df_dr = (scalar_field(r_plus) - f_center) / epsilon;
    
    // (1/r) * ∂f/∂θ
    const SphericalCoords theta_plus{point.r(), point.theta() + epsilon, point.phi()};
    const double df_dtheta = (scalar_field(theta_plus) - f_center) / (point.r() * epsilon);
    
    // (1/(r*sin(θ))) * ∂f/∂φ
    const SphericalCoords phi_plus{point.r(), point.theta(), point.phi() + epsilon};
    const double df_dphi = (scalar_field(phi_plus) - f_center) / 
                          (point.r() * std::sin(point.theta()) * epsilon);
    
    const SphericalCoords gradient{df_dr, df_dtheta, df_dphi};
    
    // Cache the result
    gradient_cache_.put(cache_key, gradient);
    
    // Update performance metrics
    const auto end_time = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration<double, std::milli>(end_time - start_time);
    
    total_calculations_.fetch_add(1, std::memory_order_relaxed);
    
    // Update average calculation time using atomic operations
    double expected = average_calculation_time_.load(std::memory_order_relaxed);
    double desired;
    do {
        const uint64_t total = total_calculations_.load(std::memory_order_relaxed);
        desired = (expected * (total - 1) + duration.count()) / total;
    } while (!average_calculation_time_.compare_exchange_weak(expected, desired, 
                                                            std::memory_order_relaxed));
    
    return gradient;
}

// Divergence calculation in spherical coordinates
double SphericalCoordinateProcessor::divergence_spherical(
    std::function<SphericalCoords(const SphericalCoords&)> vector_field,
    const SphericalCoords& point, double epsilon) const {
    
    const std::string cache_key = std::to_string(point.r()) + "," + 
                                 std::to_string(point.theta()) + "," + 
                                 std::to_string(point.phi()) + "," + 
                                 std::to_string(epsilon);
    
    if (auto cached = divergence_cache_.get(cache_key)) {
        return *cached;
    }
    
    const SphericalCoords F_center = vector_field(point);
    const double r = point.r();
    const double theta = point.theta();
    const double sin_theta = std::sin(theta);
    
    // ∂(r²Fr)/∂r term
    const SphericalCoords r_plus{r + epsilon, theta, point.phi()};
    const SphericalCoords F_r_plus = vector_field(r_plus);
    const double term1 = ((r + epsilon) * (r + epsilon) * F_r_plus.r() - r * r * F_center.r()) / 
                        (r * r * epsilon);
    
    // ∂(sin(θ)Fθ)/∂θ term
    const SphericalCoords theta_plus{r, theta + epsilon, point.phi()};
    const SphericalCoords F_theta_plus = vector_field(theta_plus);
    const double term2 = (std::sin(theta + epsilon) * F_theta_plus.theta() - 
                         sin_theta * F_center.theta()) / (r * sin_theta * epsilon);
    
    // ∂Fφ/∂φ term
    const SphericalCoords phi_plus{r, theta, point.phi() + epsilon};
    const SphericalCoords F_phi_plus = vector_field(phi_plus);
    const double term3 = (F_phi_plus.phi() - F_center.phi()) / (r * sin_theta * epsilon);
    
    const double divergence = term1 + term2 + term3;
    
    // Cache the result
    divergence_cache_.put(cache_key, divergence);
    
    total_calculations_.fetch_add(1, std::memory_order_relaxed);
    
    return divergence;
}

// Laplacian calculation in spherical coordinates
double SphericalCoordinateProcessor::laplacian_spherical(
    std::function<double(const SphericalCoords&)> scalar_field,
    const SphericalCoords& point, double epsilon) const {
    
    const std::string cache_key = std::to_string(point.r()) + "," + 
                                 std::to_string(point.theta()) + "," + 
                                 std::to_string(point.phi()) + "," + 
                                 std::to_string(epsilon);
    
    if (auto cached = laplacian_cache_.get(cache_key)) {
        return *cached;
    }
    
    const double f_center = scalar_field(point);
    const double r = point.r();
    const double theta = point.theta();
    const double sin_theta = std::sin(theta);
    const double cos_theta = std::cos(theta);
    
    // Second derivative with respect to r: (1/r²) * ∂/∂r(r² * ∂f/∂r)
    const SphericalCoords r_minus{r - epsilon, theta, point.phi()};
    const SphericalCoords r_plus{r + epsilon, theta, point.phi()};
    const double d2f_dr2 = (scalar_field(r_plus) - 2.0 * f_center + scalar_field(r_minus)) / 
                          (epsilon * epsilon);
    const double df_dr = (scalar_field(r_plus) - scalar_field(r_minus)) / (2.0 * epsilon);
    const double radial_term = d2f_dr2 + (2.0 / r) * df_dr;
    
    // Angular terms: (1/(r²sin²θ)) * [sin(θ)∂/∂θ(sin(θ)∂f/∂θ) + ∂²f/∂φ²]
    const SphericalCoords theta_minus{r, theta - epsilon, point.phi()};
    const SphericalCoords theta_plus{r, theta + epsilon, point.phi()};
    const double d2f_dtheta2 = (scalar_field(theta_plus) - 2.0 * f_center + scalar_field(theta_minus)) / 
                              (epsilon * epsilon);
    const double df_dtheta = (scalar_field(theta_plus) - scalar_field(theta_minus)) / (2.0 * epsilon);
    
    const SphericalCoords phi_minus{r, theta, point.phi() - epsilon};
    const SphericalCoords phi_plus{r, theta, point.phi() + epsilon};
    const double d2f_dphi2 = (scalar_field(phi_plus) - 2.0 * f_center + scalar_field(phi_minus)) / 
                            (epsilon * epsilon);
    
    const double angular_term = (1.0 / (r * r)) * 
                               (d2f_dtheta2 + (cos_theta / sin_theta) * df_dtheta + 
                                (1.0 / (sin_theta * sin_theta)) * d2f_dphi2);
    
    const double laplacian = radial_term + angular_term;
    
    // Cache the result
    laplacian_cache_.put(cache_key, laplacian);
    
    total_calculations_.fetch_add(1, std::memory_order_relaxed);
    
    return laplacian;
}

// SIMD-optimized solid angle calculations for multiple objects
std::vector<SolidAngle<double>> SphericalCoordinateProcessor::solid_angles_simd(
    const SphericalCoords& observer,
    std::span<const SphericalCoords> sphere_centers,
    std::span<const double> sphere_radii) const {
    
    if (sphere_centers.size() != sphere_radii.size()) {
        throw std::invalid_argument("Sphere centers and radii arrays must have the same size");
    }
    
    std::vector<SolidAngle<double>> solid_angles;
    solid_angles.reserve(sphere_centers.size());
    
    if (!avx2_available_ || sphere_centers.size() < 4) {
        // Fallback to scalar calculation
        for (size_t i = 0; i < sphere_centers.size(); ++i) {
            solid_angles.push_back(
                solid_angle_from_sphere<double>(observer, sphere_centers[i], sphere_radii[i]));
        }
        return solid_angles;
    }
    
    // SIMD processing of 4 spheres at a time
    const size_t simd_iterations = sphere_centers.size() / 4;
    
    for (size_t i = 0; i < simd_iterations; ++i) {
        const size_t base_idx = i * 4;
        
        // Calculate distances using SIMD
        std::array<SphericalCoords, 4> centers = {
            sphere_centers[base_idx], sphere_centers[base_idx + 1],
            sphere_centers[base_idx + 2], sphere_centers[base_idx + 3]
        };
        
        auto distances = spherical_distances_simd(observer, std::span{centers});
        
        // Load radii into SIMD register
        const __m256d radii = _mm256_set_pd(
            sphere_radii[base_idx], sphere_radii[base_idx + 1],
            sphere_radii[base_idx + 2], sphere_radii[base_idx + 3]
        );
        
        const __m256d distances_simd = _mm256_set_pd(
            distances[0], distances[1], distances[2], distances[3]
        );
        
        // Calculate angular radii: asin(radius / distance)
        const __m256d angular_ratio = _mm256_div_pd(radii, distances_simd);
        
        alignas(32) double ratio_array[4];
        _mm256_store_pd(ratio_array, angular_ratio);
        
        // Calculate solid angles
        for (int j = 0; j < 4; ++j) {
            if (distances[j] <= sphere_radii[base_idx + j]) {
                // Observer inside sphere
                solid_angles.emplace_back(4.0 * std::numbers::pi);
            } else {
                const double angular_radius = std::asin(clamp_compat(ratio_array[j], 0.0, 1.0));
                const double omega = 2.0 * std::numbers::pi * (1.0 - std::cos(angular_radius));
                solid_angles.emplace_back(omega);
            }
        }
    }
    
    // Handle remaining spheres
    for (size_t i = simd_iterations * 4; i < sphere_centers.size(); ++i) {
        solid_angles.push_back(
            solid_angle_from_sphere<double>(observer, sphere_centers[i], sphere_radii[i]));
    }
    
    return solid_angles;
}

// Time evolution of SDT state vector using template metaprogramming
template<SphericalMathType T>
SDTStateVector<T> SphericalCoordinateProcessor::evolve_state_vector(
    const SDTStateVector<T>& current_state, T time_step) const {
    
    SDTStateVector<T> evolved_state = current_state;
    
    // Update position using velocity (Euler integration)
    evolved_state.position = SphericalCoords{
        current_state.position.r() + current_state.velocity.v_r * time_step,
        current_state.position.theta() + current_state.velocity.v_theta * time_step,
        current_state.position.phi() + current_state.velocity.v_phi * time_step
    };
    
    // Update velocity using acceleration
    evolved_state.velocity = SphericalVelocity<T>{
        current_state.velocity.v_r + current_state.acceleration.a_r * time_step,
        current_state.velocity.v_theta + current_state.acceleration.a_theta * time_step,
        current_state.velocity.v_phi + current_state.acceleration.a_phi * time_step
    };
    
    // Update electromagnetic fields (simplified evolution)
    evolved_state.electromagnetic_field.E_r *= std::exp(-time_step * 0.1); // Damping
    evolved_state.electromagnetic_field.E_theta *= std::exp(-time_step * 0.1);
    evolved_state.electromagnetic_field.E_phi *= std::exp(-time_step * 0.1);
    
    return evolved_state;
}

// Runge-Kutta integration for spherical dynamics
template<SphericalMathType T>
SDTStateVector<T> SphericalCoordinateProcessor::runge_kutta_step(
    const SDTStateVector<T>& state, T dt,
    std::function<SDTStateVector<T>(const SDTStateVector<T>&, T)> derivative_func,
    T current_time) const {
    
    // Fourth-order Runge-Kutta method
    const auto k1 = derivative_func(state, current_time);
    const auto k2 = derivative_func(state + k1 * (dt / T{2}), current_time + dt / T{2});
    const auto k3 = derivative_func(state + k2 * (dt / T{2}), current_time + dt / T{2});
    const auto k4 = derivative_func(state + k3 * dt, current_time + dt);
    
    return state + (k1 + k2 * T{2} + k3 * T{2} + k4) * (dt / T{6});
}

// Performance monitoring
SphericalCoordinateProcessor::PerformanceMetrics 
SphericalCoordinateProcessor::get_performance_metrics() const {
    return PerformanceMetrics{
        .total_calculations = total_calculations_.load(std::memory_order_relaxed),
        .average_calculation_time_ms = average_calculation_time_.load(std::memory_order_relaxed),
        .gradient_cache_hit_ratio = gradient_cache_.hit_ratio(),
        .divergence_cache_hit_ratio = divergence_cache_.hit_ratio(),
        .laplacian_cache_hit_ratio = laplacian_cache_.hit_ratio(),
        .simd_enabled = simd_enabled_,
        .avx2_available = avx2_available_
    };
}

// Cache management
void SphericalCoordinateProcessor::clear_caches() {
    // Since we're using lock-free caches, we can't easily clear them
    // In a production implementation, we'd implement a proper clear method
    total_calculations_.store(0, std::memory_order_relaxed);
    average_calculation_time_.store(0.0, std::memory_order_relaxed);
}

void SphericalCoordinateProcessor::optimize_cache_size(size_t target_hit_ratio_percent) {
    // Cache optimization logic would go here
    // This is a placeholder for a more sophisticated implementation
    const auto metrics = get_performance_metrics();
    
    if (metrics.gradient_cache_hit_ratio * 100 < target_hit_ratio_percent) {
        // Increase cache size or adjust eviction policy
    }
}

} // namespace hsml::core

// Explicit template instantiations for common floating-point types
template hsml::core::SDTStateVector<float> 
hsml::core::SphericalCoordinateProcessor::evolve_state_vector<float>(
    const hsml::core::SDTStateVector<float>&, float) const;

template hsml::core::SDTStateVector<double> 
hsml::core::SphericalCoordinateProcessor::evolve_state_vector<double>(
    const hsml::core::SDTStateVector<double>&, double) const;

template hsml::core::SDTStateVector<float> 
hsml::core::SphericalCoordinateProcessor::runge_kutta_step<float>(
    const hsml::core::SDTStateVector<float>&, float,
    std::function<hsml::core::SDTStateVector<float>(const hsml::core::SDTStateVector<float>&, float)>,
    float) const;

template hsml::core::SDTStateVector<double> 
hsml::core::SphericalCoordinateProcessor::runge_kutta_step<double>(
    const hsml::core::SDTStateVector<double>&, double,
    std::function<hsml::core::SDTStateVector<double>(const hsml::core::SDTStateVector<double>&, double)>,
    double) const;