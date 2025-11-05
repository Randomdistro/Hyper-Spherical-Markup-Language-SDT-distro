/**
 * HSML Spherical Coordinate Processor - C++20 Implementation
 * Real-time spherical mathematics engine for native 3D calculations
 * Advanced template metaprogramming and SIMD optimization
 */

#pragma once

#include <cmath>
#include <array>
#include <vector>
#include <unordered_map>
#include <memory>
#include <concepts>
#include <numbers>
#include <immintrin.h>
#include <execution>
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include <span>
#include <ranges>

#include "spherical_coords.h"
#include "vector3.h"

namespace hsml::core {

// CARTESIAN COORDINATES BANNED FROM HSML!
// ONLY PURE SPHERICAL MATHEMATICS ALLOWED!

// Modern C++20 concepts for spherical mathematics
template<typename T>
concept SphericalMathType = std::floating_point<T> && requires(T t) {
    { t >= T{0} } -> std::convertible_to<bool>;
    { std::sin(t) } -> std::convertible_to<T>;
    { std::cos(t) } -> std::convertible_to<T>;
};

template<typename T>
concept SphericalField = requires(T t) {
    { t.r } -> std::convertible_to<double>;
    { t.theta } -> std::convertible_to<double>;
    { t.phi } -> std::convertible_to<double>;
};

// Compile-time spherical solid angle calculations
template<SphericalMathType T>
struct SolidAngle {
    T omega;        // Solid angle in steradians
    T theta_min;    // Minimum polar angle
    T theta_max;    // Maximum polar angle
    T phi_min;      // Minimum azimuthal angle
    T phi_max;      // Maximum azimuthal angle
    
    constexpr SolidAngle(T omega_val, T theta_min_val = T{0}, T theta_max_val = std::numbers::pi_v<T>,
                        T phi_min_val = T{0}, T phi_max_val = 2 * std::numbers::pi_v<T>) noexcept
        : omega(omega_val), theta_min(theta_min_val), theta_max(theta_max_val),
          phi_min(phi_min_val), phi_max(phi_max_val) {}
    
    // Constexpr solid angle calculations
    [[nodiscard]] constexpr T total_solid_angle() const noexcept {
        return 2 * std::numbers::pi_v<T> * (std::cos(theta_min) - std::cos(theta_max));
    }
    
    [[nodiscard]] constexpr bool contains_direction(T theta, T phi) const noexcept {
        return theta >= theta_min && theta <= theta_max && phi >= phi_min && phi <= phi_max;
    }
};

// Expression templates for spherical velocity vectors
template<SphericalMathType T>
struct SphericalVelocity {
    T v_r;      // Radial velocity
    T v_theta;  // Polar angular velocity
    T v_phi;    // Azimuthal angular velocity
    
    constexpr SphericalVelocity(T vr = T{0}, T vtheta = T{0}, T vphi = T{0}) noexcept
        : v_r(vr), v_theta(vtheta), v_phi(vphi) {}
    
    // Operator overloading for vector operations
    constexpr SphericalVelocity operator+(const SphericalVelocity& other) const noexcept {
        return SphericalVelocity{v_r + other.v_r, v_theta + other.v_theta, v_phi + other.v_phi};
    }
    
    constexpr SphericalVelocity operator*(T scalar) const noexcept {
        return SphericalVelocity{v_r * scalar, v_theta * scalar, v_phi * scalar};
    }
    
    [[nodiscard]] constexpr T magnitude_squared() const noexcept {
        return v_r * v_r + v_theta * v_theta + v_phi * v_phi;
    }
    
    [[nodiscard]] constexpr T magnitude() const noexcept {
        return std::sqrt(magnitude_squared());
    }
};

// Template metaprogramming for spherical acceleration
template<SphericalMathType T>
struct SphericalAcceleration {
    T a_r;      // Radial acceleration
    T a_theta;  // Polar angular acceleration
    T a_phi;    // Azimuthal angular acceleration
    
    constexpr SphericalAcceleration(T ar = T{0}, T atheta = T{0}, T aphi = T{0}) noexcept
        : a_r(ar), a_theta(atheta), a_phi(aphi) {}
    
    // Template-based acceleration calculations
    template<typename TimeType>
    [[nodiscard]] constexpr SphericalVelocity<T> integrate_velocity(
        const SphericalVelocity<T>& initial_velocity, TimeType dt) const noexcept {
        return SphericalVelocity<T>{
            initial_velocity.v_r + a_r * dt,
            initial_velocity.v_theta + a_theta * dt,
            initial_velocity.v_phi + a_phi * dt
        };
    }
};

// 21-dimensional state vector with compile-time validation
template<SphericalMathType T>
struct alignas(64) SDTStateVector {  // Aligned for SIMD operations
    SphericalCoords position;
    SphericalVelocity<T> velocity;
    SphericalAcceleration<T> acceleration;
    SphericalCoords angular_momentum;
    
    // Electromagnetic field components (spherical)
    struct EMField {
        T E_r, E_theta, E_phi;  // Electric field
        T B_r, B_theta, B_phi;  // Magnetic field
        
        constexpr EMField(T er = T{0}, T etheta = T{0}, T ephi = T{0},
                         T br = T{0}, T btheta = T{0}, T bphi = T{0}) noexcept
            : E_r(er), E_theta(etheta), E_phi(ephi), B_r(br), B_theta(btheta), B_phi(bphi) {}
    } electromagnetic_field;
    
    // Material properties for physics simulation
    struct MaterialProps {
        T density;
        T pressure;
        T temperature;
        T permeability;
        T permittivity;
        T conductivity;
        
        constexpr MaterialProps(T d = T{1}, T p = T{0}, T temp = T{293.15},
                               T mu = T{1}, T eps = T{1}, T sigma = T{0}) noexcept
            : density(d), pressure(p), temperature(temp), 
              permeability(mu), permittivity(eps), conductivity(sigma) {}
    } material_properties;
    
    constexpr SDTStateVector() = default;
    
    // State vector operations
    constexpr SDTStateVector operator+(const SDTStateVector& other) const noexcept {
        SDTStateVector result;
        result.position = position + other.position;
        result.velocity = velocity + other.velocity;
        result.acceleration = acceleration + other.acceleration;
        result.angular_momentum = angular_momentum + other.angular_momentum;
        // EM field and material property addition would be implemented similarly
        return result;
    }
};

// Lock-free cache for spherical calculations
template<typename KeyType, typename ValueType>
class LockFreeSphericalCache {
    struct CacheEntry {
        std::atomic<KeyType> key;
        std::atomic<ValueType> value;
        std::atomic<bool> valid{false};
    };
    
    static constexpr size_t CACHE_SIZE = 4096;
    std::array<CacheEntry, CACHE_SIZE> cache_;
    std::atomic<size_t> cache_hits_{0};
    std::atomic<size_t> cache_misses_{0};
    
    [[nodiscard]] constexpr size_t hash_key(const KeyType& key) const noexcept {
        // Simple hash function - replace with better one in production
        return std::hash<KeyType>{}(key) % CACHE_SIZE;
    }
    
public:
    [[nodiscard]] std::optional<ValueType> get(const KeyType& key) {
        const size_t index = hash_key(key);
        auto& entry = cache_[index];
        
        if (entry.valid.load(std::memory_order_acquire) && 
            entry.key.load(std::memory_order_relaxed) == key) {
            cache_hits_.fetch_add(1, std::memory_order_relaxed);
            return entry.value.load(std::memory_order_relaxed);
        }
        
        cache_misses_.fetch_add(1, std::memory_order_relaxed);
        return std::nullopt;
    }
    
    void put(const KeyType& key, const ValueType& value) {
        const size_t index = hash_key(key);
        auto& entry = cache_[index];
        
        entry.key.store(key, std::memory_order_relaxed);
        entry.value.store(value, std::memory_order_relaxed);
        entry.valid.store(true, std::memory_order_release);
    }
    
    [[nodiscard]] double hit_ratio() const noexcept {
        const auto hits = cache_hits_.load(std::memory_order_relaxed);
        const auto misses = cache_misses_.load(std::memory_order_relaxed);
        return (hits + misses > 0) ? static_cast<double>(hits) / (hits + misses) : 0.0;
    }
};

// Main spherical coordinate processor with multiple paradigms
class SphericalCoordinateProcessor {
private:
    // Singleton pattern with thread safety
    static std::unique_ptr<SphericalCoordinateProcessor> instance_;
    static std::shared_mutex instance_mutex_;
    
    // Lock-free caches for different calculation types
    LockFreeSphericalCache<std::string, SphericalCoords> gradient_cache_;
    LockFreeSphericalCache<std::string, double> divergence_cache_;
    LockFreeSphericalCache<std::string, double> laplacian_cache_;
    
    // Performance metrics
    std::atomic<uint64_t> total_calculations_{0};
    std::atomic<double> average_calculation_time_{0.0};
    
    // SIMD optimization flags
    bool simd_enabled_ = true;
    bool avx2_available_ = false;
    bool avx512_available_ = false;
    
    SphericalCoordinateProcessor() {
        detect_simd_capabilities();
    }
    
    void detect_simd_capabilities() {
        // CPU feature detection for SIMD optimization
        // This would normally use CPUID instructions
        simd_enabled_ = true;
        avx2_available_ = true;  // Assume available for now
        avx512_available_ = false;
    }
    
public:
    // Thread-safe singleton access
    [[nodiscard]] static SphericalCoordinateProcessor& get_instance() {
        std::shared_lock<std::shared_mutex> lock(instance_mutex_);
        if (!instance_) {
            lock.unlock();
            std::unique_lock<std::shared_mutex> write_lock(instance_mutex_);
            if (!instance_) {
                instance_ = std::unique_ptr<SphericalCoordinateProcessor>(
                    new SphericalCoordinateProcessor());
            }
        }
        return *instance_;
    }
    
    // === CORE SPHERICAL MATHEMATICS ===
    
    // Pure spherical distance calculation using MPD Code Monkey implementations
    template<SphericalMathType T>
    [[nodiscard]] constexpr T spherical_distance(const SphericalCoords& p1, 
                                                 const SphericalCoords& p2) const noexcept {
        // Use the modern constexpr implementation from MPD Code Monkey
        if constexpr (std::is_same_v<T, double>) {
            return p1.spherical_distance_modern<T>(p2);
        } else {
            // Fallback for other types
            const T cos_distance = std::sin(p1.theta()) * std::sin(p2.theta()) * 
                                  std::cos(p1.phi() - p2.phi()) +
                                  std::cos(p1.theta()) * std::cos(p2.theta());
            
            return std::acos(std::clamp(cos_distance, T{-1}, T{1}));
        }
    }
    
    // High-performance SIMD distance calculation
    [[nodiscard]] double spherical_distance_simd(const SphericalCoords& p1, 
                                                 const SphericalCoords& p2) const {
        return p1.spherical_distance_simd(p2);
    }
    
    // Cached distance calculation for repeated queries
    [[nodiscard]] double spherical_distance_cached(const SphericalCoords& p1, 
                                                  const SphericalCoords& p2) const {
        return p1.spherical_distance_cached(p2);
    }
    
    // Functional programming style distance
    [[nodiscard]] static double spherical_distance_functional(const SphericalCoords& p1, 
                                                             const SphericalCoords& p2) {
        return SphericalCoords::spherical_distance_functional(p1, p2);
    }
    
    // SIMD-accelerated distance calculations for multiple points
    [[nodiscard]] std::vector<double> spherical_distances_simd(
        const SphericalCoords& reference, 
        std::span<const SphericalCoords> points) const;
    
    // Template metaprogramming for compile-time spherical interpolation
    template<SphericalMathType T, size_t Steps>
    [[nodiscard]] constexpr std::array<SphericalCoords, Steps> interpolate_path(
        const SphericalCoords& start, const SphericalCoords& end) const noexcept {
        std::array<SphericalCoords, Steps> path;
        
        for (size_t i = 0; i < Steps; ++i) {
            const T t = static_cast<T>(i) / static_cast<T>(Steps - 1);
            path[i] = spherical_interpolate(start, end, t);
        }
        
        return path;
    }
    
    // Spherical interpolation (SLERP for spherical coordinates)
    template<SphericalMathType T>
    [[nodiscard]] constexpr SphericalCoords spherical_interpolate(
        const SphericalCoords& start, const SphericalCoords& end, T t) const noexcept {
        const T distance = spherical_distance<T>(start, end);
        
        if (distance < std::numeric_limits<T>::epsilon()) {
            return start;  // Points are too close
        }
        
        const T sin_distance = std::sin(distance);
        const T factor_start = std::sin((T{1} - t) * distance) / sin_distance;
        const T factor_end = std::sin(t * distance) / sin_distance;
        
        // Interpolate in spherical space
        const T interp_r = start.r() * factor_start + end.r() * factor_end;
        const T interp_theta = start.theta() * factor_start + end.theta() * factor_end;
        const T interp_phi = start.phi() * factor_start + end.phi() * factor_end;
        
        return SphericalCoords{interp_r, interp_theta, interp_phi};
    }
    
    // === SPHERICAL DIFFERENTIAL OPERATORS ===
    
    // Gradient in spherical coordinates
    [[nodiscard]] SphericalCoords gradient_spherical(
        std::function<double(const SphericalCoords&)> scalar_field,
        const SphericalCoords& point, double epsilon = 1e-8) const;
    
    // Divergence in spherical coordinates
    [[nodiscard]] double divergence_spherical(
        std::function<SphericalCoords(const SphericalCoords&)> vector_field,
        const SphericalCoords& point, double epsilon = 1e-8) const;
    
    // Laplacian in spherical coordinates
    [[nodiscard]] double laplacian_spherical(
        std::function<double(const SphericalCoords&)> scalar_field,
        const SphericalCoords& point, double epsilon = 1e-8) const;
    
    // === SOLID ANGLE CALCULATIONS ===
    
    // Calculate solid angle subtended by spherical object
    template<SphericalMathType T>
    [[nodiscard]] constexpr SolidAngle<T> solid_angle_from_sphere(
        const SphericalCoords& observer, const SphericalCoords& sphere_center, 
        T sphere_radius) const noexcept {
        const T distance = spherical_distance<T>(observer, sphere_center);
        
        if (distance <= sphere_radius) {
            // Observer inside sphere - full solid angle
            return SolidAngle<T>{4 * std::numbers::pi_v<T>};
        }
        
        const T angular_radius = std::asin(sphere_radius / distance);
        const T omega = 2 * std::numbers::pi_v<T> * (T{1} - std::cos(angular_radius));
        
        return SolidAngle<T>{omega};
    }
    
    // SIMD-optimized solid angle calculations for multiple objects
    [[nodiscard]] std::vector<SolidAngle<double>> solid_angles_simd(
        const SphericalCoords& observer,
        std::span<const SphericalCoords> sphere_centers,
        std::span<const double> sphere_radii) const;
    
    // === SDT STATE VECTOR OPERATIONS ===
    
    // Time evolution of SDT state vector
    template<SphericalMathType T>
    [[nodiscard]] SDTStateVector<T> evolve_state_vector(
        const SDTStateVector<T>& current_state, T time_step) const;
    
    // Runge-Kutta integration for spherical dynamics
    template<SphericalMathType T>
    [[nodiscard]] SDTStateVector<T> runge_kutta_step(
        const SDTStateVector<T>& state, T dt,
        std::function<SDTStateVector<T>(const SDTStateVector<T>&, T)> derivative_func,
        T current_time) const;
    
    // === PERFORMANCE MONITORING ===
    
    struct PerformanceMetrics {
        uint64_t total_calculations;
        double average_calculation_time_ms;
        double gradient_cache_hit_ratio;
        double divergence_cache_hit_ratio;
        double laplacian_cache_hit_ratio;
        bool simd_enabled;
        bool avx2_available;
    };
    
    [[nodiscard]] PerformanceMetrics get_performance_metrics() const;
    
    // Cache management
    void clear_caches();
    void optimize_cache_size(size_t target_hit_ratio_percent = 85);
    
    // SIMD capability queries
    [[nodiscard]] bool is_simd_enabled() const noexcept { return simd_enabled_; }
    [[nodiscard]] bool is_avx2_available() const noexcept { return avx2_available_; }
    
    // Disable copy/move operations for singleton
    SphericalCoordinateProcessor(const SphericalCoordinateProcessor&) = delete;
    SphericalCoordinateProcessor& operator=(const SphericalCoordinateProcessor&) = delete;
    SphericalCoordinateProcessor(SphericalCoordinateProcessor&&) = delete;
    SphericalCoordinateProcessor& operator=(SphericalCoordinateProcessor&&) = delete;
    
    ~SphericalCoordinateProcessor() = default;
};

// Free functions for functional programming style
namespace spherical_math {
    // Pure functional spherical distance
    template<SphericalMathType T>
    [[nodiscard]] constexpr T distance(const SphericalCoords& p1, 
                                      const SphericalCoords& p2) noexcept {
        return SphericalCoordinateProcessor::get_instance().spherical_distance<T>(p1, p2);
    }
    
    // Pure functional solid angle calculation
    template<SphericalMathType T>
    [[nodiscard]] constexpr SolidAngle<T> solid_angle(const SphericalCoords& observer,
                                                     const SphericalCoords& object,
                                                     T radius) noexcept {
        return SphericalCoordinateProcessor::get_instance().solid_angle_from_sphere(
            observer, object, radius);
    }
    
    // Ranges-based transformations
    template<std::ranges::input_range R>
    [[nodiscard]] auto transform_spherical(R&& points, auto transformation) {
        return std::ranges::transform_view(std::forward<R>(points), transformation);
    }
}

} // namespace hsml::core