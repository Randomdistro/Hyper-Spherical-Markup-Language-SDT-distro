#pragma once

#include "spherical_coords.h"
#include "vector3.h"
#include "matrix4.h"
#include "simd_math.h"
#include <concepts>
#include <coroutine>
#include <ranges>
#include <memory>
#include <unordered_map>
#include <array>
#include <span>
#include <functional>
#include <type_traits>
#include <execution>
#include <algorithm>
#include <immintrin.h>

namespace hsml {
namespace core {

// === COMPILE-TIME SPHERICAL MATHEMATICS (FUNCTIONAL PARADIGM) ===

template<typename T>
concept SphericalCoordinate = requires(T t) {
    t.radius();
    t.theta();
    t.phi();
};

template<typename T>
concept SolidAngleType = requires(T t) {
    t.omega();
    t.theta_min();
    t.theta_max();
    t.phi_min();
    t.phi_max();
};

// Compile-time solid angle calculations using constexpr templates
template<std::floating_point T>
struct SolidAngle {
    T omega_;
    T theta_min_, theta_max_;
    T phi_min_, phi_max_;
    
    constexpr SolidAngle(T omega, T th_min, T th_max, T ph_min, T ph_max)
        : omega_(omega), theta_min_(th_min), theta_max_(th_max),
          phi_min_(ph_min), phi_max_(ph_max) {}
    
    constexpr T omega() const { return omega_; }
    constexpr T theta_min() const { return theta_min_; }
    constexpr T theta_max() const { return theta_max_; }
    constexpr T phi_min() const { return phi_min_; }
    constexpr T phi_max() const { return phi_max_; }
    
    // Compile-time solid angle integration
    constexpr T integrated_omega() const {
        return (phi_max_ - phi_min_) * (std::cos(theta_min_) - std::cos(theta_max_));
    }
};

using SolidAnglef = SolidAngle<float>;
using SolidAngled = SolidAngle<double>;

// === TEMPLATE METAPROGRAMMING FOR PIXEL MAPPING ===

template<typename T, size_t Width, size_t Height>
struct PixelMappingLUT {
    static constexpr size_t size = Width * Height;
    std::array<SolidAngle<T>, size> mappings;
    
    constexpr PixelMappingLUT(T viewer_distance, T monitor_width, T monitor_height) {
        for (size_t y = 0; y < Height; ++y) {
            for (size_t x = 0; x < Width; ++x) {
                mappings[y * Width + x] = calculate_pixel_mapping(x, y, viewer_distance, 
                                                                monitor_width, monitor_height);
            }
        }
    }
    
    constexpr SolidAngle<T> get(size_t x, size_t y) const {
        return mappings[y * Width + x];
    }
    
private:
    constexpr SolidAngle<T> calculate_pixel_mapping(size_t x, size_t y, T viewer_distance,
                                                   T monitor_width, T monitor_height) const {
        T physical_x = (static_cast<T>(x) / Width) * monitor_width - (monitor_width / 2);
        T physical_y = (static_cast<T>(y) / Height) * monitor_height - (monitor_height / 2);
        
        T theta = std::atan2(physical_y, viewer_distance);
        T phi = std::atan2(physical_x, viewer_distance);
        
        T pixel_width = monitor_width / Width;
        T pixel_height = monitor_height / Height;
        T solid_angle_size = (pixel_width * pixel_height) / (viewer_distance * viewer_distance);
        
        return SolidAngle<T>(
            solid_angle_size,
            theta - (pixel_height / 2) / viewer_distance,
            theta + (pixel_height / 2) / viewer_distance,
            phi - (pixel_width / 2) / viewer_distance,
            phi + (pixel_width / 2) / viewer_distance
        );
    }
};

// === SIMD-OPTIMIZED STERADIAN-TO-PIXEL MAPPING ===

class SIMDSteradianProcessor {
public:
    struct alignas(32) SIMDSolidAngle {
        __m256 omega;
        __m256 theta_min, theta_max;
        __m256 phi_min, phi_max;
    };
    
    // Process 8 pixels simultaneously using AVX2
    static void process_pixel_batch_avx2(
        const float* pixel_coords,  // [x0,y0,x1,y1,...,x7,y7]
        float viewer_distance,
        float monitor_width, float monitor_height,
        SIMDSolidAngle* result
    ) {
        const __m256 vd = _mm256_set1_ps(viewer_distance);
        const __m256 mw = _mm256_set1_ps(monitor_width);
        const __m256 mh = _mm256_set1_ps(monitor_height);
        const __m256 half_mw = _mm256_set1_ps(monitor_width * 0.5f);
        const __m256 half_mh = _mm256_set1_ps(monitor_height * 0.5f);
        
        // Load x and y coordinates
        __m256 x_coords = _mm256_setr_ps(
            pixel_coords[0], pixel_coords[2], pixel_coords[4], pixel_coords[6],
            pixel_coords[8], pixel_coords[10], pixel_coords[12], pixel_coords[14]
        );
        __m256 y_coords = _mm256_setr_ps(
            pixel_coords[1], pixel_coords[3], pixel_coords[5], pixel_coords[7],
            pixel_coords[9], pixel_coords[11], pixel_coords[13], pixel_coords[15]
        );
        
        // Convert to physical coordinates
        __m256 physical_x = _mm256_sub_ps(_mm256_mul_ps(_mm256_div_ps(x_coords, mw), mw), half_mw);
        __m256 physical_y = _mm256_sub_ps(_mm256_mul_ps(_mm256_div_ps(y_coords, mh), mh), half_mh);
        
        // Calculate angles (using approximation for SIMD)
        __m256 theta = simd_atan2_ps(physical_y, vd);
        __m256 phi = simd_atan2_ps(physical_x, vd);
        
        // Calculate solid angle size
        __m256 pixel_area = _mm256_div_ps(mw, _mm256_set1_ps(8.0f)) * _mm256_div_ps(mh, _mm256_set1_ps(8.0f));
        __m256 distance_sq = _mm256_mul_ps(vd, vd);
        result->omega = _mm256_div_ps(pixel_area, distance_sq);
        
        // Calculate bounds
        __m256 half_pixel_size = _mm256_set1_ps(0.001f); // Approximate
        result->theta_min = _mm256_sub_ps(theta, half_pixel_size);
        result->theta_max = _mm256_add_ps(theta, half_pixel_size);
        result->phi_min = _mm256_sub_ps(phi, half_pixel_size);
        result->phi_max = _mm256_add_ps(phi, half_pixel_size);
    }
    
private:
    // Fast SIMD atan2 approximation
    static __m256 simd_atan2_ps(__m256 y, __m256 x) {
        // Simplified atan2 approximation for SIMD
        // Using polynomial approximation for better performance
        __m256 ratio = _mm256_div_ps(y, x);
        __m256 abs_ratio = _mm256_andnot_ps(_mm256_set1_ps(-0.0f), ratio);
        
        // Polynomial coefficients for atan approximation
        const __m256 a1 = _mm256_set1_ps(0.999866f);
        const __m256 a3 = _mm256_set1_ps(-0.3302995f);
        const __m256 a5 = _mm256_set1_ps(0.180141f);
        const __m256 a7 = _mm256_set1_ps(-0.085133f);
        
        __m256 r2 = _mm256_mul_ps(ratio, ratio);
        __m256 r4 = _mm256_mul_ps(r2, r2);
        __m256 r6 = _mm256_mul_ps(r4, r2);
        
        __m256 result = _mm256_add_ps(
            _mm256_mul_ps(a1, ratio),
            _mm256_add_ps(
                _mm256_mul_ps(a3, _mm256_mul_ps(ratio, r2)),
                _mm256_add_ps(
                    _mm256_mul_ps(a5, _mm256_mul_ps(ratio, r4)),
                    _mm256_mul_ps(a7, _mm256_mul_ps(ratio, r6))
                )
            )
        );
        
        return result;
    }
};

// === COROUTINE-BASED ASYNC SPHERICAL CALCULATIONS ===

struct SphericalTask {
    struct promise_type {
        SphericalCoords result;
        
        SphericalTask get_return_object() {
            return SphericalTask{std::coroutine_handle<promise_type>::from_promise(*this)};
        }
        
        std::suspend_never initial_suspend() { return {}; }
        std::suspend_always final_suspend() noexcept { return {}; }
        void unhandled_exception() {}
        
        void return_value(const SphericalCoords& value) {
            result = value;
        }
    };
    
    std::coroutine_handle<promise_type> coro;
    
    SphericalTask(std::coroutine_handle<promise_type> h) : coro(h) {}
    ~SphericalTask() { if (coro) coro.destroy(); }
    
    SphericalCoords get_result() {
        if (coro && coro.done()) {
            return coro.promise().result;
        }
        return SphericalCoords{};
    }
    
    bool is_ready() const {
        return coro && coro.done();
    }
};

// Async spherical distance calculation
SphericalTask async_spherical_distance(const SphericalCoords& p1, const SphericalCoords& p2) {
    // Simulate async work
    co_await std::suspend_always{};
    
    // Pure spherical distance calculation
    double cos_angular = std::cos(p1.theta()) * std::cos(p2.theta()) + 
                        std::sin(p1.theta()) * std::sin(p2.theta()) * 
                        std::cos(p2.phi() - p1.phi());
    
    double angular_distance = std::acos(std::clamp(cos_angular, -1.0, 1.0));
    double distance = std::sqrt(
        p1.radius() * p1.radius() + p2.radius() * p2.radius() - 
        2 * p1.radius() * p2.radius() * cos_angular
    );
    
    co_return SphericalCoords(distance, angular_distance, 0.0);
}

// === LOCK-FREE CONCURRENT DOM ELEMENT CACHE ===

template<typename Key, typename Value>
class LockFreeSpatialCache {
public:
    struct Node {
        std::atomic<Key> key;
        std::atomic<Value*> value;
        std::atomic<Node*> next;
        
        Node() : key{}, value{nullptr}, next{nullptr} {}
    };
    
private:
    static constexpr size_t TABLE_SIZE = 1024;
    std::array<std::atomic<Node*>, TABLE_SIZE> table_;
    std::atomic<size_t> size_{0};
    
    size_t hash(const Key& key) const {
        return std::hash<Key>{}(key) % TABLE_SIZE;
    }
    
public:
    LockFreeSpatialCache() {
        for (auto& bucket : table_) {
            bucket.store(nullptr, std::memory_order_relaxed);
        }
    }
    
    bool insert(const Key& key, std::unique_ptr<Value> value) {
        size_t index = hash(key);
        auto new_node = std::make_unique<Node>();
        new_node->key.store(key, std::memory_order_relaxed);
        new_node->value.store(value.release(), std::memory_order_relaxed);
        
        Node* expected = table_[index].load(std::memory_order_acquire);
        new_node->next.store(expected, std::memory_order_relaxed);
        
        while (!table_[index].compare_exchange_weak(
            expected, new_node.get(), 
            std::memory_order_release, std::memory_order_acquire)) {
            new_node->next.store(expected, std::memory_order_relaxed);
        }
        
        new_node.release();
        size_.fetch_add(1, std::memory_order_relaxed);
        return true;
    }
    
    Value* find(const Key& key) const {
        size_t index = hash(key);
        Node* current = table_[index].load(std::memory_order_acquire);
        
        while (current) {
            if (current->key.load(std::memory_order_relaxed) == key) {
                return current->value.load(std::memory_order_relaxed);
            }
            current = current->next.load(std::memory_order_acquire);
        }
        
        return nullptr;
    }
    
    size_t size() const {
        return size_.load(std::memory_order_relaxed);
    }
};

// === EXPRESSION TEMPLATES FOR SPHERICAL OPERATIONS ===

template<typename E>
struct SphericalExpression {
    const E& derived() const { return static_cast<const E&>(*this); }
    E& derived() { return static_cast<E&>(*this); }
};

template<typename LHS, typename RHS>
struct SphericalAdd : SphericalExpression<SphericalAdd<LHS, RHS>> {
    const LHS& lhs;
    const RHS& rhs;
    
    SphericalAdd(const LHS& l, const RHS& r) : lhs(l), rhs(r) {}
    
    SphericalCoords eval() const {
        auto left = lhs.eval();
        auto right = rhs.eval();
        // Spherical addition via Cartesian intermediate
        auto cart_sum = left.to_cartesian() + right.to_cartesian();
        return SphericalCoords::from_cartesian(cart_sum);
    }
};

template<typename Expr>
struct SphericalScale : SphericalExpression<SphericalScale<Expr>> {
    const Expr& expr;
    double scale;
    
    SphericalScale(const Expr& e, double s) : expr(e), scale(s) {}
    
    SphericalCoords eval() const {
        auto result = expr.eval();
        return SphericalCoords(result.radius() * scale, result.theta(), result.phi());
    }
};

// Make SphericalCoords work with expression templates
template<>
struct SphericalExpression<SphericalCoords> {
    const SphericalCoords& coord;
    
    SphericalExpression(const SphericalCoords& c) : coord(c) {}
    
    SphericalCoords eval() const { return coord; }
};

// Expression template operators
template<typename LHS, typename RHS>
auto operator+(const SphericalExpression<LHS>& lhs, const SphericalExpression<RHS>& rhs) {
    return SphericalAdd<LHS, RHS>(lhs.derived(), rhs.derived());
}

template<typename Expr>
auto operator*(const SphericalExpression<Expr>& expr, double scale) {
    return SphericalScale<Expr>(expr.derived(), scale);
}

// === ADVANCED SOLID ANGLE DOM PROCESSOR ===

class AdvancedSolidAngleDOMProcessor {
public:
    struct HSMLElement {
        std::string id;
        std::string tag_name;
        SphericalCoords spherical_position;
        SolidAngled solid_angle_size;
        std::vector<std::shared_ptr<HSMLElement>> children;
        std::weak_ptr<HSMLElement> parent;
        bool visible = true;
        bool interactive = true;
        std::string matter_state = "solid";
        std::unordered_map<std::string, std::any> material_properties;
        
        // Cache invalidation
        mutable std::atomic<uint64_t> cache_version{0};
    };
    
    struct HSMLViewport {
        uint32_t width, height;
        double viewer_distance;
        double monitor_width, monitor_height;
        double total_solid_angle;
        
        // 4-corner optimization
        SolidAngled top_left, top_right, bottom_left, bottom_right;
        
        // Field of view
        double fov_horizontal, fov_vertical, fov_diagonal;
    };
    
    struct RaycastResult {
        std::shared_ptr<HSMLElement> element;
        double distance;
        SphericalCoords intersection_point;
        SolidAngled solid_angle;
    };
    
private:
    std::unique_ptr<HSMLViewport> viewport_;
    std::shared_ptr<HSMLElement> root_element_;
    std::unordered_map<std::string, std::shared_ptr<HSMLElement>> element_registry_;
    
    // Lock-free caching systems
    LockFreeSpatialCache<std::pair<uint32_t, uint32_t>, RaycastResult> raycast_cache_;
    LockFreeSpatialCache<std::pair<uint32_t, uint32_t>, SolidAngled> pixel_mapping_cache_;
    
    // Performance tracking
    std::atomic<uint64_t> frame_count_{0};
    std::atomic<uint64_t> cache_hits_{0};
    std::atomic<uint64_t> cache_misses_{0};
    
    // Compile-time LUT for common resolutions
    static constexpr PixelMappingLUT<double, 1920, 1080> HD_LUT{650.0, 531.0, 299.0};
    static constexpr PixelMappingLUT<double, 3840, 2160> UHD_LUT{650.0, 531.0, 299.0};
    
public:
    AdvancedSolidAngleDOMProcessor() = default;
    
    // Initialize viewport with compile-time optimization
    template<uint32_t Width, uint32_t Height>
    void initialize_viewport(double viewer_distance = 650.0,
                           double monitor_width = 531.0,
                           double monitor_height = 299.0) {
        viewport_ = std::make_unique<HSMLViewport>();
        viewport_->width = Width;
        viewport_->height = Height;
        viewport_->viewer_distance = viewer_distance;
        viewport_->monitor_width = monitor_width;
        viewport_->monitor_height = monitor_height;
        
        // Calculate field of view
        viewport_->fov_horizontal = 2.0 * std::atan(monitor_width / (2.0 * viewer_distance));
        viewport_->fov_vertical = 2.0 * std::atan(monitor_height / (2.0 * viewer_distance));
        viewport_->fov_diagonal = 2.0 * std::atan(
            std::sqrt(monitor_width * monitor_width + monitor_height * monitor_height) / 
            (2.0 * viewer_distance)
        );
        
        // Calculate total solid angle
        viewport_->total_solid_angle = 
            (viewport_->fov_horizontal) * 
            (std::cos(-viewport_->fov_vertical / 2.0) - std::cos(viewport_->fov_vertical / 2.0));
        
        // Initialize 4-corner optimization
        initialize_corner_optimization<Width, Height>();
    }
    
    // SIMD-optimized pixel-to-solid-angle mapping
    SolidAngled get_pixel_solid_angle_simd(uint32_t pixel_x, uint32_t pixel_y) const {
        if (!viewport_) return SolidAngled{0, 0, 0, 0, 0};
        
        // Check cache first
        auto cache_key = std::make_pair(pixel_x, pixel_y);
        if (auto cached = pixel_mapping_cache_.find(cache_key)) {
            cache_hits_.fetch_add(1, std::memory_order_relaxed);
            return *cached;
        }
        
        // Use 4-corner bilinear interpolation for optimization
        double u = static_cast<double>(pixel_x) / (viewport_->width - 1);
        double v = static_cast<double>(pixel_y) / (viewport_->height - 1);
        
        // Bilinear interpolation of solid angles
        auto interpolate_solid_angle = [](const SolidAngled& tl, const SolidAngled& tr,
                                        const SolidAngled& bl, const SolidAngled& br,
                                        double u, double v) {
            auto omega_top = tl.omega() * (1.0 - u) + tr.omega() * u;
            auto omega_bottom = bl.omega() * (1.0 - u) + br.omega() * u;
            auto omega = omega_top * (1.0 - v) + omega_bottom * v;
            
            // Similar interpolation for bounds
            auto theta_min = (tl.theta_min() * (1.0 - u) + tr.theta_min() * u) * (1.0 - v) +
                           (bl.theta_min() * (1.0 - u) + br.theta_min() * u) * v;
            auto theta_max = (tl.theta_max() * (1.0 - u) + tr.theta_max() * u) * (1.0 - v) +
                           (bl.theta_max() * (1.0 - u) + br.theta_max() * u) * v;
            auto phi_min = (tl.phi_min() * (1.0 - u) + tr.phi_min() * u) * (1.0 - v) +
                         (bl.phi_min() * (1.0 - u) + br.phi_min() * u) * v;
            auto phi_max = (tl.phi_max() * (1.0 - u) + tr.phi_max() * u) * (1.0 - v) +
                         (bl.phi_max() * (1.0 - u) + br.phi_max() * u) * v;
            
            return SolidAngled{omega, theta_min, theta_max, phi_min, phi_max};
        };
        
        auto result = interpolate_solid_angle(
            viewport_->top_left, viewport_->top_right,
            viewport_->bottom_left, viewport_->bottom_right,
            u, v
        );
        
        // Cache the result
        auto result_ptr = std::make_unique<SolidAngled>(result);
        pixel_mapping_cache_.insert(cache_key, std::move(result_ptr));
        cache_misses_.fetch_add(1, std::memory_order_relaxed);
        
        return result;
    }
    
    // Parallel raycast using std::execution
    std::optional<RaycastResult> perform_parallel_raycast(uint32_t pixel_x, uint32_t pixel_y) {
        if (!viewport_ || !root_element_) return std::nullopt;
        
        // Get solid angle for pixel
        auto solid_angle = get_pixel_solid_angle_simd(pixel_x, pixel_y);
        
        // Create ray direction from solid angle center
        double center_theta = (solid_angle.theta_min() + solid_angle.theta_max()) / 2.0;
        double center_phi = (solid_angle.phi_min() + solid_angle.phi_max()) / 2.0;
        SphericalCoords ray_direction{1.0, center_theta, center_phi};
        
        // Collect all visible elements
        std::vector<std::shared_ptr<HSMLElement>> visible_elements;
        collect_visible_elements(root_element_, visible_elements);
        
        // Parallel intersection testing
        std::vector<RaycastResult> results(visible_elements.size());
        
        std::transform(std::execution::par_unseq,
                      visible_elements.begin(), visible_elements.end(),
                      results.begin(),
                      [&ray_direction, &solid_angle](const auto& element) -> RaycastResult {
                          return test_ray_element_intersection(ray_direction, element, solid_angle);
                      });
        
        // Find closest intersection
        auto closest = std::min_element(results.begin(), results.end(),
                                      [](const auto& a, const auto& b) {
                                          return a.distance < b.distance && a.element != nullptr;
                                      });
        
        if (closest != results.end() && closest->element) {
            return *closest;
        }
        
        return std::nullopt;
    }
    
    // Create element with modern C++ features
    template<typename... Args>
    std::shared_ptr<HSMLElement> create_element(const std::string& tag_name, Args&&... args) {
        auto element = std::make_shared<HSMLElement>();
        element->id = generate_element_id();
        element->tag_name = tag_name;
        
        // Perfect forwarding for initialization
        (initialize_property(element, std::forward<Args>(args)), ...);
        
        element_registry_[element->id] = element;
        return element;
    }
    
    // Performance statistics
    struct PerformanceStats {
        uint64_t frame_count;
        uint64_t cache_hits;
        uint64_t cache_misses;
        double cache_hit_rate;
        size_t element_count;
        size_t raycast_cache_size;
        size_t pixel_cache_size;
    };
    
    PerformanceStats get_performance_stats() const {
        auto hits = cache_hits_.load(std::memory_order_relaxed);
        auto misses = cache_misses_.load(std::memory_order_relaxed);
        auto total = hits + misses;
        
        return PerformanceStats{
            frame_count_.load(std::memory_order_relaxed),
            hits,
            misses,
            total > 0 ? static_cast<double>(hits) / total : 0.0,
            element_registry_.size(),
            raycast_cache_.size(),
            pixel_mapping_cache_.size()
        };
    }
    
private:
    template<uint32_t Width, uint32_t Height>
    void initialize_corner_optimization() {
        if constexpr (Width == 1920 && Height == 1080) {
            viewport_->top_left = HD_LUT.get(0, 0);
            viewport_->top_right = HD_LUT.get(Width - 1, 0);
            viewport_->bottom_left = HD_LUT.get(0, Height - 1);
            viewport_->bottom_right = HD_LUT.get(Width - 1, Height - 1);
        } else if constexpr (Width == 3840 && Height == 2160) {
            viewport_->top_left = UHD_LUT.get(0, 0);
            viewport_->top_right = UHD_LUT.get(Width - 1, 0);
            viewport_->bottom_left = UHD_LUT.get(0, Height - 1);
            viewport_->bottom_right = UHD_LUT.get(Width - 1, Height - 1);
        } else {
            // Runtime calculation for non-standard resolutions
            viewport_->top_left = calculate_pixel_solid_angle(0, 0);
            viewport_->top_right = calculate_pixel_solid_angle(Width - 1, 0);
            viewport_->bottom_left = calculate_pixel_solid_angle(0, Height - 1);
            viewport_->bottom_right = calculate_pixel_solid_angle(Width - 1, Height - 1);
        }
    }
    
    SolidAngled calculate_pixel_solid_angle(uint32_t x, uint32_t y) const {
        if (!viewport_) return SolidAngled{0, 0, 0, 0, 0};
        
        double physical_x = (static_cast<double>(x) / viewport_->width) * viewport_->monitor_width - 
                           (viewport_->monitor_width / 2.0);
        double physical_y = (static_cast<double>(y) / viewport_->height) * viewport_->monitor_height - 
                           (viewport_->monitor_height / 2.0);
        
        double theta = std::atan2(physical_y, viewport_->viewer_distance);
        double phi = std::atan2(physical_x, viewport_->viewer_distance);
        
        double pixel_width = viewport_->monitor_width / viewport_->width;
        double pixel_height = viewport_->monitor_height / viewport_->height;
        double solid_angle_size = (pixel_width * pixel_height) / 
                                 (viewport_->viewer_distance * viewport_->viewer_distance);
        
        return SolidAngled{
            solid_angle_size,
            theta - (pixel_height / 2.0) / viewport_->viewer_distance,
            theta + (pixel_height / 2.0) / viewport_->viewer_distance,
            phi - (pixel_width / 2.0) / viewport_->viewer_distance,
            phi + (pixel_width / 2.0) / viewport_->viewer_distance
        };
    }
    
    void collect_visible_elements(const std::shared_ptr<HSMLElement>& element,
                                 std::vector<std::shared_ptr<HSMLElement>>& result) const {
        if (!element || !element->visible) return;
        
        result.push_back(element);
        
        for (const auto& child : element->children) {
            collect_visible_elements(child, result);
        }
    }
    
    static RaycastResult test_ray_element_intersection(
        const SphericalCoords& ray,
        const std::shared_ptr<HSMLElement>& element,
        const SolidAngled& solid_angle
    ) {
        if (!element || !element->visible || !element->interactive) {
            return RaycastResult{nullptr, std::numeric_limits<double>::infinity(), {}, solid_angle};
        }
        
        // Pure spherical intersection test
        double spherical_distance = ray.angular_distance(element->spherical_position);
        double element_radius = 10.0; // Default element size
        
        if (spherical_distance <= element_radius) {
            double radial_distance = std::abs(ray.radius() - element->spherical_position.radius());
            double total_distance = std::sqrt(
                radial_distance * radial_distance + 
                (ray.radius() * spherical_distance) * (ray.radius() * spherical_distance)
            );
            
            if (total_distance <= element_radius) {
                return RaycastResult{
                    element,
                    total_distance,
                    element->spherical_position,
                    solid_angle
                };
            }
        }
        
        return RaycastResult{nullptr, std::numeric_limits<double>::infinity(), {}, solid_angle};
    }
    
    std::string generate_element_id() const {
        static std::atomic<uint64_t> counter{0};
        return "hsml_" + std::to_string(counter.fetch_add(1, std::memory_order_relaxed));
    }
    
    template<typename T>
    void initialize_property(std::shared_ptr<HSMLElement>& element, T&& property) {
        // Property initialization using perfect forwarding
        if constexpr (std::is_same_v<std::decay_t<T>, SphericalCoords>) {
            element->spherical_position = std::forward<T>(property);
        } else if constexpr (std::is_same_v<std::decay_t<T>, SolidAngled>) {
            element->solid_angle_size = std::forward<T>(property);
        } else if constexpr (std::is_same_v<std::decay_t<T>, std::string>) {
            element->matter_state = std::forward<T>(property);
        }
        // Add more property types as needed
    }
};

} // namespace core
} // namespace hsml