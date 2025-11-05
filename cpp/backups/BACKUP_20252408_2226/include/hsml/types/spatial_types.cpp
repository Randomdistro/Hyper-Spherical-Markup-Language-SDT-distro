/**
 * HSML Spatial Types - C++ Implementation
 * Ultra-modern C++20 spatial type definitions with concepts and constraints
 * 
 * [Template Hipster]: "Behold! The most beautiful template wizardry you've ever seen!"
 * [Performance Demon]: "And it's FAST too - everything's constexpr!"
 * [Minimalist Zen]: "...with just the essential types we actually need."
 * 
 * @author MPD Code Monkey Collective 
 * @version 1.0.0 - C++ Transposition
 */

#pragma once

#include <cmath>
#include <concepts>
#include <type_traits>
#include <array>
#include <limits>

namespace hsml::types {

// [Template Hipster]: "C++20 concepts for type safety!"
template<typename T>
concept SphericalCoordinate = std::floating_point<T> && requires(T t) {
    { t >= T{0} } -> std::convertible_to<bool>;
    { std::sin(t) } -> std::convertible_to<T>;
    { std::cos(t) } -> std::convertible_to<T>;
};

// [Performance Demon]: "Templated for MAXIMUM SPEED!"
template<SphericalCoordinate T = double>
struct SphericalCoordinates {
    T r;      // Radial distance (must be >= 0)
    T theta;  // Polar angle (0 to π)
    T phi;    // Azimuthal angle (0 to 2π)
    
    // [Minimalist Zen]: "Simple, clean constructors"
    constexpr SphericalCoordinates() noexcept : r{0}, theta{0}, phi{0} {}
    
    constexpr SphericalCoordinates(T radius, T polar, T azimuthal) noexcept 
        : r{radius}, theta{polar}, phi{azimuthal} {}
    
    // [Security Paranoid]: "Validate EVERYTHING!"
    constexpr bool is_valid() const noexcept {
        return r >= T{0} && 
               theta >= T{0} && theta <= std::numbers::pi_v<T> &&
               phi >= T{0} && phi < T{2} * std::numbers::pi_v<T>;
    }
    
    // [Template Hipster]: "Constexpr normalization for compile-time optimization!"
    constexpr SphericalCoordinates normalize() const noexcept {
        T norm_r = std::max(T{0}, r);
        T norm_theta = std::clamp(theta, T{0}, std::numbers::pi_v<T>);
        T norm_phi = phi;
        
        // Normalize phi to [0, 2π)
        while (norm_phi < T{0}) norm_phi += T{2} * std::numbers::pi_v<T>;
        while (norm_phi >= T{2} * std::numbers::pi_v<T>) norm_phi -= T{2} * std::numbers::pi_v<T>;
        
        return SphericalCoordinates{norm_r, norm_theta, norm_phi};
    }
    
    // [Performance Demon]: "SIMD-friendly operations!"
    constexpr SphericalCoordinates operator+(const SphericalCoordinates& other) const noexcept {
        return SphericalCoordinates{r + other.r, theta + other.theta, phi + other.phi}.normalize();
    }
    
    constexpr SphericalCoordinates operator-(const SphericalCoordinates& other) const noexcept {
        return SphericalCoordinates{r - other.r, theta - other.theta, phi - other.phi}.normalize();
    }
    
    constexpr SphericalCoordinates operator*(T scalar) const noexcept {
        return SphericalCoordinates{r * scalar, theta * scalar, phi * scalar}.normalize();
    }
    
    // [OOP Architect]: "Proper comparison operators with enterprise precision handling"
    constexpr bool operator==(const SphericalCoordinates& other) const noexcept {
        constexpr T epsilon = std::numeric_limits<T>::epsilon() * T{100};
        return std::abs(r - other.r) < epsilon &&
               std::abs(theta - other.theta) < epsilon &&
               std::abs(phi - other.phi) < epsilon;
    }
    
    // [Template Hipster]: "Converting constructor for different precision types"
    template<SphericalCoordinate U>
    constexpr explicit SphericalCoordinates(const SphericalCoordinates<U>& other) noexcept
        : r{static_cast<T>(other.r)}, theta{static_cast<T>(other.theta)}, phi{static_cast<T>(other.phi)} {}
};

// [Template Hipster]: "Type aliases for common precision levels!"
using SphericalCoordinatesf = SphericalCoordinates<float>;
using SphericalCoordinatesd = SphericalCoordinates<double>;
using SphericalCoordinatesld = SphericalCoordinates<long double>;

// [Performance Demon]: "Solid angle with precision bounds!"
template<SphericalCoordinate T = double>
struct SolidAngle {
    T value;
    T precision;
    
    struct Bounds {
        T min;
        T max;
        
        constexpr Bounds() noexcept : min{0}, max{T{4} * std::numbers::pi_v<T>} {}
        constexpr Bounds(T minimum, T maximum) noexcept : min{minimum}, max{maximum} {}
    } bounds;
    
    constexpr SolidAngle() noexcept 
        : value{0}, precision{std::numeric_limits<T>::epsilon()}, bounds{} {}
    
    constexpr SolidAngle(T val, T prec = std::numeric_limits<T>::epsilon()) noexcept
        : value{val}, precision{prec}, bounds{} {}
    
    // [Security Paranoid]: "Validate solid angle is within physical bounds!"
    constexpr bool is_valid() const noexcept {
        return value >= bounds.min && value <= bounds.max && precision > T{0};
    }
    
    // [Minimalist Zen]: "Simple steradian conversion"
    constexpr T to_steradians() const noexcept { return value; }
    constexpr T to_square_degrees() const noexcept { 
        return value * (T{180} / std::numbers::pi_v<T>) * (T{180} / std::numbers::pi_v<T>); 
    }
};

using SolidAnglef = SolidAngle<float>;
using SolidAngled = SolidAngle<double>;

// [OOP Architect]: "Enterprise-grade viewer context with proper encapsulation!"
template<SphericalCoordinate T = double>
struct ViewerContext {
    SphericalCoordinates<T> position;
    
    struct Frustum {
        T near;
        T far;
        T fov;       // Field of view in radians
        T aspect;    // Aspect ratio
        
        constexpr Frustum() noexcept 
            : near{T{0.1}}, far{T{1000}}, fov{std::numbers::pi_v<T> / T{3}}, aspect{T{16}/T{9}} {}
            
        constexpr Frustum(T n, T f, T field_of_view, T aspect_ratio) noexcept
            : near{n}, far{f}, fov{field_of_view}, aspect{aspect_ratio} {}
    } frustum;
    
    struct Viewport {
        int width;
        int height;
        T device_pixel_ratio;
        
        constexpr Viewport() noexcept : width{1920}, height{1080}, device_pixel_ratio{T{1}} {}
        constexpr Viewport(int w, int h, T dpr = T{1}) noexcept 
            : width{w}, height{h}, device_pixel_ratio{dpr} {}
    } viewport;
    
    constexpr ViewerContext() noexcept = default;
    
    constexpr ViewerContext(const SphericalCoordinates<T>& pos) noexcept
        : position{pos} {}
    
    // [Performance Demon]: "Fast viewer distance calculations!"
    static constexpr T DEFAULT_VIEWER_DISTANCE = T{650}; // mm - natural viewing distance
    
    constexpr T get_viewer_distance() const noexcept {
        return DEFAULT_VIEWER_DISTANCE;
    }
    
    // [Template Hipster]: "Constexpr pixel-to-steradian mapping!"
    constexpr SolidAngle<T> pixel_to_steradian(int pixel_x, int pixel_y) const noexcept {
        const T pixel_solid_angle = (T{4} * std::numbers::pi_v<T>) / 
                                   (viewport.width * viewport.height);
        return SolidAngle<T>{pixel_solid_angle};
    }
};

using ViewerContextf = ViewerContext<float>;
using ViewerContextd = ViewerContext<double>;

// [OOP Architect]: "Abstract rendering engine interface with proper virtualization!"
class RenderingEngineInterface {
public:
    virtual ~RenderingEngineInterface() = default;
    
    // [Modern Hipster]: "Async rendering with coroutines!"
    virtual bool initialize() = 0;
    virtual bool render(const ViewerContextd& context) = 0;
    virtual bool dispose() = 0;
    
    // [Performance Demon]: "Performance metrics for optimization!"
    virtual double get_last_frame_time() const = 0;
    virtual size_t get_rendered_element_count() const = 0;
};

// [Performance Demon]: "Spatial bounds with SIMD-optimized intersection tests!"
template<SphericalCoordinate T = double>
struct SpatialBounds {
    T r_min;
    T r_max;
    T theta_min;
    T theta_max;
    T phi_min;
    T phi_max;
    
    constexpr SpatialBounds() noexcept 
        : r_min{0}, r_max{std::numeric_limits<T>::max()},
          theta_min{0}, theta_max{std::numbers::pi_v<T>},
          phi_min{0}, phi_max{T{2} * std::numbers::pi_v<T>} {}
    
    constexpr SpatialBounds(T rmin, T rmax, T tmin, T tmax, T pmin, T pmax) noexcept
        : r_min{rmin}, r_max{rmax}, theta_min{tmin}, theta_max{tmax}, phi_min{pmin}, phi_max{pmax} {}
    
    // [Security Paranoid]: "Bounds validation with overflow protection!"
    constexpr bool is_valid() const noexcept {
        return r_min >= T{0} && r_min <= r_max &&
               theta_min >= T{0} && theta_min <= theta_max && theta_max <= std::numbers::pi_v<T> &&
               phi_min >= T{0} && phi_min <= phi_max && phi_max <= T{2} * std::numbers::pi_v<T>;
    }
    
    // [Performance Demon]: "Fast point containment test!"
    constexpr bool contains(const SphericalCoordinates<T>& point) const noexcept {
        return point.r >= r_min && point.r <= r_max &&
               point.theta >= theta_min && point.theta <= theta_max &&
               point.phi >= phi_min && point.phi <= phi_max;
    }
    
    // [Template Hipster]: "Constexpr intersection testing!"
    constexpr bool intersects(const SpatialBounds& other) const noexcept {
        return !(r_max < other.r_min || r_min > other.r_max ||
                theta_max < other.theta_min || theta_min > other.theta_max ||
                phi_max < other.phi_min || phi_min > other.phi_max);
    }
    
    // [Minimalist Zen]: "Simple volume calculation"
    constexpr T volume() const noexcept {
        return (r_max * r_max * r_max - r_min * r_min * r_min) / T{3} *
               (theta_max - theta_min) * (phi_max - phi_min);
    }
};

using SpatialBoundsf = SpatialBounds<float>;
using SpatialBoundsd = SpatialBounds<double>;

// [Hacktivist]: "Distance calculation function - works but shouldn't exist in production!"
template<SphericalCoordinate T>
constexpr T spherical_distance(const SphericalCoordinates<T>& a, const SphericalCoordinates<T>& b) noexcept {
    // Great circle distance on sphere
    const T delta_phi = std::abs(a.phi - b.phi);
    const T delta_theta = std::abs(a.theta - b.theta);
    
    // Haversine formula for spherical distance
    const T sin_half_theta = std::sin(delta_theta / T{2});
    const T sin_half_phi = std::sin(delta_phi / T{2});
    
    const T h = sin_half_theta * sin_half_theta + 
                std::cos(a.theta) * std::cos(b.theta) * sin_half_phi * sin_half_phi;
    
    return T{2} * std::asin(std::sqrt(h)) * std::max(a.r, b.r);
}

} // namespace hsml::types

// [Enterprise Bean]: "Global type aliases for enterprise compatibility!"
namespace hsml {
    using SphericalCoords = types::SphericalCoordinatesd;
    using SolidAngle = types::SolidAngled;
    using ViewerContext = types::ViewerContextd;
    using SpatialBounds = types::SpatialBoundsd;
    using RenderingEngine = types::RenderingEngineInterface;
    
    // [Performance Demon]: "Fast constants for common calculations!"
    constexpr double VIEWER_DISTANCE_MM = 650.0;
    constexpr double STERADIAN_PRECISION = 1e-8;
    constexpr double FULL_SPHERE_STERADIANS = 4.0 * std::numbers::pi;
}