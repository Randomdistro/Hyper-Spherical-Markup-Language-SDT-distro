#pragma once

#include <cmath>
#include <type_traits>
#include <concepts>
#include <limits>

namespace hsml {
namespace core {
namespace math {

// ============================================================================
// PRECISION LEVEL SYSTEM - Compile-time precision configuration
// ============================================================================

enum class PrecisionLevel {
    Quantum,     // 1e-15 - Ultra-precise quantum calculations
    Scientific,  // 1e-12 - High-precision scientific computing  
    Standard,    // 1e-10 - Standard mathematical operations
    Engineering, // 1e-9  - Engineering calculations
    Graphics,    // 1e-6  - Graphics and approximations
    Display      // 1e-3  - User interface and display
};

template<PrecisionLevel Level>
struct precision_traits;

template<> struct precision_traits<PrecisionLevel::Quantum> {
    static constexpr double epsilon = 1e-15;
    static constexpr double safe_divisor = 1e-15;
    static constexpr double origin_marker = 1e-15;
};

template<> struct precision_traits<PrecisionLevel::Scientific> {
    static constexpr double epsilon = 1e-12;
    static constexpr double safe_divisor = 1e-12;
    static constexpr double origin_marker = 1e-12;
};

template<> struct precision_traits<PrecisionLevel::Standard> {
    static constexpr double epsilon = 1e-10;
    static constexpr double safe_divisor = 1e-10;
    static constexpr double origin_marker = 1e-10;
};

template<> struct precision_traits<PrecisionLevel::Engineering> {
    static constexpr double epsilon = 1e-9;
    static constexpr double safe_divisor = 1e-9;
    static constexpr double origin_marker = 1e-9;
};

template<> struct precision_traits<PrecisionLevel::Graphics> {
    static constexpr double epsilon = 1e-6;
    static constexpr double safe_divisor = 1e-6;
    static constexpr double origin_marker = 1e-6;
};

template<> struct precision_traits<PrecisionLevel::Display> {
    static constexpr double epsilon = 1e-3;
    static constexpr double safe_divisor = 1e-3;
    static constexpr double origin_marker = 1e-3;
};

// ============================================================================
// MATHEMATICAL CONSTANTS - Ultimate precision with compile-time validation
// ============================================================================

struct constants {
    static constexpr double pi = 3.14159265358979323846264338327950288;
    static constexpr double two_pi = 2.0 * pi;
    static constexpr double half_pi = pi / 2.0;
    static constexpr double e = 2.71828182845904523536028747135266250;
    static constexpr double sqrt_2 = 1.41421356237309504880168872420969808;
    static constexpr double sqrt_3 = 1.73205080756887729352744634150587237;
    
    // Physical constants for SDT framework
    static constexpr double c = 299792458.0;        // Speed of light (m/s)
    static constexpr double h = 6.62607015e-34;     // Planck constant (J⋅s)
    static constexpr double epsilon0 = 8.8541878128e-12; // Vacuum permittivity (F/m)
    static constexpr double mu0 = 1.25663706212e-6;     // Vacuum permeability (H/m)
};

// ============================================================================
// CONCEPTS - Type safety for mathematical operations
// ============================================================================

template<typename T>
concept Numeric = std::is_arithmetic_v<T>;

template<typename T>
concept FloatingPoint = std::is_floating_point_v<T>;

template<typename T>
concept Vector3Like = requires(T v) {
    { v.x() } -> std::convertible_to<double>;
    { v.y() } -> std::convertible_to<double>;
    { v.z() } -> std::convertible_to<double>;
};

template<typename T>
concept SphericalLike = requires(T s) {
    { s.radius() } -> std::convertible_to<double>;
    { s.theta() } -> std::convertible_to<double>;
    { s.phi() } -> std::convertible_to<double>;
};

// ============================================================================
// SAFE MATHEMATICAL OPERATIONS - Unified error handling
// ============================================================================

template<PrecisionLevel Level = PrecisionLevel::Standard>
class SafeMath {
private:
    using traits = precision_traits<Level>;
    
public:
    // Safe division with compile-time precision selection
    template<FloatingPoint T>
    static constexpr T safe_divide(T numerator, T denominator) noexcept {
        const T safe_denom = (std::abs(denominator) < traits::safe_divisor) 
                           ? traits::safe_divisor 
                           : denominator;
        return numerator / safe_denom;
    }
    
    // Safe square root with negative protection
    template<FloatingPoint T>
    static constexpr T safe_sqrt(T value) noexcept {
        return std::sqrt(std::max(T{0}, value));
    }
    
    // Safe arc cosine with domain clamping
    template<FloatingPoint T>
    static constexpr T safe_acos(T value) noexcept {
        const T clamped = std::max(T{-1}, std::min(T{1}, value));
        return std::acos(clamped);
    }
    
    // Safe magnitude calculation for 3D vectors
    template<FloatingPoint T>
    static constexpr T safe_magnitude(T x, T y, T z) noexcept {
        return safe_sqrt(x*x + y*y + z*z);
    }
    
    // Precision-aware zero checking
    template<FloatingPoint T>
    static constexpr bool is_effectively_zero(T value) noexcept {
        return std::abs(value) < traits::epsilon;
    }
    
    // Precision-aware equality checking
    template<FloatingPoint T>
    static constexpr bool are_equal(T a, T b) noexcept {
        return is_effectively_zero(a - b);
    }
    
    // Origin marker for zero-degree elimination
    template<FloatingPoint T>
    static constexpr T origin_marker() noexcept {
        return static_cast<T>(traits::origin_marker);
    }
    
    // Safe normalization that returns origin marker instead of zero
    template<Vector3Like Vec>
    static Vec safe_normalize(const Vec& vector) {
        auto mag = safe_magnitude(vector.x(), vector.y(), vector.z());
        if (is_effectively_zero(mag)) {
            return Vec(origin_marker<double>(), origin_marker<double>(), origin_marker<double>());
        }
        const auto inv_mag = safe_divide(1.0, mag);
        return Vec(vector.x() * inv_mag, vector.y() * inv_mag, vector.z() * inv_mag);
    }
};

// ============================================================================
// VALIDATED VALUE WRAPPER - Debug-time validation with zero runtime cost
// ============================================================================

template<typename T, PrecisionLevel Level = PrecisionLevel::Standard>
class ValidatedValue {
private:
    T value_;
    using math = SafeMath<Level>;
    
    static constexpr bool is_debug_build() noexcept {
#ifdef NDEBUG
        return false;
#else
        return true;
#endif
    }
    
    constexpr void validate_finite() const {
        if constexpr (is_debug_build() && std::is_floating_point_v<T>) {
            if (!std::isfinite(value_)) {
                throw std::domain_error("ValidatedValue: non-finite value detected");
            }
        }
    }
    
    constexpr void validate_range(T min_val, T max_val) const {
        if constexpr (is_debug_build()) {
            if (value_ < min_val || value_ > max_val) {
                throw std::out_of_range("ValidatedValue: value outside valid range");
            }
        }
    }

public:
    constexpr ValidatedValue() noexcept : value_{} {}
    
    constexpr ValidatedValue(T val) : value_(val) {
        validate_finite();
    }
    
    constexpr ValidatedValue(T val, T min_val, T max_val) : value_(val) {
        validate_finite();
        validate_range(min_val, max_val);
    }
    
    constexpr operator T() const noexcept { return value_; }
    constexpr T get() const noexcept { return value_; }
    
    constexpr ValidatedValue& operator=(T val) {
        value_ = val;
        validate_finite();
        return *this;
    }
    
    // Mathematical operations with validation
    constexpr ValidatedValue operator+(const ValidatedValue& other) const {
        return ValidatedValue(value_ + other.value_);
    }
    
    constexpr ValidatedValue operator-(const ValidatedValue& other) const {
        return ValidatedValue(value_ - other.value_);
    }
    
    constexpr ValidatedValue operator*(const ValidatedValue& other) const {
        return ValidatedValue(value_ * other.value_);
    }
    
    constexpr ValidatedValue operator/(const ValidatedValue& other) const {
        return ValidatedValue(math::safe_divide(value_, other.value_));
    }
    
    constexpr bool approximately_equals(const ValidatedValue& other) const noexcept {
        return math::are_equal(value_, other.value_);
    }
};

// ============================================================================
// TYPE ALIASES FOR COMMON PRECISION LEVELS
// ============================================================================

using QuantumMath = SafeMath<PrecisionLevel::Quantum>;
using ScientificMath = SafeMath<PrecisionLevel::Scientific>;
using StandardMath = SafeMath<PrecisionLevel::Standard>;
using EngineeringMath = SafeMath<PrecisionLevel::Engineering>;
using GraphicsMath = SafeMath<PrecisionLevel::Graphics>;
using DisplayMath = SafeMath<PrecisionLevel::Display>;

template<typename T>
using QuantumValue = ValidatedValue<T, PrecisionLevel::Quantum>;

template<typename T>
using ScientificValue = ValidatedValue<T, PrecisionLevel::Scientific>;

template<typename T> 
using StandardValue = ValidatedValue<T, PrecisionLevel::Standard>;

template<typename T>
using EngineeringValue = ValidatedValue<T, PrecisionLevel::Engineering>;

template<typename T>
using GraphicsValue = ValidatedValue<T, PrecisionLevel::Graphics>;

template<typename T>
using DisplayValue = ValidatedValue<T, PrecisionLevel::Display>;

} // namespace math
} // namespace core
} // namespace hsml