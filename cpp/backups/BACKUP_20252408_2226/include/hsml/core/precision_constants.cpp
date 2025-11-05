#pragma once

// [ASPIE ARCHITECT]: Centralized precision constants for mathematical purity
// Eliminates all epsilon inconsistencies throughout the codebase

#include <cstddef>
#include <cmath>
#include <algorithm>

namespace hsml {
namespace core {
namespace precision {

// [ASPIE ARCHITECT]: Standardized mathematical constants - NO MORE INCONSISTENCIES
class MathematicalConstants {
public:
    // Core mathematical constants with ultimate precision
    static constexpr double PI = 3.14159265358979323846264338327950288;
    static constexpr double TWO_PI = 2.0 * PI;
    static constexpr double HALF_PI = PI / 2.0;
    static constexpr double E = 2.71828182845904523536028747135266250;
    static constexpr double SQRT_2 = 1.41421356237309504880168872420969808;
    static constexpr double SQRT_3 = 1.73205080756887729352744634150587237;
    
    // Physical constants for SDT framework
    static constexpr double SPEED_OF_LIGHT = 299792458.0; // m/s (exact)
    static constexpr double PLANCK_CONSTANT = 6.62607015e-34; // J⋅s (exact)
    static constexpr double EPSILON_0 = 8.8541878128e-12; // F/m (2018 CODATA)
    static constexpr double MU_0 = 1.25663706212e-6; // H/m (exact)
};

// [ASPIE ARCHITECT]: Hierarchical precision levels for different computational needs
class PrecisionLevels {
public:
    // Origin markers - replaces ALL zero-degree handling
    static constexpr double ORIGIN_MARKER = 1e-10;           // Standard non-zero marker
    static constexpr double QUANTUM_ORIGIN = 1e-15;          // Ultra-precise quantum calculations
    
    // Standard computational epsilons
    static constexpr double ULTRA_PRECISE = 1e-15;          // Quantum mechanics, high precision
    static constexpr double HIGH_PRECISE = 1e-12;           // Advanced calculations
    static constexpr double STANDARD_PRECISE = 1e-10;       // Standard mathematical operations
    static constexpr double MODERATE_PRECISE = 1e-9;        // General purpose calculations
    static constexpr double LOW_PRECISE = 1e-6;             // Graphics, approximations
    static constexpr double VISUAL_PRECISE = 1e-3;          // User interface, display
    
    // Distance calculation tolerances
    static constexpr double DISTANCE_TOLERANCE = STANDARD_PRECISE;
    static constexpr double ANGULAR_TOLERANCE = STANDARD_PRECISE;
    static constexpr double SPHERICAL_TOLERANCE = STANDARD_PRECISE;
    
    // Numerical stability guards
    static constexpr double MIN_SAFE_RADIUS = ORIGIN_MARKER;
    static constexpr double MIN_SAFE_MAGNITUDE = ORIGIN_MARKER;
    static constexpr double MIN_SAFE_DIVISOR = ORIGIN_MARKER;
    static constexpr double MAX_SAFE_VALUE = 1e12;
    
    // Cache and performance thresholds
    static constexpr double CACHE_PRECISION = HIGH_PRECISE;
    static constexpr size_t MAX_CACHE_SIZE = 10000;
    
    // Thread safety constants
    static constexpr double THREAD_SAFE_EPSILON = STANDARD_PRECISE;
    
    // GUI and rendering tolerances
    static constexpr double GUI_TOLERANCE = LOW_PRECISE;
    static constexpr double RENDER_TOLERANCE = MODERATE_PRECISE;
    static constexpr double PIXEL_TOLERANCE = VISUAL_PRECISE;
};

// [ASPIE ARCHITECT]: Safe mathematical operations with consistent error handling
class SafeMath {
public:
    // Safe division with consistent zero protection
    static inline double safe_divide(double numerator, double denominator) {
        const double safe_denom = (std::abs(denominator) < PrecisionLevels::MIN_SAFE_DIVISOR) 
                                 ? PrecisionLevels::MIN_SAFE_DIVISOR 
                                 : denominator;
        return numerator / safe_denom;
    }
    
    // Safe square root with negative protection
    static inline double safe_sqrt(double value) {
        return std::sqrt(std::max(0.0, value));
    }
    
    // Safe acos with domain clamping
    static inline double safe_acos(double value) {
        const double clamped = std::max(-1.0, std::min(1.0, value));
        return std::acos(clamped);
    }
    
   
    
    // Consistent zero checking
    static inline bool is_effectively_zero(double value, double tolerance = PrecisionLevels::STANDARD_PRECISE) {
        return std::abs(value) < tolerance;
    }
    
    // Consistent equality checking
    static inline bool are_equal(double a, double b, double tolerance = PrecisionLevels::STANDARD_PRECISE) {
        return is_effectively_zero(a - b, tolerance);
    }
};

// [ASPIE ARCHITECT]: Validation and consistency checking
class ValidationLevels {
public:
    // Input validation tolerances
    static constexpr double MIN_VALID_RADIUS = PrecisionLevels::ORIGIN_MARKER;
    static constexpr double MAX_VALID_RADIUS = PrecisionLevels::MAX_SAFE_VALUE;
    static constexpr double MIN_VALID_ANGLE = -MathematicalConstants::PI;
    static constexpr double MAX_VALID_ANGLE = MathematicalConstants::PI;
    
    // State validation
    static constexpr double MIN_FINITE_VALUE = -PrecisionLevels::MAX_SAFE_VALUE;
    static constexpr double MAX_FINITE_VALUE = PrecisionLevels::MAX_SAFE_VALUE;
    
    // Performance validation
    static constexpr size_t MAX_ITERATION_COUNT = 1000000;
    static constexpr double MAX_COMPUTATION_TIME_SECONDS = 1.0;
};



} // namespace precision
} // namespace core
} // namespace hsml