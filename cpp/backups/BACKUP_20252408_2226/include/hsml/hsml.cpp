#pragma once

// ============================================================================
// HSML - Hyper-Spherical Markup Language Framework
// Version 2.0 - Unified Multi-Paradigm Architecture
// 
// Core Tenet: PURE SPHERICAL COORDINATES ONLY
// NO CARTESIAN COORDINATES ARE EVER USED
// 
// This framework provides a complete implementation of spherical coordinate
// rendering with zero Cartesian compromise, featuring:
// - Compile-time spherical mathematics
// - Lock-free concurrent rendering pipeline
// - Type-safe state tensors with physical validation
// - Heterogeneous computing with automatic GPU/CPU dispatch
// - Property-based testing framework
// - Zero-overhead performance monitoring
// ============================================================================

#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <cmath>

// ============================================================================
// CORE MATHEMATICAL FOUNDATION
// ============================================================================

#include "hsml/core/mathematical_foundation.h"
#include "hsml/core/optimized_vector3.h"
#include "hsml/core/unified_spherical_coords.h"

// Compatibility layer for existing code
#include "hsml/core/vector3.h"
#include "hsml/core/matrix4.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"

// ============================================================================
// CORE DATA STRUCTURES AND SCENE GRAPH
// ============================================================================

#include "hsml/core/bubble.h"
#include "hsml/core/presence.h"
#include "hsml/core/hsml_dom.h"
#include "hsml/core/state_tensor_modern.hpp"
#include "hsml/core/object_registry.h"
#include "hsml/core/object_factory.h"

// ============================================================================
// RENDERING SYSTEM
// ============================================================================

#include "hsml/rendering/renderer_interface.h"
#include "hsml/rendering/software_renderer.h"

#ifdef HSML_USE_OPENGL
#include "hsml/rendering/opengl_renderer.h"
#endif

// ============================================================================
// LANGUAGE PARSING AND CODE GENERATION
// ============================================================================

#include "hsml/parsers/hsml_parser.h"
#include "hsml/parsers/hsml_lexer.h"
#include "hsml/parsers/hsml_semantic_analyzer.h"
#include "hsml/codegen/hsml_code_generator.h"

// ============================================================================
// OPTIONAL SIMD OPTIMIZATIONS
// ============================================================================

#ifdef HSML_ENABLE_SIMD
#include "hsml/core/simd_math.h"
#endif

// ============================================================================
// OPTIONAL SDT INTEGRATION
// ============================================================================

#ifdef HSML_SDT_INTEGRATION_ENABLED
#include "hsml/core/hcs21_state_vector.h"
#include "hsml/core/compensation_solver.h"
#include "hsml/core/spherical_physics_engine.h"
#include "hsml/core/sdt_integration_coordinator.h"
#endif

// ============================================================================
// OPTIONAL GUI COMPONENTS
// ============================================================================

#ifdef HSML_USE_QT6
#include "hsml/gui/hsml_visual_studio_gui.h"
#include "hsml/gui/unified_gui_coordinator.h"
#endif

namespace hsml {

// ============================================================================
// FRAMEWORK VERSION AND CAPABILITIES
// ============================================================================

/**
 * @brief Framework information and compile-time feature detection
 * 
 * This struct provides version information and compile-time detection
 * of available features and optimizations.
 * 
 * @example
 * ```cpp
 * std::cout << "HSML Version: " << FrameworkInfo::version_string << "\n";
 * std::cout << "SIMD Available: " << FrameworkInfo::has_simd << "\n";
 * std::cout << "OpenGL Available: " << FrameworkInfo::has_opengl << "\n";
 * ```
 */
struct FrameworkInfo {
    /// Major version number
    static constexpr int major_version = 2;
    
    /// Minor version number
    static constexpr int minor_version = 0;
    
    /// Patch version number
    static constexpr int patch_version = 0;
    
    /// Complete version string
    static constexpr const char* version_string = "2.0.0";
    
    /// Build date (compile-time)
    static constexpr const char* build_date = __DATE__;
    
    /// Build time (compile-time)
    static constexpr const char* build_time = __TIME__;
    
    // ========================================================================
    // COMPILE-TIME FEATURE DETECTION
    // ========================================================================
    
    /// True if SIMD optimizations are enabled
    static constexpr bool has_simd = 
#ifdef HSML_ENABLE_SIMD
        true;
#else
        false;
#endif
    
    /// True if OpenGL rendering backend is available
    static constexpr bool has_opengl = 
#ifdef HSML_USE_OPENGL
        true;
#else
        false;
#endif
    
    /// True if Qt6 GUI components are available
    static constexpr bool has_qt6 = 
#ifdef HSML_USE_QT6
        true;
#else
        false;
#endif
    
    /// True if SDT integration is enabled
    static constexpr bool has_sdt_integration = 
#ifdef HSML_SDT_INTEGRATION_ENABLED
        true;
#else
        false;
#endif
    
    /// True if performance monitoring is enabled
    static constexpr bool has_performance_monitoring = 
#ifdef HSML_ENABLE_PERFORMANCE_COUNTERS
        true;
#else
        false;
#endif
};

// ============================================================================
// FRAMEWORK INITIALIZATION AND CONFIGURATION
// ============================================================================

/**
 * @brief Framework initialization and configuration
 * 
 * This class manages the global initialization and configuration
 * of the HSML framework, including:
 * - Memory pool setup
 * - Thread pool configuration
 * - Performance monitoring initialization
 * - Feature detection and validation
 * 
 * @example
 * ```cpp
 * // Initialize framework with default settings
 * hsml::Framework framework;
 * framework.initialize();
 * 
 * // Configure with custom settings
 * hsml::FrameworkConfig config;
 * config.enable_simd = true;
 * config.enable_performance_monitoring = true;
 * config.thread_pool_size = 8;
 * 
 * framework.initialize(config);
 * ```
 */
class Framework {
public:
    /**
     * @brief Framework configuration options
     */
    struct Config {
        /// Enable SIMD optimizations
        bool enable_simd = false;
        
        /// Enable performance monitoring
        bool enable_performance_monitoring = false;
        
        /// Number of threads in the global thread pool
        size_t thread_pool_size = std::thread::hardware_concurrency();
        
        /// Enable memory pool optimization
        bool enable_memory_pool = true;
        
        /// Enable compile-time validation
        bool enable_compile_time_validation = true;
        
        /// Enable runtime validation
        bool enable_runtime_validation = true;
    };
    
    /**
     * @brief Initialize the framework with default configuration
     * 
     * @return true if initialization succeeded
     * @return false if initialization failed
     */
    static bool initialize();
    
    /**
     * @brief Initialize the framework with custom configuration
     * 
     * @param config Configuration options
     * @return true if initialization succeeded
     * @return false if initialization failed
     */
    static bool initialize(const Config& config);
    
    /**
     * @brief Shutdown the framework and cleanup resources
     */
    static void shutdown();
    
    /**
     * @brief Get current framework configuration
     * 
     * @return Current configuration
     */
    static Config get_config();
    
    /**
     * @brief Check if framework is initialized
     * 
     * @return true if framework is initialized
     * @return false if framework is not initialized
     */
    static bool is_initialized();
    
    /**
     * @brief Get framework information
     * 
     * @return Framework information
     */
    static FrameworkInfo get_info();
};

// ============================================================================
// GLOBAL UTILITIES AND HELPERS
// ============================================================================

/**
 * @brief Global utility functions for common operations
 */
namespace utils {
    
    /**
     * @brief Convert degrees to radians
     * 
     * @param degrees Angle in degrees
     * @return Angle in radians
     */
    template<typename T>
    constexpr T degrees_to_radians(T degrees) noexcept {
        return degrees * T{M_PI} / T{180};
    }
    
    /**
     * @brief Convert radians to degrees
     * 
     * @param radians Angle in radians
     * @return Angle in degrees
     */
    template<typename T>
    constexpr T radians_to_degrees(T radians) noexcept {
        return radians * T{180} / T{M_PI};
    }
    
    /**
     * @brief Clamp value between min and max
     * 
     * @param value Value to clamp
     * @param min Minimum value
     * @param max Maximum value
     * @return Clamped value
     */
    template<typename T>
    constexpr T clamp(T value, T min, T max) noexcept {
        return std::clamp(value, min, max);
    }
    
    /**
     * @brief Linear interpolation between two values
     * 
     * @param a First value
     * @param b Second value
     * @param t Interpolation factor [0, 1]
     * @return Interpolated value
     */
    template<typename T>
    constexpr T lerp(T a, T b, T t) noexcept {
        return a + t * (b - a);
    }
    
    /**
     * @brief Smooth step interpolation
     * 
     * @param edge0 Lower edge
     * @param edge1 Upper edge
     * @param x Input value
     * @return Smoothly interpolated value
     */
    template<typename T>
    constexpr T smoothstep(T edge0, T edge1, T x) noexcept {
        T t = clamp((x - edge0) / (edge1 - edge0), T{0}, T{1});
        return t * t * (T{3} - T{2} * t);
    }
}

// ============================================================================
// ERROR HANDLING AND EXCEPTIONS
// ============================================================================

/**
 * @brief Base exception class for HSML framework
 */
class HSMLException : public std::exception {
public:
    explicit HSMLException(const std::string& message) : message_(message) {}
    
    const char* what() const noexcept override {
        return message_.c_str();
    }
    
private:
    std::string message_;
};

/**
 * @brief Exception thrown for invalid spherical coordinates
 */
class InvalidSphericalCoordsException : public HSMLException {
public:
    explicit InvalidSphericalCoordsException(const std::string& message) 
        : HSMLException("Invalid spherical coordinates: " + message) {}
};

/**
 * @brief Exception thrown for rendering errors
 */
class RenderingException : public HSMLException {
public:
    explicit RenderingException(const std::string& message) 
        : HSMLException("Rendering error: " + message) {}
};

/**
 * @brief Exception thrown for parsing errors
 */
class ParsingException : public HSMLException {
public:
    explicit ParsingException(const std::string& message) 
        : HSMLException("Parsing error: " + message) {}
};

// ============================================================================
// DEBUGGING AND LOGGING
// ============================================================================

/**
 * @brief Logging levels
 */
enum class LogLevel {
    TRACE,    /// Detailed trace information
    DEBUG,    /// Debug information
    INFO,     /// General information
    WARNING,  /// Warning messages
    ERROR,    /// Error messages
    FATAL     /// Fatal errors
};

/**
 * @brief Global logging interface
 */
class Logger {
public:
    /**
     * @brief Set the minimum log level
     * 
     * @param level Minimum log level
     */
    static void set_level(LogLevel level);
    
    /**
     * @brief Get the current log level
     * 
     * @return Current log level
     */
    static LogLevel get_level();
    
    /**
     * @brief Log a message
     * 
     * @param level Log level
     * @param message Message to log
     */
    static void log(LogLevel level, const std::string& message);
    
    /**
     * @brief Log a trace message
     * 
     * @param message Message to log
     */
    static void trace(const std::string& message);
    
    /**
     * @brief Log a debug message
     * 
     * @param message Message to log
     */
    static void debug(const std::string& message);
    
    /**
     * @brief Log an info message
     * 
     * @param message Message to log
     */
    static void info(const std::string& message);
    
    /**
     * @brief Log a warning message
     * 
     * @param message Message to log
     */
    static void warning(const std::string& message);
    
    /**
     * @brief Log an error message
     * 
     * @param message Message to log
     */
    static void error(const std::string& message);
    
    /**
     * @brief Log a fatal message
     * 
     * @param message Message to log
     */
    static void fatal(const std::string& message);
};

// ============================================================================
// COMPILE-TIME VALIDATION MACROS
// ============================================================================

/**
 * @brief Assert that coordinates are valid spherical coordinates
 * 
 * @param coords Spherical coordinates to validate
 */
#define HSML_ASSERT_VALID_SPHERICAL_COORDS(coords) \
    static_assert((coords).is_valid(), "Invalid spherical coordinates")

/**
 * @brief Assert that state tensor is physically valid
 * 
 * @param tensor State tensor to validate
 */
#define HSML_ASSERT_VALID_STATE_TENSOR(tensor) \
    static_assert((tensor).is_valid(), "Invalid state tensor")

/**
 * @brief Assert that framework is initialized
 */
#define HSML_ASSERT_INITIALIZED() \
    assert(hsml::Framework::is_initialized() && "HSML framework not initialized")

// ============================================================================
// PERFORMANCE MONITORING MACROS
// ============================================================================

#ifdef HSML_ENABLE_PERFORMANCE_COUNTERS

/**
 * @brief Measure performance of a code block
 * 
 * @param name Name of the performance measurement
 * @param code Code block to measure
 */
#define HSML_PERFORMANCE_SCOPE(name, code) \
    do { \
        auto scope = hsml::core::performance_scope(name); \
        code; \
    } while(0)

/**
 * @brief Measure performance of a function
 * 
 * @param name Name of the performance measurement
 */
#define HSML_PERFORMANCE_FUNCTION(name) \
    hsml::core::performance_scope scope(name)

#else

#define HSML_PERFORMANCE_SCOPE(name, code) code
#define HSML_PERFORMANCE_FUNCTION(name)

#endif

} // namespace hsml

// ============================================================================
// GLOBAL NAMESPACE ALIASES FOR CONVENIENCE
// ============================================================================

/// Alias for hsml namespace
namespace HSML = hsml;

/// Alias for hsml::core namespace
namespace HSMLCore = hsml::core;

/// Alias for hsml::rendering namespace
namespace HSMLRendering = hsml::rendering;

/// Alias for hsml::compute namespace
namespace HSMLCompute = hsml::compute;

/// Alias for hsml::testing namespace
namespace HSMLTesting = hsml::testing;