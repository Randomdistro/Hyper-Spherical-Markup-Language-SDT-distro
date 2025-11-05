/**
 * Spatial JavaScript Engine - C++20 Implementation
 * Revolutionary JavaScript integration with spatial web APIs
 * Multi-paradigm V8 integration with HSML spatial bindings
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <expected>
#include <concepts>
#include <functional>
#include <atomic>
#include <chrono>
#include <future>
#include <variant>

// Existing HSML foundation
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"
#include "hsml/core/vector3.h"

// Browser components
#include "spatial_document_parser.h"
#include "navigation_manager.h"

namespace p0rt3r::javascript {

// Forward declarations
class SpatialJavaScriptContext;
class SpatialAPIBindings;
class JSObjectWrapper;

/**
 * // [The Modern Hipster]: Cutting-edge V8 integration
 * JavaScript value representation with spatial extensions
 */
class SpatialJSValue {
public:
    enum class ValueType {
        UNDEFINED,
        NULL_VALUE,
        BOOLEAN,
        NUMBER,
        STRING,
        OBJECT,
        ARRAY,
        FUNCTION,
        SPHERICAL_COORDS,    // Native spatial type
        SOLID_ANGLE,         // Native spatial type
        VECTOR3,             // Native spatial type
        SPATIAL_ELEMENT      // Reference to spatial element
    };

private:
    ValueType type_;
    std::variant<
        std::monostate,                           // UNDEFINED/NULL
        bool,                                     // BOOLEAN
        double,                                   // NUMBER
        std::string,                              // STRING
        std::shared_ptr<JSObjectWrapper>,         // OBJECT/ARRAY/FUNCTION
        hsml::core::SphericalCoords,              // SPHERICAL_COORDS
        hsml::core::SolidAngle,                   // SOLID_ANGLE
        hsml::core::Vector3,                      // VECTOR3
        std::weak_ptr<parser::SpatialElement>     // SPATIAL_ELEMENT
    > value_;

public:
    SpatialJSValue() : type_(ValueType::UNDEFINED) {}
    
    // Constructors for different types
    explicit SpatialJSValue(bool value) : type_(ValueType::BOOLEAN), value_(value) {}
    explicit SpatialJSValue(double value) : type_(ValueType::NUMBER), value_(value) {}
    explicit SpatialJSValue(std::string_view value) : type_(ValueType::STRING), value_(std::string(value)) {}
    explicit SpatialJSValue(const hsml::core::SphericalCoords& coords) : type_(ValueType::SPHERICAL_COORDS), value_(coords) {}
    explicit SpatialJSValue(const hsml::core::SolidAngle& angle) : type_(ValueType::SOLID_ANGLE), value_(angle) {}
    explicit SpatialJSValue(const hsml::core::Vector3& vec) : type_(ValueType::VECTOR3), value_(vec) {}
    
    // Type checking
    ValueType type() const noexcept { return type_; }
    bool is_spatial_type() const noexcept {
        return type_ == ValueType::SPHERICAL_COORDS || 
               type_ == ValueType::SOLID_ANGLE || 
               type_ == ValueType::VECTOR3 ||
               type_ == ValueType::SPATIAL_ELEMENT;
    }
    
    // Value access
    template<typename T>
    [[nodiscard]] auto get() const -> std::expected<T, std::string> {
        if (auto* val = std::get_if<T>(&value_)) {
            return *val;
        }
        return std::unexpected("Type mismatch in SpatialJSValue::get()");
    }
    
    // Conversion to/from V8 values
    // These would be implemented with actual V8 integration
    void* to_v8_value() const;
    static SpatialJSValue from_v8_value(void* v8_value);
    
    // String representation
    [[nodiscard]] std::string to_string() const;
    [[nodiscard]] std::string to_json() const;
};

/**
 * // [The Security Paranoid]: Sandboxed JavaScript execution context
 * Secure isolated execution environment for spatial scripts
 */
class SpatialJavaScriptContext {
public:
    struct ContextConfig {
        std::string context_id;
        std::string document_uri;
        
        // Security settings
        bool enable_network_access = true;
        bool enable_file_system_access = false;
        bool enable_native_modules = false;
        std::vector<std::string> allowed_domains;
        std::vector<std::string> blocked_apis;
        
        // Performance limits
        size_t max_heap_size_mb = 64;
        std::chrono::milliseconds execution_timeout{5000};
        size_t max_recursion_depth = 1000;
        size_t max_array_length = 1000000;
        
        // Spatial API access
        bool enable_spatial_navigation = true;
        bool enable_physics_simulation = true;
        bool enable_dom_manipulation = true;
        bool enable_rendering_control = false;  // Restricted by default
        
        ContextConfig(std::string_view id, std::string_view uri) 
            : context_id(id), document_uri(uri) {}
    };

private:
    ContextConfig config_;
    void* v8_isolate_ = nullptr;        // V8::Isolate*
    void* v8_context_ = nullptr;        // V8::Local<V8::Context>
    std::unique_ptr<SpatialAPIBindings> api_bindings_;
    
    // Execution state
    std::atomic<bool> is_initialized_{false};
    std::atomic<bool> is_executing_{false};
    std::atomic<uint64_t> scripts_executed_{0};
    std::atomic<size_t> memory_usage_bytes_{0};
    
    // Error handling
    std::vector<std::string> compilation_errors_;
    std::vector<std::string> runtime_errors_;
    
    // Security sandbox
    std::unordered_set<std::string> accessible_apis_;
    std::function<bool(std::string_view)> security_validator_;

public:
    explicit SpatialJavaScriptContext(const ContextConfig& config);
    ~SpatialJavaScriptContext();
    
    // Context lifecycle
    [[nodiscard]] auto initialize() -> std::expected<void, std::string>;
    void shutdown();
    bool is_initialized() const noexcept { return is_initialized_.load(); }
    
    // Script execution
    [[nodiscard]] auto execute_script(std::string_view script_source, 
                                     std::string_view script_name = "anonymous")
        -> std::expected<SpatialJSValue, std::string>;
    
    [[nodiscard]] auto execute_script_async(std::string_view script_source, 
                                           std::string_view script_name = "anonymous")
        -> std::future<std::expected<SpatialJSValue, std::string>>;
    
    // Function calls
    [[nodiscard]] auto call_function(std::string_view function_name,
                                    const std::vector<SpatialJSValue>& arguments)
        -> std::expected<SpatialJSValue, std::string>;
    
    // Global object management
    void set_global_property(std::string_view name, const SpatialJSValue& value);
    [[nodiscard]] auto get_global_property(std::string_view name) -> std::expected<SpatialJSValue, std::string>;
    
    // API bindings
    SpatialAPIBindings& api_bindings() { return *api_bindings_; }
    void register_native_function(std::string_view name, 
                                 std::function<SpatialJSValue(const std::vector<SpatialJSValue>&)> callback);
    
    // Context information
    const ContextConfig& config() const noexcept { return config_; }
    const std::string& context_id() const noexcept { return config_.context_id; }
    
    // Performance monitoring
    uint64_t scripts_executed() const noexcept { return scripts_executed_.load(); }
    size_t memory_usage_bytes() const noexcept { return memory_usage_bytes_.load(); }
    double memory_usage_mb() const noexcept { return static_cast<double>(memory_usage_bytes()) / (1024.0 * 1024.0); }
    
    // Error handling
    const std::vector<std::string>& compilation_errors() const noexcept { return compilation_errors_; }
    const std::vector<std::string>& runtime_errors() const noexcept { return runtime_errors_; }
    void clear_errors();
    
    // Security
    void enable_api_access(std::string_view api_name);
    void disable_api_access(std::string_view api_name);
    bool has_api_access(std::string_view api_name) const;
};

/**
 * // [The Performance Demon]: High-performance spatial API bindings
 * Native C++ functions exposed to JavaScript with SIMD optimization
 */
class SpatialAPIBindings {
private:
    SpatialJavaScriptContext* context_;
    navigation::NavigationManager* navigation_manager_ = nullptr;
    
    // Performance optimization
    mutable std::unordered_map<std::string, void*> cached_function_templates_;

public:
    explicit SpatialAPIBindings(SpatialJavaScriptContext* context) : context_(context) {}
    
    // API exposure
    void expose_spherical_coordinates_api();
    void expose_solid_angle_api();
    void expose_vector3_api();
    void expose_spatial_navigation_api();
    void expose_physics_simulation_api();
    void expose_dom_manipulation_api();
    void expose_rendering_control_api();
    
    // Navigation API implementation
    void set_navigation_manager(navigation::NavigationManager* nav_manager) { navigation_manager_ = nav_manager; }
    
    // // [The Functional Purist]: Pure spatial API functions
    [[nodiscard]] static auto js_spherical_coords_create(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    [[nodiscard]] static auto js_spherical_coords_distance(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    [[nodiscard]] static auto js_spherical_coords_interpolate(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    
    [[nodiscard]] static auto js_solid_angle_create(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    [[nodiscard]] static auto js_solid_angle_omega(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    
    [[nodiscard]] static auto js_vector3_create(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    [[nodiscard]] static auto js_vector3_dot(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    [[nodiscard]] static auto js_vector3_cross(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    
    // Navigation API functions
    [[nodiscard]] auto js_navigation_to_position(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    [[nodiscard]] auto js_navigation_teleport(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    [[nodiscard]] auto js_navigation_dive_depth(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    
    // Physics API functions
    [[nodiscard]] auto js_physics_set_gravity(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    [[nodiscard]] auto js_physics_apply_force(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
    [[nodiscard]] auto js_physics_create_constraint(const std::vector<SpatialJSValue>& args) -> SpatialJSValue;
};

/**
 * // [The OOP Architect]: Comprehensive JavaScript engine
 * Main engine coordinating multiple contexts and API bindings
 */
template<typename EngineBackend = void> // Defaulted for V8
class SpatialJavaScriptEngine {
public:
    struct EngineConfig {
        // Engine settings
        size_t max_heap_size_mb = 256;
        uint32_t max_contexts = 16;
        bool enable_debugging = false;
        bool enable_profiling = false;
        
        // Security settings
        bool enable_unsafe_eval = false;
        bool enable_wasm = true;
        std::vector<std::string> trusted_origins;
        
        // Performance settings
        bool enable_jit_compilation = true;
        bool enable_optimization = true;
        uint32_t compilation_cache_size_mb = 64;
        
        // Spatial integration
        bool auto_expose_spatial_apis = true;
        navigation::NavigationManager* navigation_manager = nullptr;
        
        EngineConfig() = default;
    };

private:
    EngineConfig config_;
    void* v8_platform_ = nullptr;     // V8::Platform*
    
    // Context management
    std::unordered_map<std::string, std::unique_ptr<SpatialJavaScriptContext>> contexts_;
    std::atomic<uint32_t> next_context_id_{1};
    
    // Performance metrics
    mutable std::atomic<uint64_t> total_scripts_executed_{0};
    mutable std::atomic<uint64_t> total_compilation_time_ms_{0};
    mutable std::atomic<uint64_t> total_execution_time_ms_{0};
    mutable std::atomic<size_t> peak_memory_usage_mb_{0};
    
    // Engine state
    std::atomic<bool> is_initialized_{false};
    std::vector<std::jthread> worker_threads_;
    std::atomic<bool> should_terminate_{false};

public:
    explicit SpatialJavaScriptEngine(const EngineConfig& config = {}) : config_(config) {}
    ~SpatialJavaScriptEngine();
    
    // Engine lifecycle
    [[nodiscard]] auto initialize() -> std::expected<void, std::string>;
    void shutdown();
    bool is_initialized() const noexcept { return is_initialized_.load(); }
    
    // Context management
    [[nodiscard]] auto create_context(std::string_view document_uri,
                                     const SpatialJavaScriptContext::ContextConfig& context_config = {})
        -> std::expected<std::string, std::string>;
    
    void destroy_context(std::string_view context_id);
    [[nodiscard]] SpatialJavaScriptContext* get_context(std::string_view context_id) const;
    [[nodiscard]] std::vector<std::string> get_all_context_ids() const;
    
    // Script execution
    [[nodiscard]] auto execute_script(std::string_view script_source, 
                                     std::string_view context_id,
                                     std::string_view script_name = "anonymous")
        -> std::expected<SpatialJSValue, std::string>;
    
    // API exposure
    void expose_spherical_coordinates_api();
    void expose_solid_angle_api();
    void expose_spatial_navigation_api();
    void expose_physics_simulation_api();
    
    // Configuration
    void update_config(const EngineConfig& new_config);
    const EngineConfig& config() const noexcept { return config_; }
    
    // Performance metrics
    struct JSEngineMetrics {
        uint64_t total_scripts_executed;
        uint64_t total_compilation_time_ms;
        uint64_t total_execution_time_ms;
        size_t peak_memory_usage_mb;
        uint32_t active_contexts;
        size_t compilation_cache_size_mb;
        
        // Spatial-specific metrics
        uint64_t spatial_api_calls;
        uint64_t navigation_operations;
        uint64_t physics_calculations;
        double average_spatial_call_time_ms;
    };
    
    [[nodiscard]] JSEngineMetrics get_metrics() const;
    
    // Debugging and profiling
    void enable_debugging(bool enable);
    void enable_profiling(bool enable);
    [[nodiscard]] std::string get_profiling_data() const;
    void dump_heap_snapshot(std::string_view filename) const;
    
    // Error handling
    using ErrorCallback = std::function<void(std::string_view context_id, std::string_view error)>;
    void set_error_callback(ErrorCallback callback);
    
    // Garbage collection
    void trigger_garbage_collection();
    void set_memory_pressure_level(int level); // 0 = none, 1 = moderate, 2 = critical
};

// Utility functions and type definitions
namespace js_utils {
    // // [The Hacktivist]: Quick JavaScript integration helpers
    template<typename T>
    [[nodiscard]] inline auto to_spatial_js_value(const T& value) -> SpatialJSValue {
        if constexpr (std::is_same_v<T, hsml::core::SphericalCoords>) {
            return SpatialJSValue(value);
        } else if constexpr (std::is_same_v<T, hsml::core::SolidAngle>) {
            return SpatialJSValue(value);
        } else if constexpr (std::is_same_v<T, hsml::core::Vector3>) {
            return SpatialJSValue(value);
        } else if constexpr (std::is_arithmetic_v<T>) {
            return SpatialJSValue(static_cast<double>(value));
        } else if constexpr (std::is_convertible_v<T, std::string>) {
            return SpatialJSValue(std::string(value));
        } else {
            static_assert(sizeof(T) == 0, "Unsupported type for SpatialJSValue conversion");
        }
    }
    
    [[nodiscard]] inline auto validate_spatial_coordinates(const SpatialJSValue& value) -> bool {
        if (value.type() != SpatialJSValue::ValueType::SPHERICAL_COORDS) return false;
        
        auto coords_result = value.get<hsml::core::SphericalCoords>();
        if (!coords_result) return false;
        
        const auto& coords = coords_result.value();
        return coords.r() >= 0.0 && coords.theta() >= 0.0 && coords.theta() <= M_PI;
    }
    
    [[nodiscard]] inline auto create_spatial_error(std::string_view message) -> SpatialJSValue {
        // This would create a proper JavaScript Error object
        return SpatialJSValue(std::string("SpatialError: ") + std::string(message));
    }
}

// Factory functions
[[nodiscard]] auto create_spatial_javascript_engine(
    const SpatialJavaScriptEngine<>::EngineConfig& config = {})
    -> std::unique_ptr<SpatialJavaScriptEngine<>>;

// Global JavaScript console for debugging
namespace spatial_console {
    void log(std::string_view message);
    void warn(std::string_view message);
    void error(std::string_view message);
    void debug(std::string_view message);
    
    // Spatial-specific console functions
    void log_coordinates(const hsml::core::SphericalCoords& coords);
    void log_solid_angle(const hsml::core::SolidAngle& angle);
    void log_navigation_state(const navigation::NavigationManager& nav_manager);
}

} // namespace p0rt3r::javascript