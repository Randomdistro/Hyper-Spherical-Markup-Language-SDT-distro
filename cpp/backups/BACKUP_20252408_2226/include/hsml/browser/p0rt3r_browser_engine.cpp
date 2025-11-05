/**
 * P0rt3r Browser Engine - C++20 Implementation
 * World's first native HSML browser engine with revolutionary spatial UI
 * Integrates existing HSML foundation components for maximum performance
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <future>
#include <expected>
#include <concepts>
#include <atomic>
#include <chrono>
#include <thread>

// Existing HSML foundation components
#include "hsml/core/solid_angle_dom_processor.h"
#include "hsml/core/spherical_coordinate_processor.h"
#include "hsml/rendering/spherical_renderer.h"
#include "hsml/parsers/hsml_parser.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"

// Forward declarations
namespace p0rt3r {
    class SpatialDocument;
    class NavigationManager;
    class ResourceLoader;
    class JavaScriptEngine;
    class DevToolsInterface;
    enum class GraphicsAPI;
}

namespace p0rt3r::engine {

// Resource types for document loading
enum class ResourceType {
    LOCAL_FILE,
    REMOTE_HTTP,
    REMOTE_HSML_PROTOCOL,
    EMBEDDED_RESOURCE,
    CACHED_CONTENT
};

// Document handle for managing loaded documents
using DocumentHandle = uint64_t;
constexpr DocumentHandle INVALID_DOCUMENT = 0;

// Navigation promise for asynchronous spatial navigation
using NavigationPromise = std::future<std::expected<void, std::string>>;

// Load error types
enum class LoadError {
    NETWORK_ERROR,
    PARSE_ERROR,
    RESOURCE_NOT_FOUND,
    PERMISSION_DENIED,
    UNSUPPORTED_FORMAT,
    OUT_OF_MEMORY
};

// Revolutionary browser engine leveraging existing HSML components
template<GraphicsAPI API = GraphicsAPI::OPENGL_4_5>
class P0rt3rBrowserEngine {
private:
    // Core HSML components (existing foundation)
    std::unique_ptr<hsml::core::SolidAngleDOMProcessor> dom_processor_;
    std::unique_ptr<hsml::core::SphericalCoordinateProcessor> coord_processor_;
    std::unique_ptr<hsml::rendering::SphericalRenderer<API>> renderer_;
    
    // New P0rt3r-specific components
    std::unique_ptr<SpatialDocumentParser> document_parser_;
    std::unique_ptr<SpatialLayoutEngine> layout_engine_;
    std::unique_ptr<NavigationManager> navigation_manager_;
    std::unique_ptr<ResourceLoader> resource_loader_;
    std::unique_ptr<JavaScriptEngine> js_engine_;
    
    // Multi-threaded rendering pipeline
    std::unique_ptr<RenderingPipeline<API>> rendering_pipeline_;
    std::unique_ptr<NetworkingLayer> networking_layer_;
    
    // Developer tools
    std::unique_ptr<DevToolsInterface> dev_tools_;
    
    // Document management
    std::unordered_map<DocumentHandle, std::unique_ptr<SpatialDocument>> loaded_documents_;
    std::atomic<DocumentHandle> next_document_handle_{1};
    DocumentHandle active_document_{INVALID_DOCUMENT};
    
    // Performance metrics
    mutable std::atomic<uint64_t> documents_loaded_{0};
    mutable std::atomic<uint64_t> spatial_objects_rendered_{0};
    mutable std::atomic<double> average_frame_time_ms_{0.0};
    mutable std::atomic<double> navigation_latency_ms_{0.0};
    mutable std::atomic<size_t> memory_usage_mb_{0};
    mutable std::atomic<uint32_t> active_javascript_contexts_{0};
    
    // Thread management
    std::vector<std::jthread> worker_threads_;
    std::atomic<bool> should_terminate_{false};
    
    // Configuration
    struct BrowserConfig {
        double default_viewer_distance = 650.0;  // mm
        uint32_t max_concurrent_documents = 16;
        uint32_t cache_size_mb = 512;
        bool enable_javascript = true;
        bool enable_webgl_fallback = true;
        GraphicsAPI preferred_api = API;
        uint32_t worker_thread_count = std::thread::hardware_concurrency();
    } config_;
    
    bool initialized_ = false;
    
    // Initialization helpers
    [[nodiscard]] auto initialize_core_components() -> std::expected<void, std::string>;
    [[nodiscard]] auto initialize_browser_components() -> std::expected<void, std::string>;
    [[nodiscard]] auto initialize_rendering_pipeline() -> std::expected<void, std::string>;
    [[nodiscard]] auto initialize_javascript_engine() -> std::expected<void, std::string>;
    
    // Document management helpers
    DocumentHandle generate_document_handle() noexcept;
    void cleanup_document(DocumentHandle handle);
    
    // Worker thread functions
    void rendering_worker_thread();
    void network_worker_thread();
    void javascript_worker_thread();
    
public:
    // Constructor
    explicit P0rt3rBrowserEngine() = default;
    
    // Disable copy/move operations for singleton-like behavior
    P0rt3rBrowserEngine(const P0rt3rBrowserEngine&) = delete;
    P0rt3rBrowserEngine& operator=(const P0rt3rBrowserEngine&) = delete;
    P0rt3rBrowserEngine(P0rt3rBrowserEngine&&) = delete;
    P0rt3rBrowserEngine& operator=(P0rt3rBrowserEngine&&) = delete;
    
    // Destructor
    ~P0rt3rBrowserEngine();
    
    // Initialize browser engine with existing HSML components
    [[nodiscard]] auto initialize(const BrowserConfig& config = {}) 
        -> std::expected<void, std::string>;
    
    // Shutdown browser engine
    void shutdown();
    
    // Load and render HSML document
    template<ResourceType Type = ResourceType::REMOTE_HTTP>
    auto load_document(std::string_view uri) 
        -> std::future<std::expected<DocumentHandle, LoadError>>;
    
    // Document management
    [[nodiscard]] auto get_active_document() const noexcept -> DocumentHandle;
    auto set_active_document(DocumentHandle handle) -> std::expected<void, std::string>;
    auto close_document(DocumentHandle handle) -> void;
    [[nodiscard]] auto get_loaded_documents() const -> std::vector<DocumentHandle>;
    
    // Spatial navigation methods
    auto navigate_to_sphere(const hsml::core::SphericalCoords& target) 
        -> NavigationPromise;
    auto dive_into_depth(double depth_level) -> NavigationPromise;
    auto teleport_to_coordinate(const hsml::core::SphericalCoords& coord) 
        -> NavigationPromise;
    auto navigate_back() -> NavigationPromise;
    auto navigate_forward() -> NavigationPromise;
    
    // Viewport management
    auto set_viewer_position(const hsml::core::SphericalCoords& position) -> void;
    [[nodiscard]] auto get_viewer_position() const -> hsml::core::SphericalCoords;
    auto set_field_of_view(const hsml::core::SolidAngle& fov) -> void;
    [[nodiscard]] auto get_field_of_view() const -> hsml::core::SolidAngle;
    
    // Rendering control
    auto render_frame() -> std::expected<void, std::string>;
    auto start_continuous_rendering() -> void;
    auto stop_continuous_rendering() -> void;
    [[nodiscard]] auto is_rendering() const noexcept -> bool;
    
    // JavaScript integration
    [[nodiscard]] auto get_javascript_engine() -> JavaScriptEngine*;
    auto execute_javascript(std::string_view script, DocumentHandle document = INVALID_DOCUMENT)
        -> std::expected<std::string, std::string>;
    auto enable_javascript_debugging(bool enable) -> void;
    
    // Developer tools integration
    auto get_developer_tools() -> DevToolsInterface&;
    auto enable_developer_mode(bool enable) -> void;
    [[nodiscard]] auto is_developer_mode_enabled() const noexcept -> bool;
    
    // Performance monitoring
    struct BrowserMetrics {
        uint64_t documents_loaded;
        uint64_t spatial_objects_rendered;
        double average_frame_time_ms;
        double navigation_latency_ms;
        size_t memory_usage_mb;
        uint32_t active_javascript_contexts;
        
        // Rendering performance
        double gpu_utilization_percentage;
        size_t gpu_memory_used_mb;
        uint32_t triangles_rendered_per_frame;
        double solid_angle_coverage_percentage;
        
        // Network performance
        size_t bytes_downloaded_mb;
        size_t bytes_cached_mb;
        double average_download_speed_mbps;
        uint32_t active_network_connections;
    };
    
    [[nodiscard]] BrowserMetrics get_metrics() const;
    
    // Configuration management
    auto update_config(const BrowserConfig& new_config) -> std::expected<void, std::string>;
    [[nodiscard]] auto get_config() const noexcept -> const BrowserConfig&;
    
    // Cache management
    auto clear_cache() -> void;
    auto optimize_cache() -> void;
    [[nodiscard]] auto get_cache_usage() const -> size_t;
    
    // Security and privacy
    auto enable_private_browsing(bool enable) -> void;
    [[nodiscard]] auto is_private_browsing_enabled() const noexcept -> bool;
    auto clear_browsing_data() -> void;
    
    // Accessibility support
    auto enable_accessibility_features(bool enable) -> void;
    auto set_accessibility_zoom_level(double zoom_factor) -> void;
    auto enable_high_contrast_mode(bool enable) -> void;
    
    // Event system
    template<typename EventHandler>
    auto register_event_handler(std::string_view event_type, EventHandler&& handler) -> void
        requires std::invocable<EventHandler, const std::string&>;
    
    // Diagnostic and debugging
    [[nodiscard]] auto get_diagnostic_info() const -> std::string;
    auto dump_performance_profile(std::string_view filename) const -> void;
    auto run_self_diagnostics() -> std::vector<std::string>;
    
    // Advanced features
    auto enable_experimental_features(bool enable) -> void;
    [[nodiscard]] auto get_supported_graphics_apis() const -> std::vector<GraphicsAPI>;
    auto switch_graphics_api(GraphicsAPI api) -> std::expected<void, std::string>;
};

// Factory function for creating P0rt3r browser engine
template<GraphicsAPI API = GraphicsAPI::OPENGL_4_5>
[[nodiscard]] auto create_p0rt3r_engine() -> std::unique_ptr<P0rt3rBrowserEngine<API>>;

// Helper functions for browser initialization
namespace utils {
    [[nodiscard]] auto detect_optimal_graphics_api() -> GraphicsAPI;
    [[nodiscard]] auto get_system_capabilities() -> SystemCapabilities;
    [[nodiscard]] auto validate_browser_environment() -> std::vector<std::string>;
}

} // namespace p0rt3r::engine

// Implementation details in separate files
#include "spatial_document_parser.h"
#include "spatial_layout_engine.h"
#include "navigation_manager.h"
#include "resource_loader.h"
#include "javascript_engine.h"
#include "dev_tools_interface.h"