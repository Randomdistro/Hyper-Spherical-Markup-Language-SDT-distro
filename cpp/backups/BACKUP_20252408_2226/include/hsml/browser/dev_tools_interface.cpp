/**
 * Developer Tools Interface - C++20 Implementation
 * Revolutionary 3D developer tools for spatial web debugging
 * Multi-paradigm debugging and profiling system
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

// Existing HSML foundation
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"

namespace p0rt3r::devtools {

/**
 * // [The Modern Hipster]: Revolutionary 3D spatial debugging
 * Developer tools interface for spatial web development
 */
class DevToolsInterface {
public:
    struct DevToolsConfig {
        bool enable_3d_inspector = true;
        bool enable_spatial_profiler = true;
        bool enable_performance_monitor = true;
        bool enable_network_inspector = true;
        bool enable_console_api = true;
        
        // 3D visualization settings
        bool show_spatial_grid = true;
        bool show_coordinate_axes = true;
        bool highlight_selected_elements = true;
        double visualization_scale = 1.0;
        
        DevToolsConfig() = default;
    };

private:
    DevToolsConfig config_;
    std::atomic<bool> is_enabled_{false};

public:
    explicit DevToolsInterface(const DevToolsConfig& config = {}) : config_(config) {}
    
    // DevTools lifecycle
    void enable() { is_enabled_.store(true); }
    void disable() { is_enabled_.store(false); }
    bool is_enabled() const noexcept { return is_enabled_.load(); }
    
    // 3D spatial debugging
    void inspect_spatial_element(std::string_view element_id);
    void highlight_coordinate_region(const hsml::core::SphericalCoords& center, double radius);
    void show_navigation_path(const std::vector<hsml::core::SphericalCoords>& path);
    
    // Performance profiling
    void start_spatial_profiling();
    void stop_spatial_profiling();
    [[nodiscard]] std::string get_profiling_report() const;
    
    // Configuration
    const DevToolsConfig& config() const noexcept { return config_; }
    void update_config(const DevToolsConfig& new_config) { config_ = new_config; }
};

} // namespace p0rt3r::devtools