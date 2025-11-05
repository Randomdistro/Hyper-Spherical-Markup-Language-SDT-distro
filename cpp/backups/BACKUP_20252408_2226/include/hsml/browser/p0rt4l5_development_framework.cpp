#pragma once

#include "../core/spherical_coords.h"
#include "../core/solid_angle.h"
#include "../core/simd_math.h"
#include "p0rt3r_engine.h"
#include <memory>
#include <vector>
#include <unordered_map>
#include <string>
#include <functional>
#include <chrono>
#include <atomic>
#include <thread>
#include <future>
#include <coroutine>

namespace hsml {
namespace browser {

// Forward declarations
class P0RT4L5TestingPortal;
class SteradianRecalculationEngine;
class HotSpotManager;
class HSMLDevelopmentAccelerator;

// Portal scaling information for steradian calculations
struct PortalScalingInfo {
    double original_radius;
    double current_radius;
    double scale_factor;
    core::SphericalCoords position;
    core::SolidAngle solid_angle;
    std::chrono::steady_clock::time_point last_update;
    bool is_minimized;
    double hot_spot_intensity;
    
    // Calculate current steradian coverage
    double calculate_steradian_coverage() const {
        return solid_angle.to_steradians() * scale_factor * scale_factor;
    }
    
    // Update scaling based on minimization state
    void update_scaling(bool minimized, double target_scale) {
        is_minimized = minimized;
        scale_factor = target_scale;
        last_update = std::chrono::steady_clock::now();
        
        // Hot spot intensity increases as portal gets smaller
        hot_spot_intensity = minimized ? (1.0 / scale_factor) : 1.0;
    }
};

// P0RT4L5 Development Testing Framework
class P0RT4L5DevelopmentFramework {
public:
    enum class PortalType {
        HSML_TESTING_PORTAL,
        PERFORMANCE_BENCHMARK_PORTAL,
        PHYSICS_SIMULATION_PORTAL,
        RENDERING_TEST_PORTAL,
        SPATIAL_NAVIGATION_PORTAL,
        INTEGRATION_TEST_PORTAL
    };
    
    enum class AccelerationMode {
        FULL_SCALE_DEVELOPMENT,
        RAPID_PROTOTYPING,
        STRESS_TESTING,
        PERFORMANCE_PROFILING,
        REAL_TIME_DEBUGGING
    };

    P0RT4L5DevelopmentFramework();
    ~P0RT4L5DevelopmentFramework();
    
    // Portal management
    AsyncTask<int> create_development_portal(PortalType type, const std::string& config_path);
    AsyncTask<bool> destroy_portal(int portal_id);
    AsyncTask<bool> minimize_portal(int portal_id, double scale_factor = 0.1);
    AsyncTask<bool> maximize_portal(int portal_id);
    AsyncTask<bool> teleport_portal(int portal_id, const core::SphericalCoords& new_position);
    
    // Steradian recalculation system
    void enable_dynamic_steradian_recalculation(bool enabled);
    void set_recalculation_frequency(double frequency_hz);
    PortalScalingInfo get_portal_scaling_info(int portal_id) const;
    std::vector<PortalScalingInfo> get_all_portal_scaling_info() const;
    
    // Hot spot management
    void enable_hot_spot_system(bool enabled);
    void set_hot_spot_sensitivity(double sensitivity);
    std::vector<core::SphericalCoords> get_active_hot_spots() const;
    AsyncTask<bool> trigger_hot_spot_interaction(const core::SphericalCoords& position);
    
    // HSML development acceleration
    void set_acceleration_mode(AccelerationMode mode);
    AccelerationMode get_acceleration_mode() const;
    AsyncTask<bool> run_accelerated_test_suite(const std::string& test_config);
    AsyncTask<bool> benchmark_hsml_performance();
    
    // Real-time testing capabilities
    AsyncTask<bool> inject_test_data(int portal_id, const std::string& test_data);
    AsyncTask<std::string> extract_test_results(int portal_id);
    void start_continuous_testing(int portal_id, std::function<void(const std::string&)> result_callback);
    void stop_continuous_testing(int portal_id);
    
    // Performance monitoring
    struct PerformanceMetrics {
        std::atomic<double> average_steradian_calc_time{0.0};
        std::atomic<size_t> portal_switches_per_second{0};
        std::atomic<double> hot_spot_response_time{0.0};
        std::atomic<size_t> active_portals{0};
        std::atomic<double> memory_usage_mb{0.0};
        std::atomic<double> cpu_usage_percent{0.0};
    };
    
    const PerformanceMetrics& get_performance_metrics() const;
    void reset_performance_metrics();
    
    // Integration with P0RT3R browser
    void attach_to_browser(std::shared_ptr<P0rt3rBrowserEngine> browser);
    void detach_from_browser();
    bool is_attached_to_browser() const;
    
    // Event callbacks
    using PortalCreatedCallback = std::function<void(int portal_id, PortalType type)>;
    using PortalScaledCallback = std::function<void(int portal_id, const PortalScalingInfo& info)>;
    using HotSpotTriggeredCallback = std::function<void(const core::SphericalCoords& position, double intensity)>;
    using TestCompletedCallback = std::function<void(int portal_id, bool success, const std::string& results)>;
    
    void set_portal_created_callback(PortalCreatedCallback callback);
    void set_portal_scaled_callback(PortalScaledCallback callback);
    void set_hot_spot_triggered_callback(HotSpotTriggeredCallback callback);
    void set_test_completed_callback(TestCompletedCallback callback);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
    
    // Core components
    std::unique_ptr<SteradianRecalculationEngine> steradian_engine_;
    std::unique_ptr<HotSpotManager> hot_spot_manager_;
    std::unique_ptr<HSMLDevelopmentAccelerator> accelerator_;
    
    // State management
    std::unordered_map<int, std::unique_ptr<P0RT4L5TestingPortal>> active_portals_;
    std::atomic<int> next_portal_id_{1};
    AccelerationMode current_acceleration_mode_;
    PerformanceMetrics metrics_;
    
    // Browser integration
    std::weak_ptr<P0rt3rBrowserEngine> attached_browser_;
    
    // Internal methods
    void initialize_components();
    void cleanup_components();
    void update_performance_metrics();
    AsyncTask<void> run_steradian_recalculation_loop();
    AsyncTask<void> run_hot_spot_monitoring_loop();
};

// Individual testing portal
class P0RT4L5TestingPortal {
public:
    P0RT4L5TestingPortal(P0RT4L5DevelopmentFramework::PortalType type, int id);
    ~P0RT4L5TestingPortal();
    
    // Portal properties
    int get_id() const { return portal_id_; }
    P0RT4L5DevelopmentFramework::PortalType get_type() const { return portal_type_; }
    
    // Spatial properties
    void set_position(const core::SphericalCoords& coords);
    core::SphericalCoords get_position() const;
    void set_scaling_info(const PortalScalingInfo& info);
    PortalScalingInfo get_scaling_info() const;
    
    // Testing capabilities
    AsyncTask<bool> load_test_configuration(const std::string& config_path);
    AsyncTask<bool> execute_test_sequence();
    AsyncTask<std::string> get_test_results();
    void set_continuous_testing(bool enabled);
    
    // Rendering integration
    void render(const core::SphericalCoords& viewer_position);
    bool is_visible_from(const core::SphericalCoords& viewer_position) const;
    double calculate_angular_size_from(const core::SphericalCoords& viewer_position) const;

private:
    int portal_id_;
    P0RT4L5DevelopmentFramework::PortalType portal_type_;
    PortalScalingInfo scaling_info_;
    std::atomic<bool> continuous_testing_enabled_{false};
    std::thread testing_thread_;
    
    // Test execution state
    struct TestState {
        bool running = false;
        std::string current_test;
        std::chrono::steady_clock::time_point start_time;
        std::vector<std::string> results;
    } test_state_;
    
    void run_continuous_testing();
};

// Steradian recalculation engine for dynamic portal scaling
class SteradianRecalculationEngine {
public:
    SteradianRecalculationEngine();
    ~SteradianRecalculationEngine();
    
    // Configuration
    void set_recalculation_frequency(double frequency_hz);
    void enable_simd_acceleration(bool enabled);
    void enable_parallel_processing(bool enabled);
    
    // Recalculation operations
    AsyncTask<PortalScalingInfo> recalculate_portal_steradians(
        const PortalScalingInfo& current_info,
        const core::SphericalCoords& viewer_position
    );
    
    AsyncTask<std::vector<PortalScalingInfo>> recalculate_all_portals(
        const std::vector<PortalScalingInfo>& portal_infos,
        const core::SphericalCoords& viewer_position
    );
    
    // Performance optimization
    void optimize_calculation_pipeline();
    void warm_up_simd_cache();
    
    // Metrics
    double get_average_calculation_time() const;
    size_t get_calculations_per_second() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
    
    // SIMD-optimized calculation kernels
    PortalScalingInfo calculate_steradians_simd(
        const PortalScalingInfo& info,
        const core::SphericalCoords& viewer_pos
    );
    
    void batch_calculate_steradians_simd(
        std::vector<PortalScalingInfo>& infos,
        const core::SphericalCoords& viewer_pos
    );
};

// Hot spot management system
class HotSpotManager {
public:
    HotSpotManager();
    ~HotSpotManager();
    
    // Hot spot detection and management
    void register_portal_hot_spot(int portal_id, const PortalScalingInfo& info);
    void unregister_portal_hot_spot(int portal_id);
    void update_portal_hot_spot(int portal_id, const PortalScalingInfo& info);
    
    // Interaction detection
    AsyncTask<bool> detect_hot_spot_interaction(const core::SphericalCoords& cursor_position);
    std::vector<int> get_hot_spots_in_radius(const core::SphericalCoords& center, double radius);
    
    // Hot spot properties
    void set_sensitivity_threshold(double threshold);
    void set_interaction_radius(double radius);
    void enable_hot_spot_visualization(bool enabled);
    
    // Event system
    using HotSpotInteractionCallback = std::function<void(int portal_id, double intensity)>;
    void set_interaction_callback(HotSpotInteractionCallback callback);

private:
    struct HotSpot {
        int portal_id;
        core::SphericalCoords position;
        double intensity;
        double radius;
        std::chrono::steady_clock::time_point last_update;
    };
    
    std::unordered_map<int, HotSpot> active_hot_spots_;
    double sensitivity_threshold_;
    double interaction_radius_;
    HotSpotInteractionCallback interaction_callback_;
    
    double calculate_interaction_intensity(
        const core::SphericalCoords& cursor,
        const HotSpot& hot_spot
    ) const;
};

// HSML development acceleration system
class HSMLDevelopmentAccelerator {
public:
    HSMLDevelopmentAccelerator();
    ~HSMLDevelopmentAccelerator();
    
    // Acceleration modes
    void set_acceleration_mode(P0RT4L5DevelopmentFramework::AccelerationMode mode);
    P0RT4L5DevelopmentFramework::AccelerationMode get_acceleration_mode() const;
    
    // Test suite execution
    AsyncTask<bool> run_full_scale_development_tests();
    AsyncTask<bool> run_rapid_prototyping_tests();
    AsyncTask<bool> run_stress_tests();
    AsyncTask<bool> run_performance_profiling();
    AsyncTask<bool> run_real_time_debugging();
    
    // HSML-specific acceleration
    AsyncTask<bool> accelerate_spherical_coordinate_tests();
    AsyncTask<bool> accelerate_solid_angle_calculations();
    AsyncTask<bool> accelerate_physics_simulation_tests();
    AsyncTask<bool> accelerate_rendering_pipeline_tests();
    
    // Code generation and testing
    AsyncTask<std::string> generate_test_hsml_code(const std::string& test_spec);
    AsyncTask<bool> validate_generated_code(const std::string& hsml_code);
    AsyncTask<std::string> optimize_hsml_performance(const std::string& hsml_code);

private:
    P0RT4L5DevelopmentFramework::AccelerationMode current_mode_;
    
    // Test execution engines
    AsyncTask<bool> execute_test_batch(const std::vector<std::string>& test_names);
    AsyncTask<std::string> generate_performance_report();
    AsyncTask<bool> validate_hsml_syntax(const std::string& code);
};

// Factory for creating P0RT4L5 instances
class P0RT4L5Factory {
public:
    static std::unique_ptr<P0RT4L5DevelopmentFramework> create_framework();
    static std::unique_ptr<P0RT4L5DevelopmentFramework> create_framework_with_browser(
        std::shared_ptr<P0rt3rBrowserEngine> browser
    );
    
    // Pre-configured setups
    static std::unique_ptr<P0RT4L5DevelopmentFramework> create_for_hsml_development();
    static std::unique_ptr<P0RT4L5DevelopmentFramework> create_for_performance_testing();
    static std::unique_ptr<P0RT4L5DevelopmentFramework> create_for_integration_testing();
};

} // namespace browser
} // namespace hsml