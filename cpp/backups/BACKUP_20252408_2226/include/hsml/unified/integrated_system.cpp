#pragma once

#include "../core/constexpr_spherical_coords.hpp"
#include "../core/state_tensor.hpp"
#include "../core/constexpr_solid_angle.hpp"
#include "../compute/heterogeneous_compute_engine.hpp"
#include "../rendering/renderer_base.hpp"
#include "../rendering/concurrent_spherical_renderer.hpp"
#include <vector>
#include <memory>
#include <unordered_map>
#include <string>

namespace hsml::unified {

// Forward declarations for SDT integration
namespace sdt {
    class State21D;
    class Node;
    class NodalNetwork;
    class SDTMath;
}

class IntegratedHSMLSystem {
private:
    // Core rendering components
    std::unique_ptr<hsml::rendering::concurrent_spherical_renderer<>> renderer_;
    std::unique_ptr<hsml::compute::heterogeneous_compute_engine<>> compute_engine_;
    
    // SDT physics components
    std::unique_ptr<sdt::NodalNetwork> sdt_network_;
    std::vector<std::shared_ptr<sdt::Node>> sdt_nodes_;
    
    // State management
    std::unordered_map<std::string, hsml::core::state_tensor<double>> state_cache_;
    std::vector<hsml::core::spherical_coords<double>> spatial_coordinates_;
    
    // Performance monitoring
    struct SystemMetrics {
        size_t total_nodes = 0;
        size_t active_renders = 0;
        double avg_frame_time = 0.0;
        double sdt_update_time = 0.0;
        size_t memory_usage = 0;
    } metrics_;
    
public:
    IntegratedHSMLSystem(size_t width = 1920, size_t height = 1080);
    ~IntegratedHSMLSystem();
    
    // System initialization
    bool initialize();
    void shutdown();
    
    // SDT Physics Integration
    void add_sdt_node(const hsml::core::spherical_coords<double>& position,
                      double mass = 1.0,
                      const std::string& type = "basic");
    
    void update_sdt_physics(double dt);
    void sync_sdt_to_rendering();
    
    // HSML Rendering Integration
    void render_frame();
    void set_viewport(size_t width, size_t height);
    void clear_scene();
    
    // Unified coordinate system
    hsml::core::spherical_coords<double> world_to_spherical(double x, double y, double z) const;
    hsml::core::vector3<double> spherical_to_world(const hsml::core::spherical_coords<double>& coord) const;
    
    // State management
    void update_state_tensor(const std::string& entity_id, const hsml::core::state_tensor<double>& state);
    hsml::core::state_tensor<double> get_state_tensor(const std::string& entity_id) const;
    
    // Performance and monitoring
    SystemMetrics get_metrics() const { return metrics_; }
    void reset_metrics();
    
    // Compute dispatch
    template<typename Operation>
    auto dispatch_compute(std::span<const hsml::core::spherical_coords<double>> coords) 
        -> std::future<std::vector<typename Operation::result_type>>;
    
    // Spatial queries
    std::vector<std::string> query_spatial_region(
        const hsml::core::spherical_coords<double>& center,
        double radius) const;
    
    double calculate_solid_angle_between(
        const hsml::core::spherical_coords<double>& a,
        const hsml::core::spherical_coords<double>& b) const;
};

class UnifiedDemoApplication {
private:
    std::unique_ptr<IntegratedHSMLSystem> system_;
    bool running_ = false;
    double simulation_time_ = 0.0;
    
    // Demo scenarios
    void setup_solar_system_demo();
    void setup_particle_physics_demo();
    void setup_electromagnetic_demo();
    void setup_performance_stress_test();
    
public:
    UnifiedDemoApplication();
    ~UnifiedDemoApplication();
    
    bool initialize();
    void run();
    void shutdown();
    
    // Demo controls
    void switch_demo(const std::string& demo_name);
    void toggle_pause();
    void reset_simulation();
    void set_time_scale(double scale);
    
    // Live parameter adjustment
    void adjust_physics_parameters(const std::unordered_map<std::string, double>& params);
    void save_simulation_state(const std::string& filename);
    void load_simulation_state(const std::string& filename);
};

// Utility functions for cross-system data conversion
namespace conversion_utils {
    hsml::core::state_tensor<double> sdt_state_to_tensor(const sdt::State21D& sdt_state);
    sdt::State21D tensor_to_sdt_state(const hsml::core::state_tensor<double>& tensor);
    
    hsml::core::spherical_coords<double> sdt_position_to_spherical(const sdt::SphericalCoord& sdt_pos);
    sdt::SphericalCoord spherical_to_sdt_position(const hsml::core::spherical_coords<double>& spherical);
    
    std::vector<hsml::rendering::rendered_fragment> render_sdt_network(
        const sdt::NodalNetwork& network,
        hsml::rendering::concurrent_spherical_renderer<>& renderer);
}

// Performance optimization utilities
namespace optimization {
    void enable_simd_optimizations();
    void configure_memory_pools(size_t node_pool_size, size_t render_pool_size);
    void setup_spatial_indexing(const std::vector<hsml::core::spherical_coords<double>>& positions);
    
    class PerformanceProfiler {
    public:
        void start_frame();
        void end_frame();
        void mark_sdt_update_start();
        void mark_sdt_update_end();
        void mark_render_start();
        void mark_render_end();
        
        struct FrameProfile {
            double total_time;
            double sdt_time;
            double render_time;
            double sync_time;
        };
        
        FrameProfile get_last_frame() const;
        std::vector<FrameProfile> get_history(size_t count) const;
    };
}

// Integration testing framework
namespace testing {
    class IntegrationValidator {
    public:
        bool validate_sdt_hsml_sync();
        bool validate_coordinate_conversions();
        bool validate_state_tensor_operations();
        bool validate_rendering_pipeline();
        bool validate_compute_dispatch();
        
        struct ValidationReport {
            bool overall_success;
            std::vector<std::string> passed_tests;
            std::vector<std::string> failed_tests;
            std::vector<std::string> warnings;
        };
        
        ValidationReport run_full_validation();
    };
    
    void run_integration_stress_test(size_t node_count, double duration);
    void benchmark_cross_system_performance();
}

}