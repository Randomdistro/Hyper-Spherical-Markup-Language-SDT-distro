#include "../../include/hsml/unified/integrated_system.hpp"
#include "../core/sdt_core.cpp"  // Include SDT implementation
#include <chrono>
#include <thread>
#include <algorithm>
#include <iostream>
#include <fstream>

namespace hsml::unified {

IntegratedHSMLSystem::IntegratedHSMLSystem(size_t width, size_t height) 
    : renderer_(std::make_unique<hsml::rendering::concurrent_spherical_renderer<>>(width, height))
    , compute_engine_(std::make_unique<hsml::compute::heterogeneous_compute_engine<>>())
    , sdt_network_(std::make_unique<HSML::NodalNetwork>()) {
}

IntegratedHSMLSystem::~IntegratedHSMLSystem() {
    shutdown();
}

bool IntegratedHSMLSystem::initialize() {
    try {
        // Initialize compute engine
        std::cout << "Initializing compute engine: " << compute_engine_->backend_name() << std::endl;
        
        // Initialize SDT network
        sdt_network_ = std::make_unique<HSML::NodalNetwork>();
        
        // Clear initial state
        spatial_coordinates_.clear();
        state_cache_.clear();
        
        metrics_ = SystemMetrics{};
        
        std::cout << "HSML Integrated System initialized successfully" << std::endl;
        return true;
    } catch (const std::exception& e) {
        std::cerr << "Failed to initialize HSML system: " << e.what() << std::endl;
        return false;
    }
}

void IntegratedHSMLSystem::shutdown() {
    if (renderer_) {
        renderer_.reset();
    }
    if (compute_engine_) {
        compute_engine_.reset();
    }
    if (sdt_network_) {
        sdt_network_.reset();
    }
    sdt_nodes_.clear();
    spatial_coordinates_.clear();
    state_cache_.clear();
}

void IntegratedHSMLSystem::add_sdt_node(const hsml::core::spherical_coords<double>& position,
                                        double mass,
                                        const std::string& type) {
    // Create SDT node
    auto sdt_node = std::make_shared<HSML::Node>();
    sdt_node->mass = mass;
    sdt_node->shape_type = type;
    
    // Convert HSML spherical coords to SDT spherical coords
    sdt_node->position = HSML::SphericalCoord(
        0.0,  // axis
        position.theta() * 180.0 / M_PI,  // polar in degrees
        position.phi() * 180.0 / M_PI,    // azimuth in degrees
        position.radius()                  // radius
    );
    
    // Add to SDT network
    sdt_network_->add_node(sdt_node);
    sdt_nodes_.push_back(sdt_node);
    
    // Add to spatial coordinates for HSML
    spatial_coordinates_.push_back(position);
    
    // Create corresponding state tensor
    hsml::core::state_tensor<double> tensor;
    tensor.set<hsml::core::state_component::density>(mass);
    tensor.set<hsml::core::state_component::energy>(mass * 9e16);  // E=mc²
    state_cache_[sdt_node->id] = tensor;
    
    metrics_.total_nodes++;
}

void IntegratedHSMLSystem::update_sdt_physics(double dt) {
    auto start = std::chrono::high_resolution_clock::now();
    
    // Update SDT network state
    sdt_network_->update_network_state(dt);
    
    // Sync updated positions back to HSML coordinates
    spatial_coordinates_.clear();
    for (const auto& node : sdt_nodes_) {
        auto cartesian = node->position.to_cartesian();
        auto spherical = hsml::core::spherical_coords<double>::from_cartesian(
            hsml::core::vector3<double>{cartesian[0], cartesian[1], cartesian[2]}
        );
        spatial_coordinates_.push_back(spherical);
        
        // Update state tensor
        auto& tensor = state_cache_[node->id];
        tensor.set<hsml::core::state_component::energy>(node->state.calculate_total_energy());
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    metrics_.sdt_update_time = std::chrono::duration<double, std::nano>(end - start).count();
}

void IntegratedHSMLSystem::sync_sdt_to_rendering() {
    // Create rendering elements from SDT nodes
    for (size_t i = 0; i < sdt_nodes_.size(); ++i) {
        // This would integrate with the concurrent renderer
        // For now, we just ensure spatial coordinates are synced
        if (i < spatial_coordinates_.size()) {
            // The rendering system can access spatial_coordinates_
        }
    }
}

void IntegratedHSMLSystem::render_frame() {
    auto start = std::chrono::high_resolution_clock::now();
    
    // Clear framebuffer
    renderer_->clear_framebuffer();
    
    // Sync SDT state to rendering
    sync_sdt_to_rendering();
    
    // For each spatial coordinate, we could submit rendering commands
    // This demonstrates the integration pattern
    
    auto end = std::chrono::high_resolution_clock::now();
    metrics_.avg_frame_time = std::chrono::duration<double, std::milli>(end - start).count();
    metrics_.active_renders++;
}

void IntegratedHSMLSystem::set_viewport(size_t width, size_t height) {
    // This would configure the renderer viewport
    // renderer_->set_viewport(width, height);
}

void IntegratedHSMLSystem::clear_scene() {
    sdt_nodes_.clear();
    spatial_coordinates_.clear();
    state_cache_.clear();
    
    // Create new SDT network
    sdt_network_ = std::make_unique<HSML::NodalNetwork>();
    
    metrics_.total_nodes = 0;
}

hsml::core::spherical_coords<double> IntegratedHSMLSystem::world_to_spherical(double x, double y, double z) const {
    return hsml::core::spherical_coords<double>::from_cartesian(
        hsml::core::vector3<double>{x, y, z}
    );
}

hsml::core::vector3<double> IntegratedHSMLSystem::spherical_to_world(
    const hsml::core::spherical_coords<double>& coord) const {
    return coord.to_cartesian();
}

void IntegratedHSMLSystem::update_state_tensor(const std::string& entity_id, 
                                               const hsml::core::state_tensor<double>& state) {
    state_cache_[entity_id] = state;
}

hsml::core::state_tensor<double> IntegratedHSMLSystem::get_state_tensor(const std::string& entity_id) const {
    auto it = state_cache_.find(entity_id);
    if (it != state_cache_.end()) {
        return it->second;
    }
    return hsml::core::state_tensor<double>::zero();
}

void IntegratedHSMLSystem::reset_metrics() {
    metrics_ = SystemMetrics{};
}

std::vector<std::string> IntegratedHSMLSystem::query_spatial_region(
    const hsml::core::spherical_coords<double>& center,
    double radius) const {
    
    std::vector<std::string> results;
    
    for (const auto& node : sdt_nodes_) {
        auto node_cartesian = node->position.to_cartesian();
        auto node_spherical = hsml::core::spherical_coords<double>::from_cartesian(
            hsml::core::vector3<double>{node_cartesian[0], node_cartesian[1], node_cartesian[2]}
        );
        
        // Calculate distance
        auto center_cart = center.to_cartesian();
        auto node_cart = node_spherical.to_cartesian();
        
        double dx = center_cart.x - node_cart.x;
        double dy = center_cart.y - node_cart.y;
        double dz = center_cart.z - node_cart.z;
        double distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        
        if (distance <= radius) {
            results.push_back(node->id);
        }
    }
    
    return results;
}

double IntegratedHSMLSystem::calculate_solid_angle_between(
    const hsml::core::spherical_coords<double>& a,
    const hsml::core::spherical_coords<double>& b) const {
    
    using engine = hsml::core::solid_angle_engine<double>;
    
    double theta_diff = std::abs(a.theta() - b.theta());
    double phi_diff = std::abs(a.phi() - b.phi());
    
    return engine::spherical_wedge_solid_angle(0, theta_diff, 0, phi_diff);
}

// UnifiedDemoApplication Implementation
UnifiedDemoApplication::UnifiedDemoApplication() 
    : system_(std::make_unique<IntegratedHSMLSystem>()) {
}

UnifiedDemoApplication::~UnifiedDemoApplication() {
    shutdown();
}

bool UnifiedDemoApplication::initialize() {
    if (!system_->initialize()) {
        return false;
    }
    
    std::cout << "Unified Demo Application initialized" << std::endl;
    return true;
}

void UnifiedDemoApplication::run() {
    if (!initialize()) {
        std::cerr << "Failed to initialize demo application" << std::endl;
        return;
    }
    
    running_ = true;
    const double dt = 1.0/60.0;  // 60 FPS
    
    std::cout << "Running unified HSML/SDT demonstration..." << std::endl;
    std::cout << "Starting with solar system demo..." << std::endl;
    
    setup_solar_system_demo();
    
    while (running_) {
        auto frame_start = std::chrono::high_resolution_clock::now();
        
        // Update physics
        system_->update_sdt_physics(dt);
        
        // Render frame
        system_->render_frame();
        
        simulation_time_ += dt;
        
        // Print metrics every 60 frames
        static int frame_count = 0;
        if (++frame_count % 60 == 0) {
            auto metrics = system_->get_metrics();
            std::cout << "Frame " << frame_count 
                      << " - Nodes: " << metrics.total_nodes
                      << ", SDT Update: " << metrics.sdt_update_time/1000000.0 << "ms"
                      << ", Frame Time: " << metrics.avg_frame_time << "ms" << std::endl;
        }
        
        auto frame_end = std::chrono::high_resolution_clock::now();
        auto frame_time = std::chrono::duration<double, std::milli>(frame_end - frame_start);
        
        // Maintain 60 FPS
        if (frame_time.count() < 16.67) {
            std::this_thread::sleep_for(
                std::chrono::milliseconds(static_cast<int>(16.67 - frame_time.count()))
            );
        }
        
        // Run for 10 seconds then stop
        if (simulation_time_ > 10.0) {
            running_ = false;
        }
    }
    
    std::cout << "Demo completed. Total simulation time: " << simulation_time_ << "s" << std::endl;
}

void UnifiedDemoApplication::shutdown() {
    running_ = false;
    if (system_) {
        system_->shutdown();
    }
}

void UnifiedDemoApplication::setup_solar_system_demo() {
    system_->clear_scene();
    
    // Add sun
    system_->add_sdt_node(hsml::core::spherical_coords<double>{0.0, 0.0, 0.0}, 1.989e30, "star");
    
    // Add planets (simplified)
    system_->add_sdt_node(hsml::core::spherical_coords<double>{5.791e10, M_PI/2, 0.0}, 3.301e23, "planet");  // Mercury
    system_->add_sdt_node(hsml::core::spherical_coords<double>{1.082e11, M_PI/2, 0.5}, 4.867e24, "planet");  // Venus
    system_->add_sdt_node(hsml::core::spherical_coords<double>{1.496e11, M_PI/2, 1.0}, 5.972e24, "planet");  // Earth
    system_->add_sdt_node(hsml::core::spherical_coords<double>{2.279e11, M_PI/2, 1.5}, 6.390e23, "planet");  // Mars
    
    std::cout << "Solar system demo setup complete with " << system_->get_metrics().total_nodes << " bodies" << std::endl;
}

void UnifiedDemoApplication::setup_particle_physics_demo() {
    system_->clear_scene();
    
    // Create a grid of particles
    for (int i = -5; i <= 5; ++i) {
        for (int j = -5; j <= 5; ++j) {
            for (int k = -5; k <= 5; ++k) {
                double x = i * 1e-10;
                double y = j * 1e-10; 
                double z = k * 1e-10;
                auto coord = system_->world_to_spherical(x, y, z);
                system_->add_sdt_node(coord, 9.109e-31, "electron");  // Electron mass
            }
        }
    }
    
    std::cout << "Particle physics demo setup complete with " << system_->get_metrics().total_nodes << " particles" << std::endl;
}

void UnifiedDemoApplication::switch_demo(const std::string& demo_name) {
    if (demo_name == "solar_system") {
        setup_solar_system_demo();
    } else if (demo_name == "particle_physics") {
        setup_particle_physics_demo();
    } else if (demo_name == "stress_test") {
        setup_performance_stress_test();
    }
}

void UnifiedDemoApplication::setup_performance_stress_test() {
    system_->clear_scene();
    
    // Create many random particles for stress testing
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> pos_dist(-100.0, 100.0);
    std::uniform_real_distribution<> mass_dist(1e-30, 1e-20);
    
    for (int i = 0; i < 1000; ++i) {
        double x = pos_dist(gen);
        double y = pos_dist(gen);
        double z = pos_dist(gen);
        double mass = mass_dist(gen);
        
        auto coord = system_->world_to_spherical(x, y, z);
        system_->add_sdt_node(coord, mass, "test_particle");
    }
    
    std::cout << "Stress test setup complete with " << system_->get_metrics().total_nodes << " particles" << std::endl;
}

// Conversion utilities implementation
namespace conversion_utils {
    hsml::core::state_tensor<double> sdt_state_to_tensor(const HSML::State21D& sdt_state) {
        hsml::core::state_tensor<double> tensor;
        
        // Map SDT 21D state to 8D tensor (simplified mapping)
        auto energy_level = sdt_state.get_level(HSML::SDTLevel::ENERGY);
        auto space_level = sdt_state.get_level(HSML::SDTLevel::SPACE_3D);
        
        if (!energy_level.empty()) {
            tensor.set<hsml::core::state_component::energy>(energy_level[0]);
        }
        if (space_level.size() >= 3) {
            tensor.set<hsml::core::state_component::velocity>(
                std::sqrt(space_level[0]*space_level[0] + space_level[1]*space_level[1] + space_level[2]*space_level[2])
            );
        }
        
        return tensor;
    }
}

} // namespace hsml::unified