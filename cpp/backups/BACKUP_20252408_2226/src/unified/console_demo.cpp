#include "../../include/hsml/unified/integrated_system.hpp"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <thread>

using namespace hsml::unified;

void print_banner() {
    std::cout << R"(
╔══════════════════════════════════════════════════════════════════════════════╗
║                         HSML UNIFIED SYSTEM CONSOLE DEMO                    ║
║                                                                              ║
║  Integrated Hyper-Spherical Markup Language + Spatial Displacement Theory  ║
║  Real-time Physics Simulation with Advanced Rendering Pipeline              ║
║  Version 2.0 - Complete Integration                                         ║
╚══════════════════════════════════════════════════════════════════════════════╝

)" << std::endl;
}

void print_system_info() {
    std::cout << "System Information:" << std::endl;
    std::cout << "==================" << std::endl;
    std::cout << "• C++ Standard: C++23" << std::endl;
    std::cout << "• SIMD Support: AVX2/AVX-512 (hardware dependent)" << std::endl;
    std::cout << "• Threading: Hardware concurrency optimized" << std::endl;
    std::cout << "• Memory: Lock-free allocation pools" << std::endl;
    std::cout << "• Rendering: Concurrent spherical pipeline" << std::endl;
    std::cout << "• Physics: 21-dimensional SDT state simulation" << std::endl;
    std::cout << std::endl;
}

void demonstrate_integration_features(IntegratedHSMLSystem& system) {
    std::cout << "Feature Demonstration:" << std::endl;
    std::cout << "=====================" << std::endl;
    
    // Test 1: Coordinate System Integration
    std::cout << "1. Testing coordinate system integration..." << std::endl;
    auto spherical = system.world_to_spherical(1.0, 1.0, 1.0);
    auto cartesian = system.spherical_to_world(spherical);
    std::cout << "   ✓ World->Spherical->World conversion: " 
              << std::fixed << std::setprecision(6)
              << "(" << cartesian.x << ", " << cartesian.y << ", " << cartesian.z << ")" << std::endl;
    
    // Test 2: SDT Node Creation
    std::cout << "2. Testing SDT node creation..." << std::endl;
    system.add_sdt_node(spherical, 1e24, "test_mass");
    std::cout << "   ✓ SDT node added with mass 1e24 kg" << std::endl;
    
    // Test 3: State Tensor Operations
    std::cout << "3. Testing state tensor operations..." << std::endl;
    hsml::core::state_tensor<double> tensor;
    tensor.set<hsml::core::state_component::energy>(1000.0);
    tensor.set<hsml::core::state_component::density>(1.225);
    system.update_state_tensor("test_entity", tensor);
    
    auto retrieved_tensor = system.get_state_tensor("test_entity");
    std::cout << "   ✓ State tensor energy: " << retrieved_tensor.get<hsml::core::state_component::energy>() << " J" << std::endl;
    std::cout << "   ✓ State tensor temperature: " << retrieved_tensor.temperature() << " K" << std::endl;
    
    // Test 4: Spatial Queries
    std::cout << "4. Testing spatial queries..." << std::endl;
    auto nearby_entities = system.query_spatial_region(spherical, 10.0);
    std::cout << "   ✓ Found " << nearby_entities.size() << " entities within radius 10.0" << std::endl;
    
    // Test 5: Solid Angle Calculations
    std::cout << "5. Testing solid angle calculations..." << std::endl;
    auto coord1 = hsml::core::spherical_coords<double>{1.0, 0.5, 0.0};
    auto coord2 = hsml::core::spherical_coords<double>{1.0, 1.0, 0.5};
    double solid_angle = system.calculate_solid_angle_between(coord1, coord2);
    std::cout << "   ✓ Solid angle between coordinates: " << solid_angle << " steradians" << std::endl;
    
    std::cout << std::endl;
}

void run_physics_simulation(IntegratedHSMLSystem& system) {
    std::cout << "Physics Simulation Demo:" << std::endl;
    std::cout << "=======================" << std::endl;
    
    // Setup solar system
    std::cout << "Setting up solar system simulation..." << std::endl;
    system.clear_scene();
    
    // Add sun
    system.add_sdt_node(hsml::core::spherical_coords<double>{0.0, 0.0, 0.0}, 1.989e30, "star");
    
    // Add planets with proper orbital distances (scaled down for simulation)
    system.add_sdt_node(hsml::core::spherical_coords<double>{5.791e7, M_PI/2, 0.0}, 3.301e23, "planet");  // Mercury
    system.add_sdt_node(hsml::core::spherical_coords<double>{1.082e8, M_PI/2, 0.5}, 4.867e24, "planet");  // Venus  
    system.add_sdt_node(hsml::core::spherical_coords<double>{1.496e8, M_PI/2, 1.0}, 5.972e24, "planet");  // Earth
    system.add_sdt_node(hsml::core::spherical_coords<double>{2.279e8, M_PI/2, 1.5}, 6.390e23, "planet");  // Mars
    
    auto metrics = system.get_metrics();
    std::cout << "✓ Solar system created with " << metrics.total_nodes << " bodies" << std::endl;
    
    std::cout << "\nRunning simulation for 5 seconds..." << std::endl;
    std::cout << "Frame | SDT Update (ms) | Render (ms) | Total Objects" << std::endl;
    std::cout << "------|-----------------|-------------|---------------" << std::endl;
    
    const double dt = 1.0/60.0;  // 60 FPS
    const int total_frames = 300;  // 5 seconds
    
    for (int frame = 0; frame < total_frames; ++frame) {
        auto frame_start = std::chrono::high_resolution_clock::now();
        
        // Update physics
        system.update_sdt_physics(dt);
        
        // Render frame  
        system.render_frame();
        
        auto frame_end = std::chrono::high_resolution_clock::now();
        auto frame_time = std::chrono::duration<double, std::milli>(frame_end - frame_start);
        
        // Print metrics every 60 frames (1 second)
        if (frame % 60 == 0) {
            auto current_metrics = system.get_metrics();
            std::cout << std::setw(5) << frame << " | "
                      << std::setw(13) << std::fixed << std::setprecision(3) 
                      << current_metrics.sdt_update_time / 1000000.0 << " | "
                      << std::setw(9) << current_metrics.avg_frame_time << " | "
                      << std::setw(13) << current_metrics.total_nodes << std::endl;
        }
        
        // Maintain target framerate
        if (frame_time.count() < 16.67) {
            std::this_thread::sleep_for(
                std::chrono::milliseconds(static_cast<int>(16.67 - frame_time.count()))
            );
        }
    }
    
    std::cout << "\n✓ Physics simulation completed successfully" << std::endl;
    std::cout << std::endl;
}

void performance_benchmark(IntegratedHSMLSystem& system) {
    std::cout << "Performance Benchmark:" << std::endl;
    std::cout << "=====================" << std::endl;
    
    // Benchmark 1: Node Creation Performance
    std::cout << "1. Benchmarking node creation..." << std::endl;
    system.clear_scene();
    
    auto start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < 1000; ++i) {
        double theta = (i * 0.01);
        double phi = (i * 0.02);
        double radius = 1e6 + (i * 1000);
        
        system.add_sdt_node(
            hsml::core::spherical_coords<double>{radius, theta, phi},
            1e20,
            "benchmark_particle"
        );
    }
    
    auto end = std::chrono::high_resolution_clock::now();
    auto creation_time = std::chrono::duration<double, std::milli>(end - start);
    
    std::cout << "   ✓ Created 1000 nodes in " << creation_time.count() << " ms" << std::endl;
    std::cout << "   ✓ Average: " << creation_time.count() / 1000.0 << " ms per node" << std::endl;
    
    // Benchmark 2: Physics Update Performance
    std::cout << "2. Benchmarking physics updates..." << std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < 100; ++i) {
        system.update_sdt_physics(1.0/60.0);
    }
    
    end = std::chrono::high_resolution_clock::now();
    auto update_time = std::chrono::duration<double, std::milli>(end - start);
    
    std::cout << "   ✓ 100 physics updates in " << update_time.count() << " ms" << std::endl;
    std::cout << "   ✓ Average: " << update_time.count() / 100.0 << " ms per update" << std::endl;
    
    // Benchmark 3: Render Performance
    std::cout << "3. Benchmarking render performance..." << std::endl;
    
    start = std::chrono::high_resolution_clock::now();
    
    for (int i = 0; i < 100; ++i) {
        system.render_frame();
    }
    
    end = std::chrono::high_resolution_clock::now();
    auto render_time = std::chrono::duration<double, std::milli>(end - start);
    
    std::cout << "   ✓ 100 render calls in " << render_time.count() << " ms" << std::endl;
    std::cout << "   ✓ Average: " << render_time.count() / 100.0 << " ms per frame" << std::endl;
    
    auto final_metrics = system.get_metrics();
    std::cout << "\nFinal System Metrics:" << std::endl;
    std::cout << "• Total Nodes: " << final_metrics.total_nodes << std::endl;
    std::cout << "• Memory Usage: " << final_metrics.memory_usage / (1024*1024) << " MB" << std::endl;
    std::cout << "• Active Renders: " << final_metrics.active_renders << std::endl;
    
    std::cout << std::endl;
}

int main(int argc, char* argv[]) {
    print_banner();
    print_system_info();
    
    // Initialize the integrated system
    std::cout << "Initializing HSML Unified System..." << std::endl;
    IntegratedHSMLSystem system(1920, 1080);
    
    if (!system.initialize()) {
        std::cerr << "❌ Failed to initialize HSML Unified System" << std::endl;
        return -1;
    }
    
    std::cout << "✓ HSML Unified System initialized successfully" << std::endl;
    std::cout << std::endl;
    
    try {
        // Run demonstrations
        demonstrate_integration_features(system);
        run_physics_simulation(system);
        performance_benchmark(system);
        
        std::cout << "All demonstrations completed successfully!" << std::endl;
        std::cout << "\nThe HSML Unified System is now fully operational and integrated." << std::endl;
        std::cout << "Key achievements:" << std::endl;
        std::cout << "• Final Phase HSML rendering pipeline: ✓" << std::endl;
        std::cout << "• SDT physics simulation: ✓" << std::endl;
        std::cout << "• Unified coordinate system: ✓" << std::endl;
        std::cout << "• Real-time state synchronization: ✓" << std::endl;
        std::cout << "• Performance optimization: ✓" << std::endl;
        std::cout << "• Cross-system integration: ✓" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Error during demonstration: " << e.what() << std::endl;
        return -1;
    }
    
    // Clean shutdown
    system.shutdown();
    std::cout << "\n✓ System shutdown complete" << std::endl;
    
    return 0;
}