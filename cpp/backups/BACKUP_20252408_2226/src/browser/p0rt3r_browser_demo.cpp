/**
 * P0rt3r Browser Demo Implementation - C++20
 * Demonstrates the revolutionary native HSML browser in action
 * Shows integration of all major components with example usage
 */

#include <iostream>
#include <memory>
#include <string>
#include <chrono>
#include <thread>

// P0rt3r Browser Components
#include "hsml/browser/p0rt3r_browser_engine.h"
#include "hsml/browser/spatial_document_parser.h"
#include "hsml/browser/spatial_layout_engine.h"
#include "hsml/browser/navigation_manager.h"
#include "hsml/browser/javascript_engine.h"

// Existing HSML Foundation
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"

using namespace p0rt3r;
using namespace hsml::core;

// Demo function prototypes
void demonstrate_browser_initialization();
void demonstrate_document_loading();
void demonstrate_spatial_navigation();
void demonstrate_orbital_browsing();
void demonstrate_javascript_integration();
void demonstrate_performance_metrics();

/**
 * Main P0rt3r Browser Demo
 */
int main() {
    std::cout << "🌐 P0rt3r Browser Demo - World's First Native HSML Browser\n";
    std::cout << "========================================================\n\n";
    
    try {
        // Demonstrate core browser functionality
        demonstrate_browser_initialization();
        demonstrate_document_loading();
        demonstrate_spatial_navigation();
        demonstrate_orbital_browsing();
        demonstrate_javascript_integration();
        demonstrate_performance_metrics();
        
        std::cout << "\n✅ P0rt3r Browser Demo completed successfully!\n";
        std::cout << "🚀 Ready for revolutionary 3D web browsing experience!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Demo error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

/**
 * Demonstrate P0rt3r browser engine initialization
 */
void demonstrate_browser_initialization() {
    std::cout << "🔧 Initializing P0rt3r Browser Engine...\n";
    
    // Create browser engine with configuration
    auto browser_engine = engine::create_p0rt3r_engine<GraphicsAPI::OPENGL_4_5>();
    
    // Configure browser settings
    engine::P0rt3rBrowserEngine<>::BrowserConfig config;
    config.default_viewer_distance = 650.0;  // 650mm from monitor
    config.max_concurrent_documents = 16;
    config.cache_size_mb = 512;
    config.enable_javascript = true;
    config.enable_webgl_fallback = true;
    
    // Initialize the browser engine
    auto init_result = browser_engine->initialize(config);
    if (!init_result) {
        throw std::runtime_error("Failed to initialize browser engine: " + init_result.error());
    }
    
    std::cout << "  ✓ Browser engine initialized successfully\n";
    std::cout << "  ✓ Graphics API: OpenGL 4.5\n";
    std::cout << "  ✓ JavaScript engine: V8 with spatial bindings\n";
    std::cout << "  ✓ Viewer distance: " << config.default_viewer_distance << "mm\n\n";
}

/**
 * Demonstrate spatial document loading and parsing
 */
void demonstrate_document_loading() {
    std::cout << "📄 Loading Spatial HSML Documents...\n";
    
    // Create document parser
    parser::SpatialDocumentParser::ParserConfig parser_config;
    parser_config.enable_javascript = true;
    parser_config.enable_css_parsing = true;
    parser_config.default_viewer_distance = 650.0;
    
    auto document_parser = std::make_unique<parser::SpatialDocumentParser>(parser_config);
    
    // Example HSML document content with spatial elements
    const std::string sample_hsml = R"hsml(
        <hsml>
            <head>
                <title>Revolutionary 3D Web Page</title>
                <style type="text/csss">
                    .orbital-element {
                        position: spherical(500, 0.5, 1.0);
                        radius: 50;
                        color: rgba(0, 150, 255, 0.8);
                        matter-state: solid;
                        orbital-velocity: 0.1;
                    }
                    
                    .depth-layer {
                        position: spherical(800, 1.0, 2.0);
                        display: depth-layer;
                        interaction-radius: 100;
                    }
                </style>
            </head>
            <body origin="spherical(650, 0, 0)" radius="2000">
                <bubble id="welcome" class="orbital-element">
                    <h1>Welcome to P0rt3r Browser!</h1>
                    <p>Experience the web in revolutionary 3D space</p>
                </bubble>
                
                <bubble id="navigation" class="depth-layer">
                    <nav>
                        <button onclick="navigator.spatial.teleportTo({r: 1000, theta: 0.5, phi: 1.5})">
                            Teleport to Galaxy View
                        </button>
                        <button onclick="navigator.spatial.diveToDepth(300)">
                            Dive Deeper
                        </button>
                    </nav>
                </bubble>
                
                <presence id="physics-demo" position="spherical(400, 1.2, 0.8)" mass="2.0">
                    <script type="text/javascript">
                        // Apply gravitational field
                        physics.setGravityField({
                            center: new SphericalCoords(0, 0, 0),
                            strength: 0.1
                        });
                        
                        // Enable orbital motion
                        this.enableOrbitalMotion(true);
                    </script>
                </presence>
            </body>
        </hsml>
    )hsml";
    
    // Parse the document
    auto parse_result = document_parser->parse_document(sample_hsml, "demo://spatial-webpage");
    if (!parse_result) {
        throw std::runtime_error("Failed to parse document: " + parse_result.error());
    }
    
    auto document = std::move(parse_result.value());
    
    std::cout << "  ✓ Parsed HSML document successfully\n";
    std::cout << "  ✓ Document URI: " << document->document_uri() << "\n";
    std::cout << "  ✓ Document origin: (" << document->document_origin().r() << ", " 
              << document->document_origin().theta() << ", " << document->document_origin().phi() << ")\n";
    std::cout << "  ✓ Found " << document->dependencies().size() << " resource dependencies\n";
    std::cout << "  ✓ Extracted " << document->inline_scripts().size() << " inline scripts\n";
    std::cout << "  ✓ Loaded " << document->stylesheets().size() << " CSSS stylesheets\n\n";
}

/**
 * Demonstrate revolutionary spatial navigation
 */
void demonstrate_spatial_navigation() {
    std::cout << "🧭 Demonstrating Spatial Navigation...\n";
    
    // Initialize navigation manager
    SphericalCoords initial_position{650.0, 0.0, 0.0};  // 650mm from origin
    auto nav_manager = std::make_unique<navigation::NavigationManager>(initial_position);
    
    // Test basic navigation
    SphericalCoords target_position{1000.0, 0.5, 1.2};
    auto nav_animation = nav_manager->navigate_to(target_position);
    
    std::cout << "  ✓ Initiated navigation to (" << target_position.r() << ", " 
              << target_position.theta() << ", " << target_position.phi() << ")\n";
    
    // Test depth diving
    auto dive_animation = nav_manager->depth_controller().dive_deeper(300.0);
    std::cout << "  ✓ Diving to depth level: 300 units\n";
    
    // Test spatial teleportation
    SphericalCoords teleport_destination{2000.0, 1.5, 2.8};
    auto teleport_animation = nav_manager->teleportation_system().teleport_to(
        teleport_destination, navigation::SpatialTeleportationSystem::TeleportationType::ANIMATED);
    
    std::cout << "  ✓ Teleporting to distant coordinates with smooth animation\n";
    
    // Create spatial bookmark
    auto bookmark_id = nav_manager->teleportation_system().create_spatial_bookmark(
        teleport_destination, "Galaxy Overview", "Perfect view of the cosmic web structure");
    
    std::cout << "  ✓ Created spatial bookmark: " << bookmark_id << "\n\n";
}

/**
 * Demonstrate orbital browsing with tabs
 */
void demonstrate_orbital_browsing() {
    std::cout << "🪐 Demonstrating Orbital Tab Browsing...\n";
    
    // Initialize orbital tab manager
    SphericalCoords user_position{0.0, 0.0, 0.0};
    auto orbital_manager = std::make_unique<navigation::OrbitalTabManager>(user_position);
    
    // Create multiple tabs in orbital arrangement
    std::vector<std::string> demo_uris = {
        "hsml://github.com/spatial-web",
        "hsml://news.3d/breaking-spatial-web",
        "hsml://wikipedia.org/spherical-coordinates",
        "hsml://youtube.com/watch?v=3d-web-revolution"
    };
    
    std::vector<std::future<std::expected<std::string, std::string>>> tab_futures;
    
    for (const auto& uri : demo_uris) {
        auto tab_future = orbital_manager->create_tab(uri);
        tab_futures.push_back(std::move(tab_future));
    }
    
    // Wait for tabs to be created and collect IDs
    std::vector<std::string> tab_ids;
    for (auto& future : tab_futures) {
        auto result = future.get();
        if (result) {
            tab_ids.push_back(result.value());
            std::cout << "  ✓ Created orbital tab: " << result.value() << "\n";
        }
    }
    
    // Arrange tabs in different orbital patterns
    orbital_manager->arrange_tabs_spherically();
    std::cout << "  ✓ Arranged tabs in spherical orbit\n";
    
    // Rotate between tabs
    if (!tab_ids.empty()) {
        auto rotation_animation = orbital_manager->rotate_to_tab(tab_ids[1]);
        std::cout << "  ✓ Rotating to tab: " << tab_ids[1] << "\n";
        
        auto focus_animation = orbital_manager->bring_tab_to_focus(tab_ids[1]);
        std::cout << "  ✓ Bringing tab into focus with smooth animation\n";
    }
    
    // Demonstrate helical arrangement
    orbital_manager->arrange_tabs_helically();
    std::cout << "  ✓ Rearranged tabs in helical pattern\n\n";
}

/**
 * Demonstrate JavaScript integration with spatial APIs
 */
void demonstrate_javascript_integration() {
    std::cout << "🔧 Demonstrating JavaScript Integration...\n";
    
    // Initialize JavaScript engine
    javascript::SpatialJavaScriptEngine::EngineConfig js_config;
    js_config.max_heap_size_mb = 256;
    js_config.enable_debugging = true;
    js_config.enable_profiling = true;
    
    auto js_engine = std::make_unique<javascript::SpatialJavaScriptEngine>(js_config);
    
    auto init_result = js_engine->initialize();
    if (!init_result) {
        throw std::runtime_error("Failed to initialize JavaScript engine: " + init_result.error());
    }
    
    // Expose spatial APIs
    js_engine->expose_spherical_coordinates_api();
    js_engine->expose_solid_angle_api();
    js_engine->expose_spatial_navigation_api();
    js_engine->expose_physics_simulation_api();
    
    std::cout << "  ✓ JavaScript engine initialized with V8\n";
    std::cout << "  ✓ Spatial APIs exposed to JavaScript\n";
    
    // Create JavaScript context
    auto context_result = js_engine->create_context("demo://javascript-integration");
    if (!context_result) {
        throw std::runtime_error("Failed to create JavaScript context: " + context_result.error());
    }
    
    std::string context_id = context_result.value();
    std::cout << "  ✓ Created JavaScript context: " << context_id << "\n";
    
    // Test spatial coordinate manipulation in JavaScript
    const std::string spatial_js_code = R"js(
        // Create spherical coordinates
        const origin = new SphericalCoords(0, 0, 0);
        const target = new SphericalCoords(1000, Math.PI/4, Math.PI/2);
        
        // Calculate distance
        const distance = origin.distanceTo(target);
        console.log('Distance between points:', distance);
        
        // Interpolate between positions
        const midpoint = origin.interpolateTo(target, 0.5);
        console.log('Midpoint coordinates:', midpoint.r, midpoint.theta, midpoint.phi);
        
        // Create solid angle
        const viewCone = new SolidAngle(Math.PI/3, 0, Math.PI/2, 0, Math.PI);
        console.log('View cone solid angle:', viewCone.omega);
        
        // Test navigation API
        navigator.spatial.navigateTo(target).then(() => {
            console.log('Navigation completed successfully');
        });
        
        // Test physics simulation
        physics.setGravityField({
            center: origin,
            strength: 9.81
        });
        
        physics.applyForce('demo-element', {
            force: new SphericalCoords(10, 0, 0),
            duration: 1000
        });
        
        'JavaScript spatial integration test completed';
    )js";
    
    // Execute the spatial JavaScript code
    auto exec_result = js_engine->execute_script(spatial_js_code, context_id, "spatial-demo.js");
    if (exec_result) {
        std::cout << "  ✓ JavaScript execution result: " << exec_result.value() << "\n";
    } else {
        std::cout << "  ❌ JavaScript execution error: " << exec_result.error() << "\n";
    }
    
    // Get JavaScript engine metrics
    auto js_metrics = js_engine->get_metrics();
    std::cout << "  ✓ Scripts executed: " << js_metrics.total_scripts_executed << "\n";
    std::cout << "  ✓ Memory usage: " << js_metrics.memory_usage_mb << " MB\n";
    std::cout << "  ✓ Active contexts: " << js_metrics.active_contexts << "\n\n";
}

/**
 * Demonstrate performance metrics and monitoring
 */
void demonstrate_performance_metrics() {
    std::cout << "📊 Performance Metrics and Monitoring...\n";
    
    // Simulate some browser operations
    std::this_thread::sleep_for(std::chrono::milliseconds(100));
    
    // Example performance metrics (would be real data from browser engine)
    engine::P0rt3rBrowserEngine<>::BrowserMetrics browser_metrics;
    browser_metrics.documents_loaded = 5;
    browser_metrics.spatial_objects_rendered = 1247;
    browser_metrics.average_frame_time_ms = 12.3;
    browser_metrics.navigation_latency_ms = 85.2;
    browser_metrics.memory_usage_mb = 342;
    browser_metrics.active_javascript_contexts = 3;
    browser_metrics.gpu_utilization_percentage = 67.8;
    browser_metrics.gpu_memory_used_mb = 128;
    browser_metrics.triangles_rendered_per_frame = 45680;
    browser_metrics.solid_angle_coverage_percentage = 23.5;
    
    std::cout << "  📈 Browser Performance Metrics:\n";
    std::cout << "    • Documents loaded: " << browser_metrics.documents_loaded << "\n";
    std::cout << "    • Spatial objects rendered: " << browser_metrics.spatial_objects_rendered << "\n";
    std::cout << "    • Average frame time: " << browser_metrics.average_frame_time_ms << " ms\n";
    std::cout << "    • Navigation latency: " << browser_metrics.navigation_latency_ms << " ms\n";
    std::cout << "    • Memory usage: " << browser_metrics.memory_usage_mb << " MB\n";
    std::cout << "    • JavaScript contexts: " << browser_metrics.active_javascript_contexts << "\n";
    std::cout << "    • GPU utilization: " << browser_metrics.gpu_utilization_percentage << "%\n";
    std::cout << "    • GPU memory: " << browser_metrics.gpu_memory_used_mb << " MB\n";
    std::cout << "    • Triangles per frame: " << browser_metrics.triangles_rendered_per_frame << "\n";
    std::cout << "    • Viewport coverage: " << browser_metrics.solid_angle_coverage_percentage << "%\n";
    
    // Performance validation
    if (browser_metrics.average_frame_time_ms < 16.67) {
        std::cout << "  ✅ Excellent! Maintaining 60+ FPS\n";
    } else if (browser_metrics.average_frame_time_ms < 33.33) {
        std::cout << "  ⚠️  Good performance at 30+ FPS\n";
    } else {
        std::cout << "  ❌ Performance optimization needed\n";
    }
    
    if (browser_metrics.navigation_latency_ms < 100.0) {
        std::cout << "  ✅ Excellent navigation responsiveness\n";
    } else {
        std::cout << "  ⚠️  Navigation latency could be improved\n";
    }
    
    if (browser_metrics.gpu_utilization_percentage > 50.0 && browser_metrics.gpu_utilization_percentage < 90.0) {
        std::cout << "  ✅ Optimal GPU utilization\n";
    } else if (browser_metrics.gpu_utilization_percentage >= 90.0) {
        std::cout << "  ⚠️  High GPU usage - consider quality reduction\n";
    } else {
        std::cout << "  ℹ️  GPU has additional capacity available\n";
    }
    
    std::cout << "\n";
}

/**
 * Demonstrate revolutionary UI paradigms
 */
void demonstrate_revolutionary_ui() {
    std::cout << "🎭 Revolutionary UI Paradigms...\n";
    
    // This would demonstrate:
    // - Physics-based window management
    // - Gesture-controlled navigation  
    // - Voice-activated spatial commands
    // - Eye-tracking viewport control
    // - Haptic feedback for spatial interactions
    
    std::cout << "  🪟 Physics-based windows: Windows behave like real 3D objects\n";
    std::cout << "  👋 Gesture navigation: Hand tracking for spatial manipulation\n";
    std::cout << "  🗣️  Voice commands: 'Navigate to bookmark Galaxy View'\n";
    std::cout << "  👁️  Eye tracking: Viewport follows natural gaze patterns\n";
    std::cout << "  🎮 Haptic feedback: Feel spatial interactions and collisions\n\n";
}

/**
 * Additional utility functions for demo
 */
namespace demo_utils {
    void print_separator(const std::string& title) {
        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << "  " << title << "\n";
        std::cout << std::string(60, '=') << "\n\n";
    }
    
    void simulate_user_interaction() {
        std::cout << "  👤 Simulating user interaction...\n";
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
    }
    
    SphericalCoords generate_random_spatial_position() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        static std::uniform_real_distribution<double> r_dist(100.0, 2000.0);
        static std::uniform_real_distribution<double> theta_dist(0.0, M_PI);
        static std::uniform_real_distribution<double> phi_dist(0.0, 2.0 * M_PI);
        
        return SphericalCoords{r_dist(gen), theta_dist(gen), phi_dist(gen)};
    }
}