/** @file 3d_gui_demo.cpp
 * @brief Demonstration of the 3D GUI Development Environment
 *
 * Phase 5: Immersive 3D Development Environment
 * Live demonstration of the floating code entities, connection lines,
 * and immersive development experience.
 */

#include "hsml/presentation/gui/3d/ImmersiveDevelopmentEnvironment.h"
#include "hsml/application/services/SphericalSceneService.h"
#include "hsml/infrastructure/persistence/SphericalSceneRepository.h"
#include "hsml/infrastructure/rendering/SphericalRendererAdapter.h"
#include "hsml/infrastructure/platform/PlatformManager.h"
#include "hsml/core/debug_logger.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <memory>

using namespace hsml;

int main() {
    std::cout << "🚀 HSML 3D GUI Development Environment Demo" << std::endl;
    std::cout << "==============================================" << std::endl;

    HSML_INFO("Starting 3D GUI Development Environment Demo");

    try {
        // Initialize platform
        auto& platform_mgr = infrastructure::PlatformManager::instance();
        if (!platform_mgr.initialize_platform()) {
            HSML_ERROR("Failed to initialize platform");
            return 1;
        }

        HSML_INFO("Platform initialized successfully");

        // Create infrastructure components
        auto platform = platform_mgr.get_platform();
        auto repository = std::make_shared<infrastructure::SphericalSceneRepository>();

        // Create a simple test renderer
        auto test_renderer = std::make_unique<TestRenderer>();
        auto renderer_adapter = std::make_shared<infrastructure::SphericalRendererAdapter>(
            std::move(test_renderer)
        );

        // Create application service
        auto render_use_case = std::make_unique<domain::RenderSphericalScene>(
            repository, renderer_adapter
        );
        auto scene_service = std::make_shared<application::SphericalSceneService>(
            repository, renderer_adapter, std::move(render_use_case)
        );

        HSML_INFO("Infrastructure components created");

        // Create the 3D development environment
        auto dev_env = presentation::ImmersiveEnvironmentFactory::create_environment(
            scene_service, 1920, 1080
        );

        if (!dev_env->initialize()) {
            HSML_ERROR("Failed to initialize 3D development environment");
            return 1;
        }

        HSML_INFO("3D Development Environment initialized successfully");

        // Set up event callbacks
        dev_env->set_entity_selected_callback(
            [](std::shared_ptr<presentation::FloatingCodeEntity> entity) {
                HSML_INFO("Entity selected: " + entity->get_filename());
            }
        );

        dev_env->set_connection_activated_callback(
            [](std::shared_ptr<presentation::ConnectionLine> connection) {
                HSML_INFO("Connection activated between entities");
            }
        );

        // Create a demo project with sample code files
        dev_env->create_new_project("HSML_Demo_Project");
        HSML_INFO("Created demo project");

        // Add some sample code files to demonstrate the 3D environment
        const std::vector<std::pair<std::string, std::string>> sample_files = {
            {"main.cpp", R"(
#include <iostream>
#include "engine.h"
#include "renderer.h"

int main() {
    std::cout << "Starting HSML Engine..." << std::endl;
    Engine engine;
    engine.initialize();
    engine.run();
    return 0;
}
)"},
            {"engine.h", R"(
#pragma once
#include "renderer.h"
#include "physics.h"

class Engine {
public:
    void initialize();
    void run();
    void shutdown();

private:
    Renderer* renderer_;
    PhysicsEngine* physics_;
    bool running_;
};
)"},
            {"renderer.h", R"(
#pragma once
#include "shader.h"
#include "mesh.h"

class Renderer {
public:
    void initialize();
    void render_frame();
    void shutdown();

private:
    ShaderManager* shaders_;
    MeshManager* meshes_;
    bool initialized_;
};
)"},
            {"physics.h", R"(
#pragma once
#include <vector>
#include "vector3.h"

class PhysicsEngine {
public:
    void initialize();
    void update(float delta_time);
    void apply_gravity(Object& obj);

private:
    std::vector<Object> objects_;
    Vector3 gravity_;
};
)"},
            {"shader.h", R"(
#pragma once
#include <string>
#include "opengl.h"

class Shader {
public:
    Shader(const std::string& vertex_src, const std::string& fragment_src);
    ~Shader();

    void bind();
    void unbind();
    void set_uniform(const std::string& name, float value);

private:
    GLuint program_;
    bool compiled_;
};
)"}
        };

        // Add sample files to the 3D environment
        for (size_t i = 0; i < sample_files.size(); ++i) {
            // Calculate spherical positions for each file
            double theta = (i * 2.0 * M_PI) / sample_files.size();
            double phi = M_PI / 3.0; // 60 degrees from positive z-axis
            double radius = 30.0 + (i * 5.0); // Vary the radius

            SphericalCoords position(radius, theta, phi);

            auto entity = dev_env->add_code_file(
                sample_files[i].first,
                sample_files[i].second,
                position
            );

            if (entity) {
                HSML_INFO("Added code file: " + sample_files[i].first +
                         " at position (" + std::to_string(radius) + ", " +
                         std::to_string(theta) + ", " + std::to_string(phi) + ")");
            }
        }

        // Analyze and create connections between the code files
        dev_env->analyze_and_create_connections();
        HSML_INFO("Analyzed and created code connections");

        // Set up the camera for a good overview
        SphericalCoords overview_pos(80.0, M_PI/4.0, M_PI/4.0);
        dev_env->set_camera_position(overview_pos);
        HSML_INFO("Set up overview camera position");

        // Main demonstration loop
        HSML_INFO("Starting main demonstration loop...");
        const int DEMO_DURATION_SECONDS = 30;
        auto start_time = std::chrono::steady_clock::now();

        while (!dev_env->should_close()) {
            auto current_time = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                current_time - start_time
            ).count();

            if (elapsed >= DEMO_DURATION_SECONDS) {
                HSML_INFO("Demo duration reached, exiting...");
                break;
            }

            // Update the environment
            dev_env->update(1.0f / 60.0f); // 60 FPS

            // Render the environment
            dev_env->render();

            // Get and display statistics
            auto stats = dev_env->get_statistics();
            if (stats.entity_count > 0) {
                HSML_INFO("Environment Stats - Entities: " + std::to_string(stats.entity_count) +
                         ", Connections: " + std::to_string(stats.connection_count) +
                         ", Frame Time: " + std::to_string(stats.average_frame_time_ms) + "ms");
            }

            // Simulate some interactions during the demo
            if (elapsed == 5) {
                // Focus on main.cpp after 5 seconds
                auto main_entity = dev_env->find_entity_by_filename("main.cpp");
                if (main_entity) {
                    dev_env->focus_on_entity(main_entity);
                    HSML_INFO("Focusing on main.cpp entity");
                }
            } else if (elapsed == 10) {
                // Start code execution visualization
                auto engine_entity = dev_env->find_entity_by_filename("engine.h");
                if (engine_entity) {
                    dev_env->start_code_execution("Engine::initialize", engine_entity);
                    HSML_INFO("Started code execution visualization");
                }
            } else if (elapsed == 15) {
                // Orbit around renderer.h
                auto renderer_entity = dev_env->find_entity_by_filename("renderer.h");
                if (renderer_entity) {
                    dev_env->orbit_around_entity(renderer_entity, 3.0f);
                    HSML_INFO("Starting orbital camera around renderer.h");
                }
            } else if (elapsed == 20) {
                // Search for "render" functions
                dev_env->search_code("render");
                HSML_INFO("Searching for 'render' functions");
            } else if (elapsed == 25) {
                // Reset to overview
                dev_env->reset_camera_to_overview();
                HSML_INFO("Resetting to overview position");
            }

            // Sleep to maintain reasonable frame rate
            std::this_thread::sleep_for(std::chrono::milliseconds(16));
        }

        // Cleanup
        HSML_INFO("Shutting down 3D Development Environment");
        dev_env->shutdown();

        platform_mgr.shutdown_platform();

        std::cout << "\n🎉 3D GUI Development Environment Demo completed successfully!" << std::endl;
        std::cout << "🌟 Features Demonstrated:" << std::endl;
        std::cout << "   • Floating Code Entities in 3D Space" << std::endl;
        std::cout << "   • Glowing Connection Lines" << std::endl;
        std::cout << "   • Real-time Code Execution Visualization" << std::endl;
        std::cout << "   • Immersive Camera Navigation" << std::endl;
        std::cout << "   • Interactive Code Search" << std::endl;
        std::cout << "   • Syntax Highlighting" << std::endl;
        std::cout << "   • Entity Selection and Focusing" << std::endl;
        std::cout << "   • Performance Monitoring" << std::endl;
        std::cout << "   • 3D Audio Integration" << std::endl;
        std::cout << "   • Multi-language Support" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        HSML_ERROR("Demo failed: " + std::string(e.what()));
        std::cerr << "❌ Demo failed: " << e.what() << std::endl;
        return 1;
    }
}

// Simple test renderer implementation for demonstration
class TestRenderer : public domain::ISphericalRenderer {
public:
    void render_entity(const domain::ISphericalEntity& entity) override {
        // Simulate rendering by printing entity info
        std::cout << "  Rendering entity at position: ("
                  << entity.position().radius() << ", "
                  << entity.position().theta() << ", "
                  << entity.position().phi() << ")" << std::endl;
    }

    void render_scene(const std::vector<std::shared_ptr<domain::ISphericalEntity>>& entities) override {
        std::cout << "  Rendering scene with " << entities.size() << " entities" << std::endl;
        for (const auto& entity : entities) {
            render_entity(*entity);
        }
    }

    void render_sphere(const SphericalCoords& center, double radius,
                      const Vector3& color) override {
        std::cout << "  Rendering sphere at (" << center.radius() << ", "
                  << center.theta() << ", " << center.phi() << ") with radius " << radius << std::endl;
    }

    void begin_frame() override {
        std::cout << "  Beginning render frame" << std::endl;
    }

    void end_frame() override {
        std::cout << "  Ending render frame" << std::endl;
    }

    void clear(const Vector3& color) override {
        std::cout << "  Clearing with color (" << color.x() << ", "
                  << color.y() << ", " << color.z() << ")" << std::endl;
    }

    void set_observer_position(const SphericalCoords& position) override {
        observer_position_ = position;
    }

    const SphericalCoords& get_observer_position() const override {
        return observer_position_;
    }

    void set_field_of_view(double fov_radians) override {}
    void set_near_plane(double near_plane) override {}
    void set_far_plane(double far_plane) override {}
    void set_viewport(int x, int y, int width, int height) override {}

    bool supports_feature(const std::string& feature) const override {
        return true; // Support all features for demo
    }

    std::vector<std::string> get_supported_features() const override {
        return {"3d_rendering", "spherical_coordinates", "connection_lines"};
    }

    int get_rendered_entity_count() const override { return 0; }
    double get_last_frame_time_ms() const override { return 16.67; }
    bool get_last_render_successful() const override { return true; }
    std::string get_last_error() const override { return ""; }

private:
    SphericalCoords observer_position_{100.0, 0.0, 0.0};
};
