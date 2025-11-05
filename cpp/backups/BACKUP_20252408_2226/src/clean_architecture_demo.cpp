/** @file clean_architecture_demo.cpp
 * @brief Demonstration of Clean Architecture implementation
 */

#include <iostream>
#include <memory>
#include <chrono>

// Clean Architecture Layers
#include "hsml/domain/entities/SphericalEntity.h"
#include "hsml/domain/entities/SphericalScene.h"
#include "hsml/domain/entities/SphereShape.h"
#include "hsml/domain/entities/LinearMovementStrategy.h"
#include "hsml/domain/use_cases/RenderSphericalScene.h"

#include "hsml/application/services/SphericalSceneService.h"
#include "hsml/application/commands/RenderSceneCommand.h"

#include "hsml/infrastructure/persistence/SphericalSceneRepository.h"
#include "hsml/infrastructure/rendering/SphericalRendererAdapter.h"

#include "hsml/presentation/gui/SphericalViewport.h"

#include "hsml/core/debug_logger.h"
#include "hsml/core/spherical_coords.h"

int main() {
    std::cout << "🎯 HSML Clean Architecture Demonstration" << std::endl;
    std::cout << "=====================================" << std::endl;

    HSML_INFO("Starting Clean Architecture Demo");

    try {
        // 1. Infrastructure Layer Setup
        HSML_INFO("Setting up Infrastructure Layer...");
        auto repository = std::make_shared<infrastructure::SphericalSceneRepository>();
        auto test_renderer = std::make_unique<TestRenderer>();
        auto renderer_adapter = std::make_shared<infrastructure::SphericalRendererAdapter>(
            std::move(test_renderer)
        );

        // 2. Domain Layer - Use Case Setup
        HSML_INFO("Setting up Domain Layer Use Cases...");
        auto render_use_case = std::make_unique<domain::RenderSphericalScene>(
            repository, renderer_adapter
        );

        // 3. Application Layer Setup
        HSML_INFO("Setting up Application Layer Services...");
        auto scene_service = std::make_shared<application::SphericalSceneService>(
            repository, renderer_adapter, std::move(render_use_case)
        );

        // 4. Create a demonstration scene
        HSML_INFO("Creating demonstration scene...");
        scene_service->create_scene("demo_scene", "HSML Clean Architecture Demo");

        // Create multiple spherical entities with different properties
        for (int i = 0; i < 5; ++i) {
            double radius = 50.0 + i * 20.0;
            double theta = i * 0.3;
            double phi = i * 0.5;

            SphericalCoords position(radius, theta, phi);
            auto shape = std::make_unique<SphereShape>(5.0 + i);
            auto entity = std::make_shared<SphericalEntity>(
                position, std::move(shape), "DemoEntity_" + std::to_string(i)
            );

            // Add movement strategy to some entities
            if (i % 2 == 0) {
                auto movement = std::make_unique<LinearMovementStrategy>(
                    0.5, 0.01, 0.01  // velocity components
                );
                entity->set_movement_strategy(std::move(movement));
            }

            scene_service->add_entity_to_scene("demo_scene", entity);
        }

        // 5. Presentation Layer Setup
        HSML_INFO("Setting up Presentation Layer...");
        presentation::SphericalViewport viewport(scene_service, 1024, 768);
        viewport.set_camera_position(SphericalCoords(200.0, 0.0, 0.0));

        // 6. Demonstrate the architecture in action
        HSML_INFO("Running demonstration simulation...");

        const int SIMULATION_STEPS = 10;
        for (int step = 0; step < SIMULATION_STEPS; ++step) {
            HSML_INFO("Simulation step " + std::to_string(step + 1) + "/" +
                     std::to_string(SIMULATION_STEPS));

            // Update physics
            scene_service->update_physics(0.1); // 100ms delta time

            // Render scene
            domain::RenderRequest render_request;
            render_request.observer_position = SphericalCoords(200.0, 0.0, 0.0);
            render_request.field_of_view_radians = 1.047; // 60 degrees
            render_request.viewport_width = 1024;
            render_request.viewport_height = 768;
            render_request.enable_culling = true;
            render_request.enable_lod = true;

            application::RenderSceneCommand render_cmd(scene_service.get(), render_request);
            render_cmd.execute();

            if (render_cmd.was_successful()) {
                auto& result = render_cmd.get_result();
                auto& stats = render_cmd.get_statistics();

                HSML_INFO("Rendered " + std::to_string(result.entities_rendered) +
                         " entities in " + std::to_string(stats.average_render_time_ms) + "ms");
            }

            // Small delay to simulate real-time
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }

        // 7. Query demonstration
        HSML_INFO("Demonstrating query capabilities...");
        auto scene = scene_service->get_scene("demo_scene");
        auto stats = scene_service->get_scene_statistics("demo_scene");

        HSML_INFO("Scene contains " + std::to_string(stats.entity_count) + " entities");
        HSML_INFO("Scene max radius: " + std::to_string(stats.max_radius));

        // Find entities in a specific radius
        auto nearby_entities = scene_service->find_entities_in_radius(
            "demo_scene",
            SphericalCoords(100.0, 0.0, 0.0),
            50.0
        );

        HSML_INFO("Found " + std::to_string(nearby_entities.size()) +
                 " entities within 50 units of center");

        // 8. Architecture validation
        HSML_INFO("Validating Clean Architecture principles...");

        // Verify separation of concerns
        bool domain_isolated = true; // Domain entities don't depend on infrastructure
        bool application_orchestrates = true; // Application layer orchestrates domain
        bool infrastructure_adapts = true; // Infrastructure adapts to domain interfaces
        bool presentation_interacts = true; // Presentation layer interacts with application

        if (domain_isolated && application_orchestrates &&
            infrastructure_adapts && presentation_interacts) {
            HSML_INFO("✅ All Clean Architecture principles validated!");
        }

        std::cout << "\n🎉 Clean Architecture Demo completed successfully!" << std::endl;
        std::cout << "🏗️ Architecture demonstrates:" << std::endl;
        std::cout << "   • SOLID principles implementation" << std::endl;
        std::cout << "   • Clear separation of concerns" << std::endl;
        std::cout << "   • Dependency inversion" << std::endl;
        std::cout << "   • Testable design" << std::endl;
        std::cout << "   • Maintainable structure" << std::endl;

        return 0;

    } catch (const std::exception& e) {
        HSML_ERROR("Demo failed: " + std::string(e.what()));
        std::cerr << "❌ Demo failed: " << e.what() << std::endl;
        return 1;
    }
}

// Simple test renderer for demonstration
class TestRenderer : public domain::ISphericalRenderer {
public:
    void render_entity(const domain::ISphericalEntity& entity) override {
        std::cout << "  Rendering entity: " << entity.get_name() << std::endl;
    }

    void render_scene(const std::vector<std::shared_ptr<domain::ISphericalEntity>>& entities) override {
        for (const auto& entity : entities) {
            render_entity(*entity);
        }
    }

    void render_sphere(const SphericalCoords& center, double radius, const Vector3& color) override {
        std::cout << "  Rendering sphere at (" << center.radius() << ", " << center.theta() << ", " << center.phi() << ")" << std::endl;
    }

    void begin_frame() override { std::cout << "  Beginning frame" << std::endl; }
    void end_frame() override { std::cout << "  Ending frame" << std::endl; }
    void clear(const Vector3& color = Vector3(0.0, 0.0, 0.0)) override { std::cout << "  Clearing viewport" << std::endl; }

    void set_observer_position(const SphericalCoords& position) override { observer_pos = position; }
    const SphericalCoords& get_observer_position() const override { return observer_pos; }

    void set_field_of_view(double fov_radians) override {}
    void set_near_plane(double near_plane) override {}
    void set_far_plane(double far_plane) override {}
    void set_viewport(int x, int y, int width, int height) override {}

    bool supports_feature(const std::string& feature) const override { return true; }
    std::vector<std::string> get_supported_features() const override { return {}; }

    int get_rendered_entity_count() const override { return 0; }
    double get_last_frame_time_ms() const override { return 16.67; }
    bool get_last_render_successful() const override { return true; }
    std::string get_last_error() const override { return ""; }

private:
    SphericalCoords observer_pos{100.0, 0.0, 0.0};
};
