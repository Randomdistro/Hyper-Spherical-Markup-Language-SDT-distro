/** @file test_clean_architecture.cpp
 * @brief Comprehensive test for Clean Architecture implementation
 */

#include <iostream>
#include <memory>
#include <cassert>

// Domain Layer
#include "hsml/domain/entities/SphericalEntity.h"
#include "hsml/domain/entities/SphericalScene.h"
#include "hsml/domain/entities/SphereShape.h"
#include "hsml/domain/entities/LinearMovementStrategy.h"
#include "hsml/domain/use_cases/RenderSphericalScene.h"

// Application Layer
#include "hsml/application/services/SphericalSceneService.h"
#include "hsml/application/commands/RenderSceneCommand.h"
#include "hsml/application/queries/SceneQueries.h"

// Infrastructure Layer
#include "hsml/infrastructure/persistence/SphericalSceneRepository.h"
#include "hsml/infrastructure/rendering/SphericalRendererAdapter.h"

// Presentation Layer
#include "hsml/presentation/gui/SphericalViewport.h"

// Core
#include "hsml/core/debug_logger.h"
#include "hsml/core/spherical_coords.h"

// Test utilities
class TestRenderer : public domain::ISphericalRenderer {
public:
    void render_entity(const domain::ISphericalEntity& entity) override {
        rendered_entities.push_back(entity.get_id());
    }

    void render_scene(const std::vector<std::shared_ptr<domain::ISphericalEntity>>& entities) override {
        for (const auto& entity : entities) {
            render_entity(*entity);
        }
    }

    void render_sphere(const SphericalCoords& center, double radius, const Vector3& color) override {
        render_calls++;
    }

    void begin_frame() override { frames_started++; }
    void end_frame() override { frames_ended++; }
    void clear(const Vector3& color = Vector3(0.0, 0.0, 0.0)) override { clears++; }

    void set_observer_position(const SphericalCoords& position) override {
        observer_position = position;
    }

    const SphericalCoords& get_observer_position() const override {
        return observer_position;
    }

    void set_field_of_view(double fov_radians) override { fov = fov_radians; }
    void set_near_plane(double near_plane) override { near = near_plane; }
    void set_far_plane(double far_plane) override { far = far_plane; }
    void set_viewport(int x, int y, int width, int height) override {
        viewport_x = x; viewport_y = y;
        viewport_width = width; viewport_height = height;
    }

    bool supports_feature(const std::string& feature) const override { return true; }
    std::vector<std::string> get_supported_features() const override { return {}; }

    int get_rendered_entity_count() const override { return rendered_entities.size(); }
    double get_last_frame_time_ms() const override { return 16.67; }
    bool get_last_render_successful() const override { return true; }
    std::string get_last_error() const override { return ""; }

    // Test inspection
    std::vector<std::string> rendered_entities;
    int render_calls = 0;
    int frames_started = 0;
    int frames_ended = 0;
    int clears = 0;
    SphericalCoords observer_position{100.0, 0.0, 0.0};
    double fov = 1.0;
    double near = 0.1;
    double far = 1000.0;
    int viewport_x = 0, viewport_y = 0, viewport_width = 800, viewport_height = 600;
};

void test_domain_layer() {
    std::cout << "🧪 Testing Domain Layer..." << std::endl;

    // Test SphericalEntity
    SphericalCoords position(100.0, 0.5, 1.0);
    auto shape = std::make_unique<SphereShape>(10.0);
    SphericalEntity entity(position, std::move(shape), "TestEntity");

    assert(entity.get_name() == "TestEntity");
    assert(entity.position().radius() == 100.0);
    assert(entity.shape().get_type() == "sphere");
    assert(entity.is_active() == true);

    // Test movement strategy
    auto movement_strategy = std::make_unique<LinearMovementStrategy>(1.0, 0.1, 0.1);
    entity.set_movement_strategy(std::move(movement_strategy));

    entity.move_by(10.0, 0.1, 0.1);
    assert(entity.position().radius() == 110.0);

    // Test SphericalScene
    SphericalScene scene;
    scene.add_entity(std::make_shared<SphericalEntity>(entity));
    assert(scene.get_entity_count() == 1);

    std::cout << "✅ Domain Layer tests passed!" << std::endl;
}

void test_application_layer() {
    std::cout << "🧪 Testing Application Layer..." << std::endl;

    // Create test infrastructure
    auto repository = std::make_shared<infrastructure::SphericalSceneRepository>();
    auto renderer = std::make_shared<TestRenderer>();
    auto render_use_case = std::make_unique<domain::RenderSphericalScene>(repository, renderer);

    // Test SphericalSceneService
    application::SphericalSceneService service(repository, renderer, std::move(render_use_case));

    // Create test scene
    service.create_scene("test_scene", "Test Scene");

    // Create test entity
    SphericalCoords position(100.0, 0.0, 0.0);
    auto shape = std::make_unique<SphereShape>(10.0);
    auto entity = std::make_shared<SphericalEntity>(position, std::move(shape), "TestEntity");

    service.add_entity_to_scene("test_scene", entity);

    // Test scene retrieval
    auto scene = service.get_scene("test_scene");
    assert(scene.get_entity_count() == 1);

    // Test render command
    domain::RenderRequest request;
    request.observer_position = SphericalCoords(200.0, 0.0, 0.0);
    request.viewport_width = 800;
    request.viewport_height = 600;

    application::RenderSceneCommand render_cmd(&service, request);
    assert(render_cmd.can_execute());

    render_cmd.execute();
    assert(render_cmd.was_successful());

    std::cout << "✅ Application Layer tests passed!" << std::endl;
}

void test_infrastructure_layer() {
    std::cout << "🧪 Testing Infrastructure Layer..." << std::endl;

    // Test repository
    auto repository = std::make_shared<infrastructure::SphericalSceneRepository>();

    // Create and save scene
    domain::SphericalScene scene;
    repository->save(scene);

    // Retrieve scene
    auto retrieved_scene = repository->get_by_id("default_scene");
    assert(retrieved_scene.get_entity_count() == 0);

    // Test renderer adapter
    auto test_renderer = std::make_unique<TestRenderer>();
    auto adapter = std::make_shared<infrastructure::SphericalRendererAdapter>(std::move(test_renderer));

    SphericalCoords observer_pos(100.0, 0.0, 0.0);
    adapter->set_observer_position(observer_pos);
    assert(adapter->get_observer_position().radius() == 100.0);

    adapter->begin_frame();
    adapter->clear();
    adapter->end_frame();

    std::cout << "✅ Infrastructure Layer tests passed!" << std::endl;
}

void test_presentation_layer() {
    std::cout << "🧪 Testing Presentation Layer..." << std::endl;

    // Create application service
    auto repository = std::make_shared<infrastructure::SphericalSceneRepository>();
    auto renderer = std::make_shared<TestRenderer>();
    auto render_use_case = std::make_unique<domain::RenderSphericalScene>(repository, renderer);
    auto service = std::make_shared<application::SphericalSceneService>(
        repository, renderer, std::move(render_use_case)
    );

    // Test SphericalViewport
    presentation::SphericalViewport viewport(service, 800, 600);

    // Test viewport operations
    viewport.set_camera_position(SphericalCoords(100.0, 0.0, 0.0));
    int width, height;
    viewport.get_size(width, height);
    assert(width == 800 && height == 600);

    std::cout << "✅ Presentation Layer tests passed!" << std::endl;
}

void test_clean_architecture_integration() {
    std::cout << "🧪 Testing Clean Architecture Integration..." << std::endl;

    // Create all layers
    auto repository = std::make_shared<infrastructure::SphericalSceneRepository>();
    auto renderer = std::make_shared<TestRenderer>();
    auto render_use_case = std::make_unique<domain::RenderSphericalScene>(repository, renderer);
    auto service = std::make_shared<application::SphericalSceneService>(
        repository, renderer, std::move(render_use_case)
    );

    // Create presentation layer
    presentation::SphericalViewport viewport(service, 800, 600);

    // Test end-to-end flow: Create entity -> Add to scene -> Render
    service->create_scene("integration_test", "Integration Test Scene");

    SphericalCoords position(50.0, 0.0, 0.0);
    auto shape = std::make_unique<SphereShape>(5.0);
    auto entity = std::make_shared<SphericalEntity>(position, std::move(shape), "IntegrationEntity");

    service->add_entity_to_scene("integration_test", entity);

    // Render the scene
    domain::RenderRequest request;
    request.observer_position = SphericalCoords(100.0, 0.0, 0.0);
    request.viewport_width = 800;
    request.viewport_height = 600;

    auto rendered_scene = service->render_scene(request);
    assert(rendered_scene.get_entity_count() == 1);

    // Verify entity was rendered
    auto test_renderer = dynamic_cast<TestRenderer*>(renderer.get());
    assert(test_renderer->get_rendered_entity_count() == 1);

    std::cout << "✅ Clean Architecture Integration tests passed!" << std::endl;
}

int main() {
    std::cout << "🎯 Clean Architecture Comprehensive Test Suite" << std::endl;
    std::cout << "===============================================" << std::endl;

    try {
        test_domain_layer();
        test_application_layer();
        test_infrastructure_layer();
        test_presentation_layer();
        test_clean_architecture_integration();

        std::cout << "🎉 All Clean Architecture tests passed!" << std::endl;
        std::cout << "🏗️ Architecture is SOLID and well-structured!" << std::endl;
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "❌ Test failed: " << e.what() << std::endl;
        return 1;
    }
}
