/** @file SphericalSceneService.h
 * @brief Application service for spherical scene operations
 *
 * Clean Architecture: Application Layer Service
 * Orchestrates domain operations for spherical scene management.
 */

#pragma once

#include "hsml/domain/interfaces/ISceneRepository.h"
#include "hsml/domain/interfaces/ISphericalRenderer.h"
#include "hsml/domain/interfaces/IRenderSphericalScene.h"
#include <memory>
#include <string>
#include <vector>

namespace hsml {
namespace application {

class SphericalSceneService {
public:
    /**
     * @brief Construct the spherical scene service
     *
     * @param scene_repo Repository for scene data access
     * @param renderer Spherical renderer for scene operations
     * @param render_use_case Use case for rendering operations
     */
    SphericalSceneService(
        std::shared_ptr<domain::ISceneRepository> scene_repo,
        std::shared_ptr<domain::ISphericalRenderer> renderer,
        std::shared_ptr<domain::IRenderSphericalScene> render_use_case
    );

    // Scene management
    domain::SphericalScene get_scene(const std::string& scene_id);
    void save_scene(const domain::SphericalScene& scene);
    void create_scene(const std::string& scene_id, const std::string& name = "");
    void delete_scene(const std::string& scene_id);

    // Entity operations
    std::vector<std::shared_ptr<domain::ISphericalEntity>> get_visible_entities(
        const domain::SphericalCoords& observer_position
    );

    void update_entity_position(
        const std::string& entity_id,
        const domain::SphericalCoords& position
    );

    void add_entity_to_scene(
        const std::string& scene_id,
        std::shared_ptr<domain::ISphericalEntity> entity
    );

    void remove_entity_from_scene(
        const std::string& scene_id,
        const std::string& entity_id
    );

    // Rendering operations
    domain::SphericalScene render_scene(const domain::RenderRequest& request);

    // Query operations
    std::vector<std::shared_ptr<domain::ISphericalEntity>> find_entities_in_radius(
        const std::string& scene_id,
        const domain::SphericalCoords& center,
        double radius
    );

    std::vector<std::shared_ptr<domain::ISphericalEntity>> find_entities_by_type(
        const std::string& scene_id,
        const std::string& type
    );

    // Physics operations
    void update_physics(double delta_time);
    void detect_collisions(const std::string& scene_id);

    // Scene statistics
    struct SceneStatistics {
        size_t entity_count{0};
        double max_radius{0.0};
        domain::SphericalCoords centroid;
        std::vector<std::string> entity_types;
    };

    SceneStatistics get_scene_statistics(const std::string& scene_id);

private:
    std::shared_ptr<domain::ISceneRepository> scene_repository_;
    std::shared_ptr<domain::ISphericalRenderer> renderer_;
    std::shared_ptr<domain::IRenderSphericalScene> render_use_case_;

    // Cached current scene
    std::string current_scene_id_;
    domain::SphericalScene current_scene_;

    // Helper methods
    void ensure_scene_loaded(const std::string& scene_id);
    void validate_entity_id(const std::string& entity_id);
    void validate_scene_id(const std::string& scene_id);

    // Physics helpers
    void resolve_collisions(
        std::vector<std::shared_ptr<domain::ISphericalEntity>>& entities
    );

    void apply_gravity(
        std::vector<std::shared_ptr<domain::ISphericalEntity>>& entities,
        double delta_time
    );
};

} // namespace application
} // namespace hsml
