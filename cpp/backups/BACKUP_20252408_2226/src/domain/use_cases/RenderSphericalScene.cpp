/** @file RenderSphericalScene.h
 * @brief Use case for rendering spherical scenes
 *
 * Clean Architecture: Domain Layer Use Case
 * Implements the IRenderSphericalScene interface for scene rendering.
 */

#pragma once

#include "hsml/domain/interfaces/IRenderSphericalScene.h"
#include <memory>
#include <vector>
#include <string>

namespace hsml {
namespace domain {

// Forward declarations
class ISphericalEntity;
class SphericalScene;
class ISceneRepository;
class ISphericalRenderer;

/**
 * @brief Use case for rendering spherical scenes
 *
 * This class implements the business logic for rendering spherical scenes.
 * It coordinates between the scene repository and renderer while maintaining
 * domain rules and validation.
 */
class RenderSphericalScene : public IRenderSphericalScene {
public:
    /**
     * @brief Construct the render spherical scene use case
     *
     * @param scene_repo Repository for accessing scene data
     * @param renderer Spherical renderer for the scene
     */
    RenderSphericalScene(
        std::shared_ptr<ISceneRepository> scene_repo,
        std::shared_ptr<ISphericalRenderer> renderer
    );

    // IRenderSphericalScene implementation
    SphericalScene execute(const RenderRequest& request) override;
    bool validate_request(const RenderRequest& request) const override;
    RenderStatistics get_statistics() const override;

private:
    std::shared_ptr<ISceneRepository> scene_repository_;
    std::shared_ptr<ISphericalRenderer> renderer_;
    RenderStatistics statistics_;

    // Business logic methods
    std::vector<std::shared_ptr<ISphericalEntity>> filter_visible_entities(
        const SphericalScene& scene,
        const RenderRequest& request
    ) const;

    std::vector<std::shared_ptr<ISphericalEntity>> perform_frustum_culling(
        const std::vector<std::shared_ptr<ISphericalEntity>>& entities,
        const RenderRequest& request
    ) const;

    void apply_level_of_detail(
        std::vector<std::shared_ptr<ISphericalEntity>>& entities,
        const RenderRequest& request
    ) const;

    bool is_entity_visible(
        const ISphericalEntity& entity,
        const RenderRequest& request
    ) const;

    double calculate_lod_factor(
        const ISphericalEntity& entity,
        const SphericalCoords& observer
    ) const;

    // Validation helpers
    bool validate_observer_position(const SphericalCoords& position) const;
    bool validate_render_parameters(const RenderRequest& request) const;
    bool validate_viewport_dimensions(int width, int height) const;

    // Statistics tracking
    void update_statistics(const RenderResult& result);
    void reset_statistics();
};

} // namespace domain
} // namespace hsml
