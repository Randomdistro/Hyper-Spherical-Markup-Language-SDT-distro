/** @file IRenderSphericalScene.h
 * @brief Use case interface for rendering spherical scenes
 *
 * Clean Architecture: Domain Layer Use Case Interface
 * Defines the contract for rendering operations in the domain.
 */

#pragma once

#include <memory>
#include <vector>

namespace hsml {
namespace domain {

// Forward declarations
class ISphericalEntity;
class SphericalScene;
class RenderRequest;
class RenderResult;

/**
 * @brief Interface for rendering spherical scenes
 *
 * This use case defines the contract for rendering operations
 * in the spherical coordinate system. It encapsulates the
 * business logic for scene rendering without infrastructure details.
 */
class IRenderSphericalScene {
public:
    virtual ~IRenderSphericalScene() = default;

    /**
     * @brief Execute the render spherical scene use case
     *
     * @param request The render request containing scene parameters
     * @return RenderResult containing the rendered scene data
     */
    virtual SphericalScene execute(const RenderRequest& request) = 0;

    /**
     * @brief Validate if a render request is valid
     *
     * @param request The render request to validate
     * @return true if valid, false otherwise
     */
    virtual bool validate_request(const RenderRequest& request) const = 0;

    /**
     * @brief Get rendering statistics for the last operation
     *
     * @return Statistics about the rendering process
     */
    virtual RenderStatistics get_statistics() const = 0;
};

/**
 * @brief Render request structure
 */
struct RenderRequest {
    SphericalCoords observer_position;
    double field_of_view_radians{1.047};  // 60 degrees default
    double near_plane{0.1};
    double far_plane{10000.0};
    int viewport_width{800};
    int viewport_height{600};
    bool enable_culling{true};
    bool enable_lod{true};
    double max_render_distance{1000.0};
    std::vector<std::string> entity_filter;  // Empty = render all
};

/**
 * @brief Render result structure
 */
struct RenderResult {
    bool success{false};
    std::string error_message;
    std::vector<std::shared_ptr<ISphericalEntity>> visible_entities;
    std::vector<std::shared_ptr<ISphericalEntity>> culled_entities;
    double render_time_ms{0.0};
    int triangles_rendered{0};
    int entities_rendered{0};
    int lod_levels_used{0};
};

/**
 * @brief Rendering statistics
 */
struct RenderStatistics {
    int total_renders{0};
    double average_render_time_ms{0.0};
    double min_render_time_ms{0.0};
    double max_render_time_ms{0.0};
    int total_entities_rendered{0};
    int total_triangles_rendered{0};
    int total_culled_entities{0};
    double average_cull_ratio{0.0};
};

} // namespace domain
} // namespace hsml
