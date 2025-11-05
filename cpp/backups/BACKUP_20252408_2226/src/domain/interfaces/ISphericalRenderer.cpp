/** @file ISphericalRenderer.h
 * @brief Spherical renderer interface for domain layer
 *
 * Clean Architecture: Domain Layer Interface
 * Defines the contract for rendering spherical entities.
 */

#pragma once

#include <memory>
#include <vector>

namespace hsml {
namespace domain {

class ISphericalEntity;
class SphericalCoords;
class Vector3;

/**
 * @brief Interface for spherical rendering operations
 *
 * Defines the contract for rendering spherical entities in the domain.
 * This interface abstracts rendering details while providing domain-specific
 * rendering operations.
 */
class ISphericalRenderer {
public:
    virtual ~ISphericalRenderer() = default;

    // Rendering operations
    virtual void render_entity(const ISphericalEntity& entity) = 0;
    virtual void render_scene(const std::vector<std::shared_ptr<ISphericalEntity>>& entities) = 0;

    // Primitive rendering (for infrastructure layer to implement)
    virtual void render_sphere(
        const SphericalCoords& center,
        double radius,
        const Vector3& color
    ) = 0;

    // Frame management
    virtual void begin_frame() = 0;
    virtual void end_frame() = 0;
    virtual void clear(const Vector3& color = Vector3(0.0, 0.0, 0.0)) = 0;

    // State management
    virtual void set_observer_position(const SphericalCoords& position) = 0;
    virtual const SphericalCoords& get_observer_position() const = 0;

    // Rendering parameters
    virtual void set_field_of_view(double fov_radians) = 0;
    virtual void set_near_plane(double near_plane) = 0;
    virtual void set_far_plane(double far_plane) = 0;

    // Viewport management
    virtual void set_viewport(int x, int y, int width, int height) = 0;

    // Rendering capabilities
    virtual bool supports_feature(const std::string& feature) const = 0;
    virtual std::vector<std::string> get_supported_features() const = 0;

    // Statistics and debugging
    virtual int get_rendered_entity_count() const = 0;
    virtual double get_last_frame_time_ms() const = 0;
    virtual bool get_last_render_successful() const = 0;

    // Error handling
    virtual std::string get_last_error() const = 0;
};

/**
 * @brief Rendering feature flags
 */
struct RenderingFeatures {
    static const std::string LEVEL_OF_DETAIL;
    static const std::string FRUSTUM_CULLING;
    static const std::string ANTIALIASING;
    static const std::string SHADOWS;
    static const std::string REFLECTIONS;
    static const std::string TRANSPARENCY;
    static const std::string POST_PROCESSING;
};

/**
 * @brief Render result structure for domain operations
 */
struct DomainRenderResult {
    bool success{false};
    std::string error_message;
    int entities_rendered{0};
    double render_time_ms{0.0};
    int triangles_rendered{0};
};

} // namespace domain
} // namespace hsml
