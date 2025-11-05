/** @file SphericalRendererAdapter.h
 * @brief Infrastructure adapter for spherical rendering
 *
 * Clean Architecture: Infrastructure Layer
 * Adapts existing rendering infrastructure to domain interfaces.
 */

#pragma once

#include "hsml/domain/interfaces/ISphericalRenderer.h"
#include "hsml/rendering/renderer_interface.h"
#include <memory>
#include <unordered_map>

namespace hsml {
namespace infrastructure {

/**
 * @brief Adapter that bridges domain renderer interface with infrastructure
 *
 * Implements ISphericalRenderer using the existing RendererInterface
 * while maintaining clean architecture boundaries.
 */
class SphericalRendererAdapter : public domain::ISphericalRenderer {
public:
    /**
     * @brief Construct the renderer adapter
     *
     * @param renderer The underlying rendering infrastructure
     */
    explicit SphericalRendererAdapter(std::unique_ptr<rendering::RendererInterface> renderer);

    // ISphericalRenderer implementation
    void render_entity(const domain::ISphericalEntity& entity) override;
    void render_scene(const std::vector<std::shared_ptr<domain::ISphericalEntity>>& entities) override;

    void render_sphere(
        const SphericalCoords& center,
        double radius,
        const Vector3& color
    ) override;

    void begin_frame() override;
    void end_frame() override;
    void clear(const Vector3& color = Vector3(0.0, 0.0, 0.0)) override;

    void set_observer_position(const SphericalCoords& position) override;
    const SphericalCoords& get_observer_position() const override;

    void set_field_of_view(double fov_radians) override;
    void set_near_plane(double near_plane) override;
    void set_far_plane(double far_plane) override;

    void set_viewport(int x, int y, int width, int height) override;

    bool supports_feature(const std::string& feature) const override;
    std::vector<std::string> get_supported_features() const override;

    int get_rendered_entity_count() const override;
    double get_last_frame_time_ms() const override;
    bool get_last_render_successful() const override;

    std::string get_last_error() const override;

private:
    std::unique_ptr<rendering::RendererInterface> renderer_;
    SphericalCoords observer_position_;
    int rendered_entity_count_;
    double last_frame_time_ms_;
    bool last_render_successful_;
    std::string last_error_;

    // Entity tracking for statistics
    std::unordered_map<std::string, int> entity_render_counts_;

    // Helper methods
    void update_render_state();
    void convert_domain_to_rendering(const domain::ISphericalEntity& entity);
    rendering::RenderState create_render_state() const;
    bool validate_rendering_state() const;

    // Error handling
    void set_error(const std::string& error);
    void clear_error();

    // Performance monitoring
    void start_frame_timer();
    void end_frame_timer();
};

/**
 * @brief Factory for creating renderer adapters
 */
class RendererAdapterFactory {
public:
    static std::shared_ptr<domain::ISphericalRenderer> create_renderer_adapter(
        rendering::Backend backend = rendering::Backend::SOFTWARE
    );
};

/**
 * @brief Renderer registry for managing different backends
 */
class RendererRegistry {
public:
    static RendererRegistry& instance();

    void register_adapter_factory(
        rendering::Backend backend,
        std::function<std::shared_ptr<domain::ISphericalRenderer>()> factory
    );

    std::shared_ptr<domain::ISphericalRenderer> create_adapter(rendering::Backend backend);

    std::vector<rendering::Backend> get_available_backends() const;

private:
    RendererRegistry() = default;

    std::unordered_map<rendering::Backend,
                      std::function<std::shared_ptr<domain::ISphericalRenderer>()>> factories_;
};

} // namespace infrastructure
} // namespace hsml
