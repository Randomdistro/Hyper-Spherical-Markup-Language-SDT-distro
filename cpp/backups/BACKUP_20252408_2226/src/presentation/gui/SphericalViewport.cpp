/** @file SphericalViewport.h
 * @brief Presentation layer GUI component for spherical viewports
 *
 * Clean Architecture: Presentation Layer
 * Provides GUI component for displaying spherical scenes.
 */

#pragma once

#include "hsml/application/services/SphericalSceneService.h"
#include "hsml/domain/entities/SphericalCoords.h"
#include <memory>
#include <functional>
#include <string>

namespace hsml {
namespace presentation {

// Forward declarations
class SphericalViewportController;

/**
 * @brief GUI component for displaying spherical scenes
 *
 * This is a presentation layer component that handles user interaction
 * and delegates business logic to the application layer.
 */
class SphericalViewport {
public:
    /**
     * @brief Construct a spherical viewport
     *
     * @param service The application service for scene operations
     * @param width Initial viewport width
     * @param height Initial viewport height
     */
    SphericalViewport(
        std::shared_ptr<application::SphericalSceneService> service,
        int width = 800,
        int height = 600
    );

    ~SphericalViewport();

    // Viewport management
    void set_size(int width, int height);
    void get_size(int& width, int& height) const;
    double get_aspect_ratio() const;

    // Camera control
    void set_camera_position(const domain::SphericalCoords& position);
    const domain::SphericalCoords& get_camera_position() const;

    void set_camera_orientation(double heading, double pitch);
    void get_camera_orientation(double& heading, double& pitch) const;

    // Rendering
    void render_frame();
    void clear_viewport();

    // Event handling
    void handle_mouse_press(int x, int y, int button);
    void handle_mouse_release(int x, int y, int button);
    void handle_mouse_move(int x, int y);
    void handle_mouse_wheel(int delta);
    void handle_key_press(int key);
    void handle_key_release(int key);

    // Entity interaction
    void select_entity_at_position(int x, int y);
    std::string get_selected_entity_id() const;

    // Viewport state
    bool is_initialized() const;
    bool initialize();
    void shutdown();

    // Event callbacks
    using EntitySelectedCallback = std::function<void(const std::string& entity_id)>;
    using CameraChangedCallback = std::function<void(const domain::SphericalCoords& position)>;

    void set_entity_selected_callback(EntitySelectedCallback callback);
    void set_camera_changed_callback(CameraChangedCallback callback);

private:
    std::shared_ptr<application::SphericalSceneService> service_;
    std::unique_ptr<SphericalViewportController> controller_;

    // Viewport state
    int width_;
    int height_;
    bool initialized_;

    // Camera state
    domain::SphericalCoords camera_position_;
    double camera_heading_;   // Left/right rotation
    double camera_pitch_;     // Up/down rotation

    // Interaction state
    bool mouse_pressed_;
    int last_mouse_x_;
    int last_mouse_y_;
    std::string selected_entity_id_;

    // Callbacks
    EntitySelectedCallback entity_selected_callback_;
    CameraChangedCallback camera_changed_callback_;

    // Helper methods
    void update_camera_from_mouse(int delta_x, int delta_y);
    void zoom_camera(double factor);
    domain::SphericalCoords screen_to_spherical(int x, int y) const;
    std::pair<int, int> spherical_to_screen(const domain::SphericalCoords& coords) const;

    void notify_camera_changed();
    void notify_entity_selected(const std::string& entity_id);
};

/**
 * @brief Controller for spherical viewport behavior
 */
class SphericalViewportController {
public:
    explicit SphericalViewportController(SphericalViewport* viewport);

    void handle_camera_movement(const domain::SphericalCoords& new_position);
    void handle_entity_selection(const std::string& entity_id);
    void handle_render_request();

    // Control modes
    enum class ControlMode {
        ORBIT,
        FPS,
        FREE
    };

    void set_control_mode(ControlMode mode);
    ControlMode get_control_mode() const;

private:
    SphericalViewport* viewport_;
    ControlMode control_mode_;

    // Movement parameters
    double movement_speed_;
    double rotation_speed_;
    double zoom_speed_;
};

} // namespace presentation
} // namespace hsml
