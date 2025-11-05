/** @file port3r_studio_system.h
 * @brief P0RT3R Studio System Header
 *
 * Header declarations for the P0RT3R Studio System implementation.
 */

#pragma once

#include "../core/math/spherical_coords.h"
#include "../core/math/vector3.h"
#include <string>
#include <vector>
#include <memory>
#include <functional>

namespace hsml {
namespace core {
    // Forward declarations for core types
    struct SphericalCoords;
    struct Vector3;
}

namespace port3r {

// Forward declarations
class StudioElement;
class StudioPanel;
class StudioWindow;
class Port3rStudioSystem;

// ============================================================================
// StudioElement - Base UI Element
// ============================================================================

class StudioElement {
public:
    StudioElement(const std::string& id, const core::SphericalCoords& position);
    virtual ~StudioElement() = default;

    // Properties
    const std::string& getId() const { return id_; }
    const core::SphericalCoords& getPosition() const { return position_; }
    const core::Vector3& getColor() const { return color_; }
    bool isVisible() const { return visible_; }
    bool isActive() const { return active_; }

    // Setters
    void setPosition(const core::SphericalCoords& position);
    void setColor(const core::Vector3& color);
    void setVisible(bool visible);
    void setActive(bool active);

    // Animation
    void moveTo(const core::SphericalCoords& target, double duration);
    void update(double delta_time);

    // Interaction
    virtual bool handleClick(const core::SphericalCoords& click_point);
    virtual bool handleHover(const core::SphericalCoords& hover_point);

    // Rendering
    virtual void render();

protected:
    std::string id_;
    core::SphericalCoords position_;
    core::Vector3 color_;
    bool visible_ = true;
    bool active_ = true;

    // Animation state
    core::SphericalCoords target_position_;
    double animation_duration_ = 0.0;
    double animation_timer_ = 0.0;
    bool animating_ = false;

    // Interaction
    double interaction_radius_ = 0.1; // radians
};

// ============================================================================
// StudioPanel - Container Element
// ============================================================================

class StudioPanel : public StudioElement {
public:
    StudioPanel(const std::string& id, const core::SphericalCoords& position,
                double width, double height);

    void addChild(std::shared_ptr<StudioElement> child);
    void removeChild(const std::string& child_id);
    std::shared_ptr<StudioElement> getChild(const std::string& child_id) const;

    // Panel-specific properties
    void setBackgroundColor(const core::Vector3& color);
    void setBorderColor(const core::Vector3& color);
    void setBorderWidth(double width);

    // Overrides
    bool handleClick(const core::SphericalCoords& click_point) override;
    bool handleHover(const core::SphericalCoords& hover_point) override;
    void render() override;

private:
    double width_;
    double height_;
    core::Vector3 background_color_;
    core::Vector3 border_color_;
    double border_width_;

    std::vector<std::shared_ptr<StudioElement>> children_;

    bool pointInPanel(const core::SphericalCoords& point) const;
};

// ============================================================================
// StudioWindow - Window Element
// ============================================================================

class StudioWindow : public StudioPanel {
public:
    StudioWindow(const std::string& id, const core::SphericalCoords& position,
                 double width, double height, const std::string& title);

    const std::string& getTitle() const { return title_; }
    void setTitle(const std::string& title) { title_ = title; }

    void setDraggable(bool draggable) { draggable_ = draggable; }
    void setResizable(bool resizable) { resizable_ = resizable; }
    void setClosable(bool closable) { closable_ = closable; }

    // Window state
    void minimize();
    void maximize();
    void close();
    bool isMinimized() const { return minimized_; }

    // Overrides
    bool handleClick(const core::SphericalCoords& click_point) override;
    void render() override;

private:
    std::string title_;
    bool draggable_ = true;
    bool resizable_ = true;
    bool closable_ = true;
    bool minimized_ = false;

    // Window chrome
    double title_bar_height_ = 0.05;
    core::Vector3 title_bar_color_;

    bool pointInTitleBar(const core::SphericalCoords& point) const;
};

// ============================================================================
// Port3rStudioSystem - Main System
// ============================================================================

class Port3rStudioSystem {
public:
    static Port3rStudioSystem& instance();

    bool initialize();
    void shutdown();
    void update(double delta_time);
    void render();

    // Window management
    std::shared_ptr<StudioWindow> createWindow(const std::string& id,
                                              const core::SphericalCoords& position,
                                              double width, double height,
                                              const std::string& title);
    void destroyWindow(const std::string& id);
    std::shared_ptr<StudioWindow> getWindow(const std::string& id) const;

    // Panel management
    std::shared_ptr<StudioPanel> createPanel(const std::string& id,
                                            const core::SphericalCoords& position,
                                            double width, double height);
    void destroyPanel(const std::string& id);
    std::shared_ptr<StudioPanel> getPanel(const std::string& id) const;

    // Input handling
    void handleMouseClick(const core::SphericalCoords& position);
    void handleMouseHover(const core::SphericalCoords& position);
    void handleMouseMove(const core::SphericalCoords& position);

    // Camera/viewport
    void setCameraPosition(const core::SphericalCoords& position);
    const core::SphericalCoords& getCameraPosition() const { return camera_position_; }

private:
    Port3rStudioSystem() = default;
    ~Port3rStudioSystem() = default;

    // Prevent copying
    Port3rStudioSystem(const Port3rStudioSystem&) = delete;
    Port3rStudioSystem& operator=(const Port3rStudioSystem&) = delete;

    bool initialized_ = false;
    core::SphericalCoords camera_position_;

    // Element collections
    std::vector<std::shared_ptr<StudioWindow>> windows_;
    std::vector<std::shared_ptr<StudioPanel>> panels_;
    std::vector<std::shared_ptr<StudioElement>> elements_;

    // Utility methods
    std::shared_ptr<StudioElement> findElementAt(const core::SphericalCoords& position) const;
    void sortElementsByDepth();
    void updateAnimations(double delta_time);
};

} // namespace port3r
} // namespace hsml
