/**
 * Spatial Interface System (SIS) - 3D-Native GUI Architecture
 * NO CARTESIAN COORDINATES - PURE SPHERICAL REALITY
 * 
 * Revolutionary concept: All GUI controls are 3D spatial objects within the virtual world
 * - Controls can be moved, rotated, scaled in 3D space
 * - Collapsible/expandable 3D control panels
 - Spatial positioning for optimal workflow
 * - No traditional 2D overlays - everything is a 3D component
 */

#pragma once

#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include "hsml/rendering/software_renderer.h"
#include "hsml/core/state_tensor.h"
#include "hsml/rendering/concurrent_spherical_renderer.h"

#include <memory>
#include <vector>
#include <map>
#include <functional>
#include <string>
#include <variant>

namespace hsml::gui {

// Forward declarations
class SpatialControl;
class SpatialPanel;
class SpatialButton;
class SpatialSlider;
class SpatialToggle;
class SpatialDropdown;
class SpatialTextInput;
class SpatialColorPicker;
class SpatialViewport;

/**
 * Spatial Interface Element - Base class for all 3D GUI components
 * Every control is a spatial object with position, orientation, and behavior
 */
class SpatialElement {
public:
    struct SpatialState {
        core::SphericalCoords position{100.0, 0.0, 0.0};
        core::SphericalCoords orientation{0.0, 0.0, 0.0};
        double scale{1.0};
        bool visible{true};
        bool interactive{true};
        bool collapsed{false};
        core::Vector3 color{1.0, 1.0, 1.0};
        double opacity{1.0};
    };

    explicit SpatialElement(const std::string& id, const SpatialState& state = {});
    virtual ~SpatialElement() = default;

    // Spatial positioning and manipulation
    void setPosition(const core::SphericalCoords& pos);
    void setOrientation(const core::SphericalCoords& orient);
    void setScale(double scale);
    void moveTo(const core::SphericalCoords& target, double duration = 0.5);
    void rotateTo(const core::SphericalCoords& target, double duration = 0.5);
    
    // Visibility and interaction
    void show();
    void hide();
    void collapse();
    void expand();
    void setInteractive(bool interactive);
    
    // Spatial relationships
    void attachTo(std::shared_ptr<SpatialElement> parent);
    void detach();
    std::vector<std::shared_ptr<SpatialElement>> getChildren() const;
    
    // Rendering
    virtual void render(rendering::SoftwareRenderer& renderer) = 0;
    virtual void update(double delta_time) = 0;
    
    // Interaction
    virtual bool handleSpatialInteraction(const core::SphericalCoords& interaction_point);
    virtual void onSpatialHover();
    virtual void onSpatialClick();
    virtual void onSpatialDrag(const core::SphericalCoords& delta);

    // Getters
    const std::string& getId() const { return id_; }
    const SpatialState& getState() const { return state_; }
    bool isVisible() const { return state_.visible; }
    bool isInteractive() const { return state_.interactive; }

protected:
    std::string id_;
    SpatialState state_;
    std::shared_ptr<SpatialElement> parent_;
    std::vector<std::shared_ptr<SpatialElement>> children_;
    
    // Animation state
    struct Animation {
        core::SphericalCoords start_pos;
        core::SphericalCoords end_pos;
        double duration;
        double elapsed{0.0};
        bool active{false};
    } animation_;
};

/**
 * Spatial Button - 3D button that exists in spherical space
 */
class SpatialButton : public SpatialElement {
public:
    explicit SpatialButton(const std::string& id, const std::string& label = "");
    
    void setLabel(const std::string& label);
    void setAction(std::function<void()> action);
    void setButtonStyle(const std::string& style); // "flat", "raised", "glowing"
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    void onSpatialClick() override;

private:
    std::string label_;
    std::function<void()> action_;
    std::string style_{"flat"};
    bool pressed_{false};
    double press_animation_{0.0};
};

/**
 * Spatial Slider - 3D slider that can be dragged in spherical space
 */
class SpatialSlider : public SpatialElement {
public:
    explicit SpatialSlider(const std::string& id, double min_val = 0.0, double max_val = 1.0);
    
    void setValue(double value);
    void setRange(double min_val, double max_val);
    void setOnValueChanged(std::function<void(double)> callback);
    void setSliderStyle(const std::string& style); // "linear", "radial", "spherical"
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    void onSpatialDrag(const core::SphericalCoords& delta) override;

private:
    double min_value_{0.0};
    double max_value_{1.0};
    double current_value_{0.5};
    std::function<void(double)> on_value_changed_;
    std::string style_{"linear"};
    bool dragging_{false};
    core::SphericalCoords drag_start_;
};

/**
 * Spatial Toggle - 3D toggle switch in spherical space
 */
class SpatialToggle : public SpatialElement {
public:
    explicit SpatialToggle(const std::string& id, bool initial_state = false);
    
    void setState(bool state);
    void setOnStateChanged(std::function<void(bool)> callback);
    void setToggleStyle(const std::string& style); // "switch", "checkbox", "radio"
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    void onSpatialClick() override;

private:
    bool state_{false};
    std::function<void(bool)> on_state_changed_;
    std::string style_{"switch"};
    double animation_progress_{0.0};
};

/**
 * Spatial Dropdown - 3D dropdown menu in spherical space
 */
class SpatialDropdown : public SpatialElement {
public:
    explicit SpatialDropdown(const std::string& id);
    
    void addOption(const std::string& option);
    void setSelectedOption(const std::string& option);
    void setOnSelectionChanged(std::function<void(const std::string&)> callback);
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    void onSpatialClick() override;

private:
    std::vector<std::string> options_;
    std::string selected_option_;
    std::function<void(const std::string&)> on_selection_changed_;
    bool expanded_{false};
    double expand_animation_{0.0};
};

/**
 * Spatial Text Input - 3D text input field in spherical space
 */
class SpatialTextInput : public SpatialElement {
public:
    explicit SpatialTextInput(const std::string& id, const std::string& placeholder = "");
    
    void setText(const std::string& text);
    void setPlaceholder(const std::string& placeholder);
    void setOnTextChanged(std::function<void(const std::string&)> callback);
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    void onSpatialClick() override;

private:
    std::string text_;
    std::string placeholder_;
    std::function<void(const std::string&)> on_text_changed_;
    bool focused_{false};
    double cursor_blink_{0.0};
};

/**
 * Spatial Color Picker - 3D color picker in spherical space
 */
class SpatialColorPicker : public SpatialElement {
public:
    explicit SpatialColorPicker(const std::string& id, const core::Vector3& initial_color = {1.0, 1.0, 1.0});
    
    void setColor(const core::Vector3& color);
    void setOnColorChanged(std::function<void(const core::Vector3&)> callback);
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;
    void onSpatialDrag(const core::SphericalCoords& delta) override;

private:
    core::Vector3 color_{1.0, 1.0, 1.0};
    std::function<void(const core::Vector3&)> on_color_changed_;
    bool color_wheel_expanded_{false};
    double hue_{0.0};
    double saturation_{1.0};
    double value_{1.0};
};

/**
 * Spatial Panel - 3D panel that can contain other spatial controls
 */
class SpatialPanel : public SpatialElement {
public:
    explicit SpatialPanel(const std::string& id, const std::string& title = "");
    
    void setTitle(const std::string& title);
    void addControl(std::shared_ptr<SpatialElement> control);
    void removeControl(const std::string& control_id);
    void setPanelStyle(const std::string& style); // "floating", "docked", "collapsible"
    void setLayout(const std::string& layout); // "vertical", "horizontal", "grid", "radial"
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;

private:
    std::string title_;
    std::string style_{"floating"};
    std::string layout_{"vertical"};
    std::vector<std::shared_ptr<SpatialElement>> controls_;
    bool title_bar_visible_{true};
    double panel_animation_{0.0};
};

/**
 * Spatial Viewport - 3D viewport for rendering content
 */
class SpatialViewport : public SpatialElement {
public:
    explicit SpatialViewport(const std::string& id, size_t width = 800, size_t height = 600);
    
    void setSize(size_t width, size_t height);
    void setContent(std::shared_ptr<SpatialElement> content);
    void setViewportStyle(const std::string& style); // "window", "fullscreen", "floating"
    
    void render(rendering::SoftwareRenderer& renderer) override;
    void update(double delta_time) override;
    bool handleSpatialInteraction(const core::SphericalCoords& interaction_point) override;

private:
    size_t width_{800};
    size_t height_{600};
    std::shared_ptr<SpatialElement> content_;
    std::string style_{"window"};
    rendering::FrameBuffer viewport_buffer_;
};

/**
 * Spatial Interface Manager - Manages all 3D GUI elements
 */
class SpatialInterfaceManager {
public:
    static SpatialInterfaceManager& instance();
    
    // Control management
    std::shared_ptr<SpatialElement> createControl(const std::string& type, const std::string& id);
    void destroyControl(const std::string& id);
    std::shared_ptr<SpatialElement> getControl(const std::string& id);
    
    // Panel management
    std::shared_ptr<SpatialPanel> createPanel(const std::string& id, const std::string& title = "");
    void destroyPanel(const std::string& id);
    
    // Spatial interaction
    void handleSpatialInteraction(const core::SphericalCoords& interaction_point);
    void update(double delta_time);
    void render(rendering::SoftwareRenderer& renderer);
    
    // Layout management
    void arrangeControls(const std::string& layout_type = "radial");
    void snapToGrid();
    void alignControls(const std::string& alignment = "center");
    
    // Workspace management
    void saveWorkspace(const std::string& name);
    void loadWorkspace(const std::string& name);
    void clearWorkspace();

private:
    SpatialInterfaceManager() = default;
    
    std::map<std::string, std::shared_ptr<SpatialElement>> controls_;
    std::map<std::string, std::shared_ptr<SpatialPanel>> panels_;
    std::shared_ptr<SpatialElement> focused_control_;
    std::shared_ptr<SpatialElement> hovered_control_;
    
    // Workspace state
    struct WorkspaceState {
        std::map<std::string, SpatialElement::SpatialState> control_states;
        std::vector<std::string> control_order;
    };
    std::map<std::string, WorkspaceState> saved_workspaces_;
};

/**
 * Spatial Interface Presets - Pre-configured 3D GUI layouts
 */
namespace SpatialPresets {
    
    // Development workspace preset
    std::shared_ptr<SpatialPanel> createDevelopmentWorkspace();
    
    // Material editing preset
    std::shared_ptr<SpatialPanel> createMaterialEditor();
    
    // Testing framework preset
    std::shared_ptr<SpatialPanel> createTestingWorkspace();
    
    // Performance monitoring preset
    std::shared_ptr<SpatialPanel> createPerformanceMonitor();
    
    // Debug tools preset
    std::shared_ptr<SpatialPanel> createDebugTools();
    
    // Plugin management preset
    std::shared_ptr<SpatialPanel> createPluginManager();
}

} // namespace hsml::gui
