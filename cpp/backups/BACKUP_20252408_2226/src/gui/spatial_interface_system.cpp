/**
 * Spatial Interface System (SIS) - Implementation
 * NO CARTESIAN COORDINATES - PURE SPHERICAL REALITY
 * 
 * Revolutionary 3D-native GUI where every control is a spatial object
 * - Controls exist as 3D elements in spherical space
 * - Can be moved, rotated, scaled, collapsed, expanded
 * - No traditional 2D overlays - everything is spatial
 */

#include "hsml/gui/spatial_interface_system.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include "hsml/rendering/software_renderer.h"
#include "hsml/core/state_tensor.h"

#include <cmath>
#include <algorithm>
#include <iostream>

namespace hsml::gui {

// =============================================================================
// SpatialElement Implementation
// =============================================================================

SpatialElement::SpatialElement(const std::string& id, const SpatialState& state)
    : id_(id), state_(state) {
}

void SpatialElement::setPosition(const core::SphericalCoords& pos) {
    state_.position = pos;
}

void SpatialElement::setOrientation(const core::SphericalCoords& orient) {
    state_.orientation = orient;
}

void SpatialElement::setScale(double scale) {
    state_.scale = std::max(0.1, scale);
}

void SpatialElement::moveTo(const core::SphericalCoords& target, double duration) {
    animation_.start_pos = state_.position;
    animation_.end_pos = target;
    animation_.duration = duration;
    animation_.elapsed = 0.0;
    animation_.active = true;
}

void SpatialElement::rotateTo(const core::SphericalCoords& target, double duration) {
    // Implement smooth rotation animation
    // This would interpolate between current and target orientation
}

void SpatialElement::show() {
    state_.visible = true;
}

void SpatialElement::hide() {
    state_.visible = false;
}

void SpatialElement::collapse() {
    state_.collapsed = true;
    setScale(0.1); // Shrink to almost invisible
}

void SpatialElement::expand() {
    state_.collapsed = false;
    setScale(1.0); // Return to normal size
}

void SpatialElement::setInteractive(bool interactive) {
    state_.interactive = interactive;
}

void SpatialElement::attachTo(std::shared_ptr<SpatialElement> parent) {
    if (parent_) {
        parent_->children_.erase(
            std::remove(parent_->children_.begin(), parent_->children_.end(), shared_from_this()),
            parent_->children_.end()
        );
    }
    
    parent_ = parent;
    if (parent_) {
        parent_->children_.push_back(shared_from_this());
    }
}

void SpatialElement::detach() {
    if (parent_) {
        parent_->children_.erase(
            std::remove(parent_->children_.begin(), parent_->children_.end(), shared_from_this()),
            parent_->children_.end()
        );
        parent_.reset();
    }
}

std::vector<std::shared_ptr<SpatialElement>> SpatialElement::getChildren() const {
    return children_;
}

bool SpatialElement::handleSpatialInteraction(const core::SphericalCoords& interaction_point) {
    if (!state_.visible || !state_.interactive) {
        return false;
    }
    
    // Calculate distance in spherical space
    double distance = state_.position.angular_distance(interaction_point);
    double interaction_radius = 10.0 * state_.scale; // Interaction radius based on scale
    
    return distance <= interaction_radius;
}

void SpatialElement::onSpatialHover() {
    // Default hover behavior - could change color, scale, etc.
    state_.color = core::Vector3(1.2, 1.2, 1.2); // Brighten slightly
}

void SpatialElement::onSpatialClick() {
    // Default click behavior
}

void SpatialElement::onSpatialDrag(const core::SphericalCoords& delta) {
    // Default drag behavior - move the element
    state_.position.r += delta.r;
    state_.position.theta += delta.theta;
    state_.position.phi += delta.phi;
}

// =============================================================================
// SpatialButton Implementation
// =============================================================================

SpatialButton::SpatialButton(const std::string& id, const std::string& label)
    : SpatialElement(id), label_(label) {
}

void SpatialButton::setLabel(const std::string& label) {
    label_ = label;
}

void SpatialButton::setAction(std::function<void()> action) {
    action_ = action;
}

void SpatialButton::setButtonStyle(const std::string& style) {
    style_ = style;
}

void SpatialButton::render(rendering::SoftwareRenderer& renderer) {
    if (!state_.visible) return;
    
    // Render button as a 3D sphere or cube in spherical space
    double button_radius = 5.0 * state_.scale;
    
    // Convert spherical position to renderer coordinates
    auto render_pos = renderer.world_to_screen_spherical(state_.position);
    
    // Render based on style
    if (style_ == "flat") {
        renderer.render_sphere_steradian(state_.position, button_radius, state_.color);
    } else if (style_ == "raised") {
        // Render with depth effect
        auto raised_pos = state_.position;
        raised_pos.r += 2.0 * state_.scale;
        renderer.render_sphere_steradian(raised_pos, button_radius, state_.color);
    } else if (style_ == "glowing") {
        // Render with glow effect
        renderer.render_sphere_steradian(state_.position, button_radius * 1.2, 
                                       core::Vector3(0.5, 0.5, 1.0)); // Glow
        renderer.render_sphere_steradian(state_.position, button_radius, state_.color);
    }
    
    // Render label if present
    if (!label_.empty()) {
        // Render 3D text label (simplified - would need proper text rendering)
        auto label_pos = state_.position;
        label_pos.r += button_radius + 2.0;
        // renderer.render_text_3d(label_pos, label_, state_.color);
    }
}

void SpatialButton::update(double delta_time) {
    SpatialElement::update(delta_time);
    
    // Update press animation
    if (pressed_) {
        press_animation_ += delta_time * 8.0; // Fast press animation
        if (press_animation_ >= 1.0) {
            press_animation_ = 1.0;
            pressed_ = false;
        }
    } else {
        press_animation_ -= delta_time * 4.0; // Slower release animation
        if (press_animation_ <= 0.0) {
            press_animation_ = 0.0;
        }
    }
    
    // Apply press animation to scale
    double press_scale = 1.0 - (press_animation_ * 0.1);
    setScale(state_.scale * press_scale);
}

bool SpatialButton::handleSpatialInteraction(const core::SphericalCoords& interaction_point) {
    if (!SpatialElement::handleSpatialInteraction(interaction_point)) {
        return false;
    }
    
    // Button-specific interaction logic
    return true;
}

void SpatialButton::onSpatialClick() {
    pressed_ = true;
    press_animation_ = 0.0;
    
    if (action_) {
        action_();
    }
}

// =============================================================================
// SpatialSlider Implementation
// =============================================================================

SpatialSlider::SpatialSlider(const std::string& id, double min_val, double max_val)
    : SpatialElement(id), min_value_(min_val), max_value_(max_val), current_value_((min_val + max_val) / 2.0) {
}

void SpatialSlider::setValue(double value) {
    current_value_ = std::clamp(value, min_value_, max_value_);
    if (on_value_changed_) {
        on_value_changed_(current_value_);
    }
}

void SpatialSlider::setRange(double min_val, double max_val) {
    min_value_ = min_val;
    max_value_ = max_val;
    setValue(current_value_); // Clamp to new range
}

void SpatialSlider::setOnValueChanged(std::function<void(double)> callback) {
    on_value_changed_ = callback;
}

void SpatialSlider::setSliderStyle(const std::string& style) {
    style_ = style;
}

void SpatialSlider::render(rendering::SoftwareRenderer& renderer) {
    if (!state_.visible) return;
    
    double slider_length = 20.0 * state_.scale;
    double slider_radius = 1.0 * state_.scale;
    
    if (style_ == "linear") {
        // Render linear slider as a cylinder in spherical space
        auto start_pos = state_.position;
        auto end_pos = state_.position;
        end_pos.phi += 0.2; // Extend in phi direction
        
        // Render slider track
        renderer.render_sphere_steradian(start_pos, slider_radius, 
                                       core::Vector3(0.3, 0.3, 0.3));
        renderer.render_sphere_steradian(end_pos, slider_radius, 
                                       core::Vector3(0.3, 0.3, 0.3));
        
        // Render slider handle at current value position
        double t = (current_value_ - min_value_) / (max_value_ - min_value_);
        auto handle_pos = start_pos;
        handle_pos.phi += t * 0.2;
        renderer.render_sphere_steradian(handle_pos, slider_radius * 1.5, state_.color);
        
    } else if (style_ == "radial") {
        // Render radial slider as an arc
        double start_angle = 0.0;
        double end_angle = M_PI; // 180 degrees
        
        double t = (current_value_ - min_value_) / (max_value_ - min_value_);
        double current_angle = start_angle + t * (end_angle - start_angle);
        
        // Render arc track
        for (int i = 0; i <= 10; ++i) {
            double angle = start_angle + (i / 10.0) * (end_angle - start_angle);
            auto track_pos = state_.position;
            track_pos.theta += angle;
            renderer.render_sphere_steradian(track_pos, slider_radius, 
                                           core::Vector3(0.3, 0.3, 0.3));
        }
        
        // Render handle at current angle
        auto handle_pos = state_.position;
        handle_pos.theta += current_angle;
        renderer.render_sphere_steradian(handle_pos, slider_radius * 1.5, state_.color);
    }
}

void SpatialSlider::update(double delta_time) {
    SpatialElement::update(delta_time);
    
    // Update dragging state
    if (dragging_) {
        // Handle drag updates
    }
}

bool SpatialSlider::handleSpatialInteraction(const core::SphericalCoords& interaction_point) {
    if (!SpatialElement::handleSpatialInteraction(interaction_point)) {
        return false;
    }
    
    // Slider-specific interaction logic
    return true;
}

void SpatialSlider::onSpatialDrag(const core::SphericalCoords& delta) {
    if (!dragging_) {
        dragging_ = true;
        drag_start_ = state_.position;
    }
    
    // Update value based on drag
    if (style_ == "linear") {
        double drag_distance = delta.phi;
        double value_range = max_value_ - min_value_;
        double new_value = current_value_ + (drag_distance * value_range / 0.2);
        setValue(new_value);
    } else if (style_ == "radial") {
        double drag_angle = delta.theta;
        double value_range = max_value_ - min_value_;
        double new_value = current_value_ + (drag_angle * value_range / M_PI);
        setValue(new_value);
    }
}

// =============================================================================
// SpatialToggle Implementation
// =============================================================================

SpatialToggle::SpatialToggle(const std::string& id, bool initial_state)
    : SpatialElement(id), state_(initial_state) {
}

void SpatialToggle::setState(bool state) {
    if (state_ != state) {
        state_ = state;
        if (on_state_changed_) {
            on_state_changed_(state_);
        }
    }
}

void SpatialToggle::setOnStateChanged(std::function<void(bool)> callback) {
    on_state_changed_ = callback;
}

void SpatialToggle::setToggleStyle(const std::string& style) {
    style_ = style;
}

void SpatialToggle::render(rendering::SoftwareRenderer& renderer) {
    if (!state_.visible) return;
    
    double toggle_radius = 3.0 * state_.scale;
    
    if (style_ == "switch") {
        // Render as a switch with two positions
        auto off_pos = state_.position;
        auto on_pos = state_.position;
        on_pos.phi += 0.1; // Offset for on position
        
        // Render track
        renderer.render_sphere_steradian(off_pos, toggle_radius, 
                                       core::Vector3(0.3, 0.3, 0.3));
        renderer.render_sphere_steradian(on_pos, toggle_radius, 
                                       core::Vector3(0.3, 0.3, 0.3));
        
        // Render handle at current state
        auto handle_pos = state_ ? on_pos : off_pos;
        core::Vector3 handle_color = state_ ? core::Vector3(0.2, 1.0, 0.2) : core::Vector3(0.8, 0.8, 0.8);
        renderer.render_sphere_steradian(handle_pos, toggle_radius * 0.8, handle_color);
        
    } else if (style_ == "checkbox") {
        // Render as a checkbox
        core::Vector3 box_color = state_ ? core::Vector3(0.2, 1.0, 0.2) : core::Vector3(0.8, 0.8, 0.8);
        renderer.render_sphere_steradian(state_.position, toggle_radius, box_color);
        
        // Render checkmark if checked
        if (state_) {
            // Simplified checkmark rendering
            auto check_pos = state_.position;
            check_pos.r += toggle_radius * 0.3;
            renderer.render_sphere_steradian(check_pos, toggle_radius * 0.3, 
                                           core::Vector3(1.0, 1.0, 1.0));
        }
    }
}

void SpatialToggle::update(double delta_time) {
    SpatialElement::update(delta_time);
    
    // Update animation
    double target_animation = state_ ? 1.0 : 0.0;
    double animation_speed = 8.0;
    
    if (animation_progress_ < target_animation) {
        animation_progress_ += delta_time * animation_speed;
        if (animation_progress_ > target_animation) {
            animation_progress_ = target_animation;
        }
    } else if (animation_progress_ > target_animation) {
        animation_progress_ -= delta_time * animation_speed;
        if (animation_progress_ < target_animation) {
            animation_progress_ = target_animation;
        }
    }
}

bool SpatialToggle::handleSpatialInteraction(const core::SphericalCoords& interaction_point) {
    if (!SpatialElement::handleSpatialInteraction(interaction_point)) {
        return false;
    }
    
    return true;
}

void SpatialToggle::onSpatialClick() {
    setState(!state_);
}

// =============================================================================
// SpatialPanel Implementation
// =============================================================================

SpatialPanel::SpatialPanel(const std::string& id, const std::string& title)
    : SpatialElement(id), title_(title) {
}

void SpatialPanel::setTitle(const std::string& title) {
    title_ = title;
}

void SpatialPanel::addControl(std::shared_ptr<SpatialElement> control) {
    controls_.push_back(control);
    control->attachTo(shared_from_this());
    
    // Arrange controls based on layout
    arrangeControls();
}

void SpatialPanel::removeControl(const std::string& control_id) {
    controls_.erase(
        std::remove_if(controls_.begin(), controls_.end(),
                      [&](const std::shared_ptr<SpatialElement>& control) {
                          return control->getId() == control_id;
                      }),
        controls_.end()
    );
}

void SpatialPanel::setPanelStyle(const std::string& style) {
    style_ = style;
}

void SpatialPanel::setLayout(const std::string& layout) {
    layout_ = layout;
    arrangeControls();
}

void SpatialPanel::render(rendering::SoftwareRenderer& renderer) {
    if (!state_.visible) return;
    
    // Render panel background
    double panel_radius = 15.0 * state_.scale;
    renderer.render_sphere_steradian(state_.position, panel_radius, 
                                   core::Vector3(0.1, 0.1, 0.15));
    
    // Render title if present
    if (!title_.empty() && title_bar_visible_) {
        auto title_pos = state_.position;
        title_pos.r += panel_radius + 2.0;
        // renderer.render_text_3d(title_pos, title_, core::Vector3(1.0, 1.0, 1.0));
    }
    
    // Render all child controls
    for (auto& control : controls_) {
        control->render(renderer);
    }
}

void SpatialPanel::update(double delta_time) {
    SpatialElement::update(delta_time);
    
    // Update panel animation
    if (state_.collapsed) {
        panel_animation_ -= delta_time * 4.0;
        if (panel_animation_ < 0.0) panel_animation_ = 0.0;
    } else {
        panel_animation_ += delta_time * 4.0;
        if (panel_animation_ > 1.0) panel_animation_ = 1.0;
    }
    
    // Update all child controls
    for (auto& control : controls_) {
        control->update(delta_time);
    }
}

bool SpatialPanel::handleSpatialInteraction(const core::SphericalCoords& interaction_point) {
    if (!SpatialElement::handleSpatialInteraction(interaction_point)) {
        return false;
    }
    
    // Check child controls first (they're on top)
    for (auto& control : controls_) {
        if (control->handleSpatialInteraction(interaction_point)) {
            return true;
        }
    }
    
    return true; // Panel itself is interactive
}

void SpatialPanel::arrangeControls() {
    if (controls_.empty()) return;
    
    double spacing = 8.0 * state_.scale;
    double start_angle = 0.0;
    
    if (layout_ == "vertical") {
        // Arrange controls vertically in spherical space
        for (size_t i = 0; i < controls_.size(); ++i) {
            auto& control = controls_[i];
            auto pos = state_.position;
            pos.theta += (i * spacing) / 100.0; // Adjust theta for vertical spacing
            control->setPosition(pos);
        }
        
    } else if (layout_ == "horizontal") {
        // Arrange controls horizontally in spherical space
        for (size_t i = 0; i < controls_.size(); ++i) {
            auto& control = controls_[i];
            auto pos = state_.position;
            pos.phi += (i * spacing) / 100.0; // Adjust phi for horizontal spacing
            control->setPosition(pos);
        }
        
    } else if (layout_ == "radial") {
        // Arrange controls in a circle around the panel center
        double angle_step = 2.0 * M_PI / controls_.size();
        for (size_t i = 0; i < controls_.size(); ++i) {
            auto& control = controls_[i];
            auto pos = state_.position;
            pos.phi += start_angle + (i * angle_step);
            pos.r += 12.0 * state_.scale; // Offset from center
            control->setPosition(pos);
        }
        
    } else if (layout_ == "grid") {
        // Arrange controls in a grid pattern
        int cols = static_cast<int>(std::sqrt(controls_.size()));
        int rows = (controls_.size() + cols - 1) / cols;
        
        for (size_t i = 0; i < controls_.size(); ++i) {
            auto& control = controls_[i];
            int row = i / cols;
            int col = i % cols;
            
            auto pos = state_.position;
            pos.phi += (col * spacing) / 100.0;
            pos.theta += (row * spacing) / 100.0;
            control->setPosition(pos);
        }
    }
}

// =============================================================================
// SpatialInterfaceManager Implementation
// =============================================================================

SpatialInterfaceManager& SpatialInterfaceManager::instance() {
    static SpatialInterfaceManager instance;
    return instance;
}

std::shared_ptr<SpatialElement> SpatialInterfaceManager::createControl(const std::string& type, const std::string& id) {
    std::shared_ptr<SpatialElement> control;
    
    if (type == "button") {
        control = std::make_shared<SpatialButton>(id);
    } else if (type == "slider") {
        control = std::make_shared<SpatialSlider>(id);
    } else if (type == "toggle") {
        control = std::make_shared<SpatialToggle>(id);
    } else if (type == "dropdown") {
        control = std::make_shared<SpatialDropdown>(id);
    } else if (type == "textinput") {
        control = std::make_shared<SpatialTextInput>(id);
    } else if (type == "colorpicker") {
        control = std::make_shared<SpatialColorPicker>(id);
    } else if (type == "panel") {
        control = std::make_shared<SpatialPanel>(id);
    } else if (type == "viewport") {
        control = std::make_shared<SpatialViewport>(id);
    }
    
    if (control) {
        controls_[id] = control;
    }
    
    return control;
}

void SpatialInterfaceManager::destroyControl(const std::string& id) {
    controls_.erase(id);
}

std::shared_ptr<SpatialElement> SpatialInterfaceManager::getControl(const std::string& id) {
    auto it = controls_.find(id);
    return (it != controls_.end()) ? it->second : nullptr;
}

std::shared_ptr<SpatialPanel> SpatialInterfaceManager::createPanel(const std::string& id, const std::string& title) {
    auto panel = std::make_shared<SpatialPanel>(id, title);
    panels_[id] = panel;
    controls_[id] = panel;
    return panel;
}

void SpatialInterfaceManager::destroyPanel(const std::string& id) {
    panels_.erase(id);
    controls_.erase(id);
}

void SpatialInterfaceManager::handleSpatialInteraction(const core::SphericalCoords& interaction_point) {
    // Find the topmost interactive element
    std::shared_ptr<SpatialElement> top_element = nullptr;
    double min_distance = std::numeric_limits<double>::max();
    
    for (auto& [id, control] : controls_) {
        if (control->isVisible() && control->isInteractive()) {
            double distance = control->getState().position.angular_distance(interaction_point);
            if (distance < min_distance) {
                min_distance = distance;
                top_element = control;
            }
        }
    }
    
    // Handle interaction with top element
    if (top_element) {
        if (top_element != hovered_control_) {
            if (hovered_control_) {
                // Exit hover state
            }
            hovered_control_ = top_element;
            top_element->onSpatialHover();
        }
        
        // Handle click/drag based on interaction type
        // This would be determined by input system
    }
}

void SpatialInterfaceManager::update(double delta_time) {
    for (auto& [id, control] : controls_) {
        control->update(delta_time);
    }
}

void SpatialInterfaceManager::render(rendering::SoftwareRenderer& renderer) {
    for (auto& [id, control] : controls_) {
        control->render(renderer);
    }
}

void SpatialInterfaceManager::arrangeControls(const std::string& layout_type) {
    // Arrange all controls in the specified layout
    std::vector<std::shared_ptr<SpatialElement>> root_controls;
    
    for (auto& [id, control] : controls_) {
        if (!control->getState().parent) {
            root_controls.push_back(control);
        }
    }
    
    if (layout_type == "radial") {
        double angle_step = 2.0 * M_PI / root_controls.size();
        for (size_t i = 0; i < root_controls.size(); ++i) {
            auto& control = root_controls[i];
            auto pos = control->getState().position;
            pos.phi = i * angle_step;
            pos.r = 50.0; // Fixed radius
            control->setPosition(pos);
        }
    }
}

void SpatialInterfaceManager::saveWorkspace(const std::string& name) {
    WorkspaceState state;
    
    for (auto& [id, control] : controls_) {
        state.control_states[id] = control->getState();
        state.control_order.push_back(id);
    }
    
    saved_workspaces_[name] = state;
}

void SpatialInterfaceManager::loadWorkspace(const std::string& name) {
    auto it = saved_workspaces_.find(name);
    if (it != saved_workspaces_.end()) {
        const auto& state = it->second;
        
        for (const auto& id : state.control_order) {
            auto control = getControl(id);
            if (control) {
                auto state_it = state.control_states.find(id);
                if (state_it != state.control_states.end()) {
                    // Restore control state
                    // This would require setters for all state properties
                }
            }
        }
    }
}

void SpatialInterfaceManager::clearWorkspace() {
    controls_.clear();
    panels_.clear();
    focused_control_.reset();
    hovered_control_.reset();
}

// =============================================================================
// SpatialPresets Implementation
// =============================================================================

namespace SpatialPresets {

std::shared_ptr<SpatialPanel> createDevelopmentWorkspace() {
    auto& manager = SpatialInterfaceManager::instance();
    auto panel = manager.createPanel("dev_workspace", "Development Workspace");
    
    // Add development tools
    auto compile_btn = std::make_shared<SpatialButton>("compile_btn", "Compile");
    auto run_btn = std::make_shared<SpatialButton>("run_btn", "Run");
    auto debug_btn = std::make_shared<SpatialButton>("debug_btn", "Debug");
    
    panel->addControl(compile_btn);
    panel->addControl(run_btn);
    panel->addControl(debug_btn);
    
    panel->setLayout("radial");
    return panel;
}

std::shared_ptr<SpatialPanel> createMaterialEditor() {
    auto& manager = SpatialInterfaceManager::instance();
    auto panel = manager.createPanel("material_editor", "Material Editor");
    
    // Add material controls
    auto albedo_picker = std::make_shared<SpatialColorPicker>("albedo_picker");
    auto metallic_slider = std::make_shared<SpatialSlider>("metallic_slider", 0.0, 1.0);
    auto roughness_slider = std::make_shared<SpatialSlider>("roughness_slider", 0.0, 1.0);
    auto emission_toggle = std::make_shared<SpatialToggle>("emission_toggle");
    
    panel->addControl(albedo_picker);
    panel->addControl(metallic_slider);
    panel->addControl(roughness_slider);
    panel->addControl(emission_toggle);
    
    panel->setLayout("vertical");
    return panel;
}

std::shared_ptr<SpatialPanel> createTestingWorkspace() {
    auto& manager = SpatialInterfaceManager::instance();
    auto panel = manager.createPanel("testing_workspace", "Testing Framework");
    
    // Add testing controls
    auto run_tests_btn = std::make_shared<SpatialButton>("run_tests_btn", "Run Tests");
    auto stop_tests_btn = std::make_shared<SpatialButton>("stop_tests_btn", "Stop Tests");
    auto coverage_toggle = std::make_shared<SpatialToggle>("coverage_toggle");
    
    panel->addControl(run_tests_btn);
    panel->addControl(stop_tests_btn);
    panel->addControl(coverage_toggle);
    
    panel->setLayout("horizontal");
    return panel;
}

std::shared_ptr<SpatialPanel> createPerformanceMonitor() {
    auto& manager = SpatialInterfaceManager::instance();
    auto panel = manager.createPanel("perf_monitor", "Performance Monitor");
    
    // Add performance controls
    auto fps_slider = std::make_shared<SpatialSlider>("fps_slider", 0.0, 120.0);
    auto memory_slider = std::make_shared<SpatialSlider>("memory_slider", 0.0, 100.0);
    auto profiling_toggle = std::make_shared<SpatialToggle>("profiling_toggle");
    
    panel->addControl(fps_slider);
    panel->addControl(memory_slider);
    panel->addControl(profiling_toggle);
    
    panel->setLayout("grid");
    return panel;
}

std::shared_ptr<SpatialPanel> createDebugTools() {
    auto& manager = SpatialInterfaceManager::instance();
    auto panel = manager.createPanel("debug_tools", "Debug Tools");
    
    // Add debug controls
    auto wireframe_toggle = std::make_shared<SpatialToggle>("wireframe_toggle");
    auto bounding_box_toggle = std::make_shared<SpatialToggle>("bounding_box_toggle");
    auto coordinate_display_toggle = std::make_shared<SpatialToggle>("coord_display_toggle");
    
    panel->addControl(wireframe_toggle);
    panel->addControl(bounding_box_toggle);
    panel->addControl(coordinate_display_toggle);
    
    panel->setLayout("vertical");
    return panel;
}

std::shared_ptr<SpatialPanel> createPluginManager() {
    auto& manager = SpatialInterfaceManager::instance();
    auto panel = manager.createPanel("plugin_manager", "Plugin Manager");
    
    // Add plugin controls
    auto load_plugin_btn = std::make_shared<SpatialButton>("load_plugin_btn", "Load Plugin");
    auto unload_plugin_btn = std::make_shared<SpatialButton>("unload_plugin_btn", "Unload Plugin");
    auto plugin_list_dropdown = std::make_shared<SpatialDropdown>("plugin_list");
    
    panel->addControl(load_plugin_btn);
    panel->addControl(unload_plugin_btn);
    panel->addControl(plugin_list_dropdown);
    
    panel->setLayout("radial");
    return panel;
}

} // namespace SpatialPresets

} // namespace hsml::gui
