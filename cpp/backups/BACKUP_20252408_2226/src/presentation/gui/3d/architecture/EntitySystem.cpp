/** @file EntitySystem.h
 * @brief Futureproofed Entity System Implementation
 *
 * Phase 6: Perfect Architecture - Concrete Entity Implementation
 * Thread-safe, memory-efficient entity management with ECS architecture.
 */

#pragma once

#include "src/presentation/gui/3d/architecture/ArchitectureCore.h"
#include <vector>
#include <unordered_map>
#include <memory>
#include <shared_mutex>
#include <atomic>
#include <string>

namespace hsml {
namespace gui3d {
namespace architecture {

/**
 * @brief Concrete Entity Implementation
 *
 * Futureproofed entity with:
 * - Dynamic component management
 * - Thread-safe operations
 * - Memory-efficient storage
 * - Serialization support
 * - Event handling
 */
class Entity : public IEntity, public std::enable_shared_from_this<Entity> {
public:
    explicit Entity(const std::string& name = "");
    ~Entity() override;

    // IEntity implementation
    uint64_t get_id() const override { return id_; }
    const std::string& get_name() const override { return name_; }
    void set_name(const std::string& name) override { name_ = name; }

    bool has_component(const std::string& type) const override;
    IComponent* get_component(const std::string& type) override;
    const IComponent* get_component(const std::string& type) const override;

    bool add_component(std::unique_ptr<IComponent> component) override;
    bool remove_component(const std::string& type) override;
    void clear_components() override;

    std::vector<std::string> get_component_types() const override;
    size_t get_component_count() const override { return components_.size(); }

    bool is_active() const override { return active_; }
    void set_active(bool active) override { active_ = active; }

    bool is_visible() const override { return visible_; }
    void set_visible(bool visible) override { visible_ = visible; }

    void initialize() override;
    void update(double delta_time) override;
    void render() override;
    void shutdown() override;

    std::string serialize() const override;
    bool deserialize(const std::string& data) override;

    void send_event(std::shared_ptr<IEvent> event) override;
    void receive_event(std::shared_ptr<IEvent> event) override;

    // Extended functionality
    template<typename T>
    T* get_component() {
        return dynamic_cast<T*>(get_component(T::COMPONENT_TYPE));
    }

    template<typename T>
    const T* get_component() const {
        return dynamic_cast<const T*>(get_component(T::COMPONENT_TYPE));
    }

    template<typename T>
    bool add_component() {
        return add_component(std::make_unique<T>());
    }

    // Performance monitoring
    size_t get_memory_usage() const;
    void optimize_memory();

    // Component iteration
    template<typename Func>
    void for_each_component(Func&& func) {
        std::shared_lock<std::shared_mutex> lock(components_mutex_);
        for (const auto& [type, component] : components_) {
            func(component.get());
        }
    }

    template<typename Func>
    void for_each_component(Func&& func) const {
        std::shared_lock<std::shared_mutex> lock(components_mutex_);
        for (const auto& [type, component] : components_) {
            func(component.get());
        }
    }

private:
    uint64_t id_;
    std::string name_;
    bool active_{true};
    bool visible_{true};

    // Component storage with thread safety
    std::unordered_map<std::string, std::unique_ptr<IComponent>> components_;
    mutable std::shared_mutex components_mutex_;

    // Event handling
    std::vector<std::shared_ptr<IEvent>> pending_events_;
    mutable std::mutex events_mutex_;

    // Component type registry for validation
    static const std::unordered_set<std::string> VALID_COMPONENT_TYPES;

    // Helper methods
    bool is_valid_component_type(const std::string& type) const;
    void process_pending_events();
    void notify_component_added(const std::string& type);
    void notify_component_removed(const std::string& type);

    // Memory optimization
    void defragment_components();
    void remove_unused_components();
};

/**
 * @brief Component Base Class
 *
 * Provides common functionality for all components with futureproofing.
 */
class ComponentBase : public IComponent {
public:
    explicit ComponentBase(const std::string& type) : type_(type) {}
    ~ComponentBase() override = default;

    const std::string& get_type() const override { return type_; }
    uint64_t get_entity_id() const override { return entity_id_; }
    void set_entity_id(uint64_t id) override { entity_id_ = id; }

    bool is_enabled() const override { return enabled_; }
    void set_enabled(bool enabled) override { enabled_ = enabled; }

    void initialize() override {}
    void update(double delta_time) override {}
    void render() override {}
    void shutdown() override {}

    std::string serialize() const override;
    bool deserialize(const std::string& data) override;

    size_t get_memory_usage() const override { return sizeof(*this); }
    void optimize_memory() override {}

protected:
    std::string type_;
    uint64_t entity_id_{0};
    bool enabled_{true};

    // Serialization helpers
    virtual void serialize_data(std::ostream& stream) const {}
    virtual void deserialize_data(std::istream& stream) {}
};

/**
 * @brief Transform Component - Position, Rotation, Scale in 3D Space
 */
class TransformComponent : public ComponentBase {
public:
    static const std::string COMPONENT_TYPE;

    TransformComponent();
    ~TransformComponent() override = default;

    // Position
    const SphericalCoords& get_position() const { return position_; }
    void set_position(const SphericalCoords& position) { position_ = position; }

    // Rotation (Euler angles in radians)
    const Vector3& get_rotation() const { return rotation_; }
    void set_rotation(const Vector3& rotation) { rotation_ = rotation; }

    // Scale
    const Vector3& get_scale() const { return scale_; }
    void set_scale(const Vector3& scale) { scale_ = scale; }

    // Transform operations
    void translate(const SphericalCoords& delta);
    void rotate(const Vector3& delta);
    void scale_by(const Vector3& factor);

    // Transform matrix (for rendering) - simplified for this implementation
    void get_transform_data(float* position, float* rotation, float* scale) const;

    // Component overrides
    void update(double delta_time) override;
    std::string serialize() const override;
    bool deserialize(const std::string& data) override;
    size_t get_memory_usage() const override;

protected:
    void serialize_data(std::ostream& stream) const override;
    void deserialize_data(std::istream& stream) override;

private:
    SphericalCoords position_{0.0, 0.0, 0.0};
    Vector3 rotation_{0.0f, 0.0f, 0.0f};
    Vector3 scale_{1.0f, 1.0f, 1.0f};

    // Animation state
    SphericalCoords target_position_;
    Vector3 target_rotation_;
    Vector3 target_scale_;
    float animation_speed_{5.0f};
    bool animating_{false};
};

/**
 * @brief Visual Component - Rendering and Appearance
 */
class VisualComponent : public ComponentBase {
public:
    static const std::string COMPONENT_TYPE;

    VisualComponent();
    ~VisualComponent() override = default;

    // Appearance
    const Vector3& get_color() const { return color_; }
    void set_color(const Vector3& color) { color_ = color; }

    float get_opacity() const { return opacity_; }
    void set_opacity(float opacity) { opacity_ = std::clamp(opacity, 0.0f, 1.0f); }

    // Glow effect
    float get_glow_intensity() const { return glow_intensity_; }
    void set_glow_intensity(float intensity) { glow_intensity_ = std::max(0.0f, intensity); }

    const Vector3& get_glow_color() const { return glow_color_; }
    void set_glow_color(const Vector3& color) { glow_color_ = color; }

    // Size and shape
    float get_radius() const { return radius_; }
    void set_radius(float radius) { radius_ = std::max(0.1f, radius); }

    // Visual effects
    bool is_selected() const { return selected_; }
    void set_selected(bool selected) { selected_ = selected; }

    bool is_hovered() const { return hovered_; }
    void set_hovered(bool hovered) { hovered_ = hovered; }

    // Animation
    void start_pulse_animation(float frequency = 1.0f, float amplitude = 0.3f);
    void stop_animation();
    bool is_animating() const { return animating_; }

    // Component overrides
    void update(double delta_time) override;
    void render() override;
    std::string serialize() const override;
    bool deserialize(const std::string& data) override;
    size_t get_memory_usage() const override;

protected:
    void serialize_data(std::ostream& stream) const override;
    void deserialize_data(std::istream& stream) override;

private:
    Vector3 color_{1.0f, 1.0f, 1.0f};
    float opacity_{1.0f};
    float glow_intensity_{0.0f};
    Vector3 glow_color_{1.0f, 1.0f, 0.5f};
    float radius_{1.0f};

    bool selected_{false};
    bool hovered_{false};

    // Animation state
    bool animating_{false};
    float pulse_frequency_{1.0f};
    float pulse_amplitude_{0.3f};
    float animation_time_{0.0f};
    Vector3 original_color_;
    float original_glow_;

    void update_animation(double delta_time);
};

/**
 * @brief Code Component - Source Code Management
 */
class CodeComponent : public ComponentBase {
public:
    static const std::string COMPONENT_TYPE;

    CodeComponent();
    ~CodeComponent() override = default;

    // File information
    const std::string& get_filename() const { return filename_; }
    void set_filename(const std::string& filename) { filename_ = filename; }

    const std::string& get_file_path() const { return file_path_; }
    void set_file_path(const std::string& path) { file_path_ = path; }

    // Content management
    const std::string& get_content() const { return content_; }
    void set_content(const std::string& content);

    const std::vector<std::string>& get_lines() const { return lines_; }
    std::string get_line(size_t line_number) const;
    size_t get_line_count() const { return lines_.size(); }

    // Language detection
    std::string get_language() const { return language_; }
    void set_language(const std::string& language) { language_ = language; }

    // Syntax highlighting
    bool is_syntax_highlighting_enabled() const { return syntax_highlighting_; }
    void set_syntax_highlighting_enabled(bool enabled) { syntax_highlighting_ = enabled; }

    // Code analysis
    std::vector<std::pair<std::string, size_t>> find_functions() const;
    std::vector<std::pair<std::string, size_t>> find_classes() const;
    std::vector<std::pair<std::string, size_t>> find_includes() const;

    // Line highlighting
    void highlight_line(size_t line_number, const Vector3& color = Vector3(1.0f, 1.0f, 0.0f));
    void clear_highlights();
    const std::vector<std::pair<size_t, Vector3>>& get_highlights() const { return highlights_; }

    // Component overrides
    void update(double delta_time) override;
    std::string serialize() const override;
    bool deserialize(const std::string& data) override;
    size_t get_memory_usage() const override;

protected:
    void serialize_data(std::ostream& stream) const override;
    void deserialize_data(std::istream& stream) override;

private:
    std::string filename_;
    std::string file_path_;
    std::string content_;
    std::vector<std::string> lines_;
    std::string language_;
    bool syntax_highlighting_{true};

    std::vector<std::pair<size_t, Vector3>> highlights_;

    // Helper methods
    void parse_content();
    void detect_language();
    std::string extract_function_name(const std::string& line) const;
    std::string extract_class_name(const std::string& line) const;
    std::string extract_include_name(const std::string& line) const;
};

/**
 * @brief Connection Component - Manages Relationships
 */
class ConnectionComponent : public ComponentBase {
public:
    static const std::string COMPONENT_TYPE;

    ConnectionComponent();
    ~ConnectionComponent() override = default;

    // Connection endpoints
    uint64_t get_target_entity_id() const { return target_entity_id_; }
    void set_target_entity_id(uint64_t id) { target_entity_id_ = id; }

    // Connection properties
    ConnectionType get_connection_type() const { return connection_type_; }
    void set_connection_type(ConnectionType type) { connection_type_ = type; }

    const std::string& get_description() const { return description_; }
    void set_description(const std::string& desc) { description_ = desc; }

    // Visual properties
    const Vector3& get_color() const { return color_; }
    void set_color(const Vector3& color) { color_ = color; }

    float get_thickness() const { return thickness_; }
    void set_thickness(float thickness) { thickness_ = std::max(0.1f, thickness); }

    // Animation
    void set_pulse_frequency(float frequency) { pulse_frequency_ = frequency; }
    float get_pulse_frequency() const { return pulse_frequency_; }

    void set_pulse_amplitude(float amplitude) { pulse_amplitude_ = amplitude; }
    float get_pulse_amplitude() const { return pulse_amplitude_; }

    // Activity monitoring
    void set_activity_level(float level) { activity_level_ = std::clamp(level, 0.0f, 1.0f); }
    float get_activity_level() const { return activity_level_; }

    // State
    bool is_active() const { return active_; }
    void set_active(bool active) { active_ = active; }

    // Component overrides
    void update(double delta_time) override;
    void render() override;
    std::string serialize() const override;
    bool deserialize(const std::string& data) override;
    size_t get_memory_usage() const override;

protected:
    void serialize_data(std::ostream& stream) const override;
    void deserialize_data(std::istream& stream) override;

private:
    uint64_t target_entity_id_{0};
    ConnectionType connection_type_{ConnectionType::FUNCTION_CALL};
    std::string description_;
    Vector3 color_{0.0f, 1.0f, 0.0f};
    float thickness_{2.0f};
    bool active_{true};

    // Animation
    float pulse_frequency_{0.0f};
    float pulse_amplitude_{0.0f};
    float animation_time_{0.0f};
    float activity_level_{0.0f};
};

/**
 * @brief Entity Factory - Futureproofed Entity Creation
 */
class EntityFactory {
public:
    static std::shared_ptr<Entity> create_code_entity(
        const std::string& filename,
        const std::string& content,
        const SphericalCoords& position
    );

    static std::shared_ptr<Entity> create_connection_entity(
        uint64_t from_entity_id,
        uint64_t to_entity_id,
        ConnectionType type,
        const std::string& description
    );

    static std::shared_ptr<Entity> create_camera_entity(
        const SphericalCoords& position
    );

    static std::shared_ptr<Entity> create_ui_entity(
        const std::string& ui_type,
        const SphericalCoords& position
    );

    // Template-based entity creation
    template<typename... Components>
    static std::shared_ptr<Entity> create_entity_with_components(
        const std::string& name,
        Components&&... components
    ) {
        auto entity = std::make_shared<Entity>(name);

        // Add each component
        (entity->add_component(std::forward<Components>(components)), ...);

        return entity;
    }

private:
    static uint64_t next_entity_id_;
};

} // namespace architecture
} // namespace gui3d
} // namespace hsml
