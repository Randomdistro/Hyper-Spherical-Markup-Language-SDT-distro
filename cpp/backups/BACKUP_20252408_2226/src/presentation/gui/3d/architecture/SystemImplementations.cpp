/** @file SystemImplementations.h
 * @brief Futureproofed System Implementations
 *
 * Phase 6: Perfect Architecture - Concrete System Implementations
 * Plugin-based systems for entity management, rendering, physics, and interactions.
 */

#pragma once

#include "src/presentation/gui/3d/architecture/ArchitectureCore.h"
#include <vector>
#include <unordered_map>
#include <memory>
#include <shared_mutex>
#include <atomic>
#include <chrono>

namespace hsml {
namespace gui3d {
namespace architecture {

/**
 * @brief Rendering System - Futureproofed Graphics Pipeline
 *
 * Handles all visual rendering with support for:
 * - Multiple rendering backends
 * - Dynamic shader management
 * - Instanced rendering
 * - Level-of-detail (LOD)
 * - Post-processing effects
 */
class RenderingSystem : public ISystem {
public:
    static const std::string SYSTEM_NAME;

    RenderingSystem();
    ~RenderingSystem() override = default;

    // ISystem implementation
    const std::string& get_name() const override { return SYSTEM_NAME; }
    int get_priority() const override { return 100; } // Render last
    bool is_enabled() const override { return enabled_; }
    void set_enabled(bool enabled) override { enabled_ = enabled; }

    bool initialize() override;
    void update(double delta_time) override;
    void render() override;
    void shutdown() override;

    void on_entity_created(std::shared_ptr<IEntity> entity) override;
    void on_entity_destroyed(uint64_t entity_id) override;
    void on_component_added(uint64_t entity_id, const std::string& component_type) override;
    void on_component_removed(uint64_t entity_id, const std::string& component_type) override;

    std::vector<std::shared_ptr<IEntity>> query_entities(
        const std::function<bool(const IEntity&)>& predicate) const override;

    double get_average_update_time() const override { return avg_update_time_; }
    double get_average_render_time() const override { return avg_render_time_; }
    size_t get_processed_entities_count() const override { return processed_entities_; }

    // Rendering system specific functionality
    void set_rendering_backend(const std::string& backend);
    void set_target_frame_rate(double fps);
    void enable_feature(const std::string& feature);
    void disable_feature(const std::string& feature);
    bool is_feature_enabled(const std::string& feature) const;

    // Performance monitoring
    void begin_frame();
    void end_frame();
    double get_frame_rate() const { return current_frame_rate_; }
    double get_frame_time() const { return current_frame_time_; }

    // Camera and viewport management
    void set_camera_entity(uint64_t entity_id);
    void set_viewport_size(int width, int height);
    void set_field_of_view(double fov_degrees);

private:
    bool enabled_{true};
    double avg_update_time_{0.0};
    double avg_render_time_{0.0};
    size_t processed_entities_{0};

    // Rendering state
    std::string rendering_backend_{"opengl"};
    double target_frame_rate_{60.0};
    double current_frame_rate_{0.0};
    double current_frame_time_{0.0};

    // Camera and viewport
    uint64_t camera_entity_id_{0};
    int viewport_width_{1920};
    int viewport_height_{1080};
    double field_of_view_{60.0};

    // Performance tracking
    std::chrono::steady_clock::time_point frame_start_time_;
    std::vector<double> frame_times_;
    size_t max_frame_samples_{60};

    // Entity tracking
    std::unordered_map<uint64_t, std::shared_ptr<IEntity>> renderable_entities_;
    mutable std::shared_mutex entities_mutex_;

    // Rendering features
    std::unordered_set<std::string> enabled_features_;

    // Helper methods
    void update_renderable_entities();
    void sort_entities_by_depth();
    void batch_render_entities();
    void apply_post_processing();
    void update_performance_stats();
};

/**
 * @brief Physics System - Futureproofed Physics Simulation
 *
 * Handles all physical interactions with support for:
 * - Spherical coordinate physics
 * - Collision detection and response
 * - Force calculations
 * - Motion integration
 * - Multi-threaded simulation
 */
class PhysicsSystem : public ISystem {
public:
    static const std::string SYSTEM_NAME;

    PhysicsSystem();
    ~PhysicsSystem() override = default;

    // ISystem implementation
    const std::string& get_name() const override { return SYSTEM_NAME; }
    int get_priority() const override { return 50; }
    bool is_enabled() const override { return enabled_; }
    void set_enabled(bool enabled) override { enabled_ = enabled; }

    bool initialize() override;
    void update(double delta_time) override;
    void render() override {} // Physics system doesn't render
    void shutdown() override;

    void on_entity_created(std::shared_ptr<IEntity> entity) override;
    void on_entity_destroyed(uint64_t entity_id) override;
    void on_component_added(uint64_t entity_id, const std::string& component_type) override;
    void on_component_removed(uint64_t entity_id, const std::string& component_type) override;

    std::vector<std::shared_ptr<IEntity>> query_entities(
        const std::function<bool(const IEntity&)>& predicate) const override;

    double get_average_update_time() const override { return avg_update_time_; }
    double get_average_render_time() const override { return 0.0; }
    size_t get_processed_entities_count() const override { return processed_entities_; }

    // Physics system specific functionality
    void set_gravity_enabled(bool enabled);
    void set_collision_detection_enabled(bool enabled);
    void set_simulation_speed(float multiplier);
    void set_fixed_time_step(double seconds);

    // Force management
    void apply_force(uint64_t entity_id, const Vector3& force);
    void apply_impulse(uint64_t entity_id, const Vector3& impulse);
    void clear_forces(uint64_t entity_id);

    // Collision queries
    bool is_colliding(uint64_t entity_id) const;
    std::vector<uint64_t> get_colliding_entities(uint64_t entity_id) const;
    void resolve_collision(uint64_t entity1_id, uint64_t entity2_id);

private:
    bool enabled_{true};
    double avg_update_time_{0.0};
    size_t processed_entities_{0};

    // Physics settings
    bool gravity_enabled_{true};
    bool collision_detection_enabled_{true};
    float simulation_speed_{1.0f};
    double fixed_time_step_{1.0 / 60.0}; // 60 Hz
    double accumulator_{0.0};

    // Physical constants
    double gravity_acceleration_{9.81};
    double air_density_{1.225};
    double drag_coefficient_{0.47};

    // Entity tracking
    std::unordered_map<uint64_t, std::shared_ptr<IEntity>> physics_entities_;
    mutable std::shared_mutex entities_mutex_;

    // Force and impulse tracking
    struct AppliedForce {
        Vector3 force;
        double duration;
        std::chrono::steady_clock::time_point start_time;
    };
    std::unordered_map<uint64_t, std::vector<AppliedForce>> applied_forces_;
    mutable std::shared_mutex forces_mutex_;

    // Collision state
    struct CollisionInfo {
        uint64_t entity1;
        uint64_t entity2;
        Vector3 normal;
        double penetration;
        std::chrono::steady_clock::time_point timestamp;
    };
    std::vector<CollisionInfo> active_collisions_;
    mutable std::shared_mutex collisions_mutex_;

    // Helper methods
    void update_physics(double delta_time);
    void apply_gravity(std::shared_ptr<IEntity> entity, double delta_time);
    void apply_drag(std::shared_ptr<IEntity> entity, double delta_time);
    void apply_forces(std::shared_ptr<IEntity> entity, double delta_time);
    void integrate_motion(std::shared_ptr<IEntity> entity, double delta_time);
    void detect_collisions();
    void resolve_collisions();
    void update_collision_state();
};

/**
 * @brief Interaction System - Futureproofed User Interaction
 *
 * Handles all user interactions with support for:
 * - Mouse and keyboard input
 * - Touch and gesture recognition
 * - Voice commands
 * - Haptic feedback
 * - Multi-user collaboration
 */
class InteractionSystem : public ISystem {
public:
    static const std::string SYSTEM_NAME;

    InteractionSystem();
    ~InteractionSystem() override = default;

    // ISystem implementation
    const std::string& get_name() const override { return SYSTEM_NAME; }
    int get_priority() const override { return 10; } // Process input early
    bool is_enabled() const override { return enabled_; }
    void set_enabled(bool enabled) override { enabled_ = enabled; }

    bool initialize() override;
    void update(double delta_time) override;
    void render() override {} // Interaction system doesn't render
    void shutdown() override;

    void on_entity_created(std::shared_ptr<IEntity> entity) override;
    void on_entity_destroyed(uint64_t entity_id) override;
    void on_component_added(uint64_t entity_id, const std::string& component_type) override;
    void on_component_removed(uint64_t entity_id, const std::string& component_type) override;

    std::vector<std::shared_ptr<IEntity>> query_entities(
        const std::function<bool(const IEntity&)>& predicate) const override;

    double get_average_update_time() const override { return avg_update_time_; }
    double get_average_render_time() const override { return 0.0; }
    size_t get_processed_entities_count() const override { return processed_entities_; }

    // Interaction system specific functionality
    void set_mouse_sensitivity(float sensitivity);
    void set_keyboard_repeat_rate(float rate);
    void enable_gesture_recognition(bool enabled);
    void enable_voice_commands(bool enabled);
    void enable_haptic_feedback(bool enabled);

    // Input state queries
    bool is_key_pressed(int key_code) const;
    bool is_mouse_button_pressed(int button) const;
    void get_mouse_position(double& x, double& y) const;
    void get_mouse_delta(double& delta_x, double& delta_y) const;

    // Selection and focus
    uint64_t get_selected_entity() const { return selected_entity_id_; }
    uint64_t get_hovered_entity() const { return hovered_entity_id_; }
    void select_entity(uint64_t entity_id);
    void deselect_entity();
    void hover_entity(uint64_t entity_id);

    // Interaction callbacks
    using EntitySelectedCallback = std::function<void(uint64_t entity_id)>;
    using EntityHoveredCallback = std::function<void(uint64_t entity_id)>;
    using KeyPressedCallback = std::function<void(int key_code)>;
    using MouseClickedCallback = std::function<void(int button, double x, double y)>;

    void set_entity_selected_callback(EntitySelectedCallback callback);
    void set_entity_hovered_callback(EntityHoveredCallback callback);
    void set_key_pressed_callback(KeyPressedCallback callback);
    void set_mouse_clicked_callback(MouseClickedCallback callback);

private:
    bool enabled_{true};
    double avg_update_time_{0.0};
    size_t processed_entities_{0};

    // Input settings
    float mouse_sensitivity_{1.0f};
    float keyboard_repeat_rate_{30.0f};
    bool gesture_recognition_enabled_{false};
    bool voice_commands_enabled_{false};
    bool haptic_feedback_enabled_{false};

    // Selection state
    uint64_t selected_entity_id_{0};
    uint64_t hovered_entity_id_{0};
    std::chrono::steady_clock::time_point last_selection_time_;

    // Input state tracking
    std::unordered_set<int> pressed_keys_;
    std::unordered_set<int> pressed_mouse_buttons_;
    double mouse_x_{0.0};
    double mouse_y_{0.0};
    double mouse_delta_x_{0.0};
    double mouse_delta_y_{0.0};
    double scroll_delta_{0.0};

    // Entity tracking
    std::unordered_map<uint64_t, std::shared_ptr<IEntity>> interactive_entities_;
    mutable std::shared_mutex entities_mutex_;

    // Callbacks
    EntitySelectedCallback entity_selected_callback_;
    EntityHoveredCallback entity_hovered_callback_;
    KeyPressedCallback key_pressed_callback_;
    MouseClickedCallback mouse_clicked_callback_;

    // Helper methods
    void process_keyboard_input();
    void process_mouse_input();
    void process_gesture_input();
    void process_voice_input();
    void update_entity_interactions();
    void perform_ray_casting(double x, double y);
    void handle_entity_selection(uint64_t entity_id);
    void handle_entity_hover(uint64_t entity_id);
};

/**
 * @brief Animation System - Futureproofed Animation Engine
 *
 * Handles all entity animations with support for:
 * - Keyframe animation
 * - Skeletal animation
 * - Particle systems
 * - Procedural animation
 * - Physics-based animation
 */
class AnimationSystem : public ISystem {
public:
    static const std::string SYSTEM_NAME;

    AnimationSystem();
    ~AnimationSystem() override = default;

    // ISystem implementation
    const std::string& get_name() const override { return SYSTEM_NAME; }
    int get_priority() const override { return 40; }
    bool is_enabled() const override { return enabled_; }
    void set_enabled(bool enabled) override { enabled_ = enabled; }

    bool initialize() override;
    void update(double delta_time) override;
    void render() override;
    void shutdown() override;

    void on_entity_created(std::shared_ptr<IEntity> entity) override;
    void on_entity_destroyed(uint64_t entity_id) override;
    void on_component_added(uint64_t entity_id, const std::string& component_type) override;
    void on_component_removed(uint64_t entity_id, const std::string& component_type) override;

    std::vector<std::shared_ptr<IEntity>> query_entities(
        const std::function<bool(const IEntity&)>& predicate) const override;

    double get_average_update_time() const override { return avg_update_time_; }
    double get_average_render_time() const override { return avg_render_time_; }
    size_t get_processed_entities_count() const override { return processed_entities_; }

    // Animation system specific functionality
    void play_animation(uint64_t entity_id, const std::string& animation_name);
    void stop_animation(uint64_t entity_id, const std::string& animation_name);
    void pause_animation(uint64_t entity_id, const std::string& animation_name);
    void resume_animation(uint64_t entity_id, const std::string& animation_name);

    void set_animation_speed(uint64_t entity_id, float speed);
    void set_animation_loop(uint64_t entity_id, bool loop);
    void set_animation_blend(uint64_t entity_id, const std::string& from_anim,
                           const std::string& to_anim, float blend_time);

    // Animation state queries
    bool is_animation_playing(uint64_t entity_id, const std::string& animation_name) const;
    float get_animation_time(uint64_t entity_id, const std::string& animation_name) const;
    float get_animation_duration(uint64_t entity_id, const std::string& animation_name) const;

    // Particle system management
    void create_particle_system(uint64_t entity_id, const std::string& system_name);
    void destroy_particle_system(uint64_t entity_id, const std::string& system_name);
    void set_particle_emission_rate(uint64_t entity_id, float rate);
    void set_particle_lifetime(uint64_t entity_id, float lifetime);

private:
    bool enabled_{true};
    double avg_update_time_{0.0};
    double avg_render_time_{0.0};
    size_t processed_entities_{0};

    // Animation state tracking
    struct AnimationState {
        std::string name;
        float time{0.0f};
        float speed{1.0f};
        float duration{0.0f};
        bool playing{false};
        bool looping{false};
        bool paused{false};
    };

    struct EntityAnimations {
        std::unordered_map<std::string, AnimationState> animations;
        std::string current_animation;
        std::string next_animation;
        float blend_time{0.0f};
        float blend_progress{0.0f};
    };

    std::unordered_map<uint64_t, EntityAnimations> entity_animations_;
    mutable std::shared_mutex animations_mutex_;

    // Entity tracking
    std::unordered_map<uint64_t, std::shared_ptr<IEntity>> animated_entities_;
    mutable std::shared_mutex entities_mutex_;

    // Helper methods
    void update_animations(double delta_time);
    void update_particle_systems(double delta_time);
    void blend_animations(uint64_t entity_id, double delta_time);
    void evaluate_animation(uint64_t entity_id, const std::string& animation_name, float time);
    void apply_animation_result(uint64_t entity_id, const AnimationResult& result);
};

/**
 * @brief System Factory - Futureproofed System Creation
 */
class SystemFactory {
public:
    static std::shared_ptr<ISystem> create_system(const std::string& system_type);

    static std::shared_ptr<RenderingSystem> create_rendering_system();
    static std::shared_ptr<PhysicsSystem> create_physics_system();
    static std::shared_ptr<InteractionSystem> create_interaction_system();
    static std::shared_ptr<AnimationSystem> create_animation_system();

    // Plugin-based system creation
    static std::shared_ptr<ISystem> create_plugin_system(const std::string& plugin_name,
                                                        const std::string& system_name);

private:
    static std::unordered_map<std::string, std::function<std::shared_ptr<ISystem>()>> system_creators_;
};

/**
 * @brief System Manager - Coordinates All Systems
 */
class SystemManager {
public:
    static SystemManager& instance();

    bool register_system(std::shared_ptr<ISystem> system);
    bool unregister_system(const std::string& name);
    std::shared_ptr<ISystem> get_system(const std::string& name) const;

    void initialize_all_systems();
    void update_all_systems(double delta_time);
    void render_all_systems();
    void shutdown_all_systems();

    std::vector<std::shared_ptr<ISystem>> get_systems_by_priority() const;
    std::vector<std::shared_ptr<ISystem>> get_enabled_systems() const;

    // System dependencies and ordering
    void set_system_dependency(const std::string& system, const std::string& depends_on);
    bool check_system_dependencies() const;

    // Performance monitoring
    double get_total_update_time() const;
    double get_total_render_time() const;
    size_t get_total_processed_entities() const;

private:
    SystemManager() = default;
    ~SystemManager();

    std::unordered_map<std::string, std::shared_ptr<ISystem>> systems_;
    std::unordered_map<std::string, std::vector<std::string>> dependencies_;
    mutable std::shared_mutex systems_mutex_;

    bool systems_initialized_{false};
};

} // namespace architecture
} // namespace gui3d
} // namespace hsml
