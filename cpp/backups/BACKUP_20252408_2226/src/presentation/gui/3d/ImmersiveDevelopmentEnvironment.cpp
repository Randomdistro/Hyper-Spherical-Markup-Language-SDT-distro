/** @file ImmersiveDevelopmentEnvironment.h
 * @brief Main 3D GUI Environment for Immersive Development
 *
 * Phase 5: Immersive 3D Development Environment
 * The core orchestrator for the 3D development experience.
 */

#pragma once

#include "hsml/presentation/gui/SphericalViewport.h"
#include "hsml/presentation/gui/3d/entities/FloatingCodeEntity.h"
#include "hsml/presentation/gui/3d/connections/ConnectionLine.h"
#include "hsml/application/services/SphericalSceneService.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <chrono>

namespace hsml {
namespace presentation {

// Forward declarations
class CodeExecutionVisualizer;
class NavigationController;
class InteractionManager;
class SpatialAudioSystem;

/**
 * @brief Main 3D Development Environment
 *
 * This is the core orchestrator for the immersive 3D development experience.
 * It manages floating code entities, connection lines, real-time execution
 * visualization, and all user interactions within the 3D space.
 */
class ImmersiveDevelopmentEnvironment {
public:
    /**
     * @brief Initialize the 3D development environment
     *
     * @param scene_service The application service for scene management
     * @param width Initial viewport width
     * @param height Initial viewport height
     */
    ImmersiveDevelopmentEnvironment(
        std::shared_ptr<application::SphericalSceneService> scene_service,
        int width = 1920,
        int height = 1080
    );

    ~ImmersiveDevelopmentEnvironment();

    // Environment lifecycle
    bool initialize();
    void shutdown();
    bool is_initialized() const { return initialized_; }

    // Project loading and management
    bool load_project(const std::string& project_path);
    bool save_project(const std::string& project_path = "");
    void create_new_project(const std::string& project_name);
    const std::string& get_current_project() const { return current_project_; }

    // Entity management
    std::shared_ptr<FloatingCodeEntity> add_code_file(
        const std::string& filepath,
        const SphericalCoords& position
    );

    void remove_code_entity(std::shared_ptr<FloatingCodeEntity> entity);
    void move_code_entity(std::shared_ptr<FloatingCodeEntity> entity, const SphericalCoords& new_position);
    std::shared_ptr<FloatingCodeEntity> find_entity_by_filename(const std::string& filename) const;

    // Connection management
    std::shared_ptr<ConnectionLine> create_connection(
        std::shared_ptr<FloatingCodeEntity> from,
        std::shared_ptr<FloatingCodeEntity> to,
        ConnectionType type,
        const std::string& description = ""
    );

    void remove_connection(std::shared_ptr<ConnectionLine> connection);
    void analyze_and_create_connections();

    // Navigation and camera control
    void set_camera_position(const SphericalCoords& position);
    const SphericalCoords& get_camera_position() const;
    void orbit_around_entity(std::shared_ptr<FloatingCodeEntity> entity, float duration = 2.0f);
    void focus_on_entity(std::shared_ptr<FloatingCodeEntity> entity);
    void reset_camera_to_overview();

    // Interaction modes
    enum class InteractionMode {
        NAVIGATION,     // General navigation and exploration
        EDITING,        // Code editing mode
        DEBUGGING,      // Debug visualization
        PROFILING,      // Performance profiling
        COLLABORATION   // Multi-user collaboration
    };

    void set_interaction_mode(InteractionMode mode);
    InteractionMode get_interaction_mode() const { return interaction_mode_; }

    // Real-time code execution visualization
    void start_code_execution(const std::string& function_name,
                             std::shared_ptr<FloatingCodeEntity> entity);
    void stop_code_execution();
    bool is_code_executing() const;
    void set_execution_speed(float multiplier) { execution_speed_ = multiplier; }
    float get_execution_speed() const { return execution_speed_; }

    // Search and filtering
    void search_code(const std::string& query);
    void filter_by_language(const std::string& language);
    void filter_by_file_type(const std::string& extension);
    void show_only_connected_entities();
    void clear_filters();

    // Visual customization
    void set_environment_theme(const std::string& theme_name);
    void set_background_color(const Vector3& color);
    void set_ambient_light_intensity(float intensity);
    void enable_depth_of_field(bool enabled);
    void set_field_of_view(float fov_degrees);

    // Performance and statistics
    struct EnvironmentStatistics {
        size_t entity_count;
        size_t connection_count;
        size_t active_connections;
        float average_frame_time_ms;
        float memory_usage_mb;
        size_t triangles_rendered;
        size_t draw_calls;
    };

    EnvironmentStatistics get_statistics() const;
    void enable_performance_monitoring(bool enabled);
    bool is_performance_monitoring_enabled() const { return performance_monitoring_; }

    // Main loop
    void update(float delta_time);
    void render();
    bool should_close() const;

    // Event handling
    void handle_mouse_press(int x, int y, int button);
    void handle_mouse_release(int x, int y, int button);
    void handle_mouse_move(int x, int y);
    void handle_mouse_wheel(int delta);
    void handle_key_press(int key);
    void handle_key_release(int key);
    void handle_text_input(const std::string& text);

    // Event callbacks
    using EntitySelectedCallback = std::function<void(std::shared_ptr<FloatingCodeEntity>)>;
    using ConnectionActivatedCallback = std::function<void(std::shared_ptr<ConnectionLine>)>;
    using ProjectLoadedCallback = std::function<void(const std::string&)>;
    using ExecutionStartedCallback = std::function<void(const std::string&, std::shared_ptr<FloatingCodeEntity>)>;

    void set_entity_selected_callback(EntitySelectedCallback callback) {
        entity_selected_callback_ = callback;
    }
    void set_connection_activated_callback(ConnectionActivatedCallback callback) {
        connection_activated_callback_ = callback;
    }
    void set_project_loaded_callback(ProjectLoadedCallback callback) {
        project_loaded_callback_ = callback;
    }
    void set_execution_started_callback(ExecutionStartedCallback callback) {
        execution_started_callback_ = callback;
    }

private:
    // Core systems
    std::shared_ptr<application::SphericalSceneService> scene_service_;
    std::unique_ptr<SphericalViewport> viewport_;
    std::unique_ptr<CodeExecutionVisualizer> execution_visualizer_;
    std::unique_ptr<NavigationController> navigation_controller_;
    std::unique_ptr<InteractionManager> interaction_manager_;
    std::unique_ptr<SpatialAudioSystem> audio_system_;

    // Environment state
    bool initialized_{false};
    std::string current_project_;
    InteractionMode interaction_mode_{InteractionMode::NAVIGATION};
    float execution_speed_{1.0f};
    bool performance_monitoring_{false};

    // Entities and connections
    std::vector<std::shared_ptr<FloatingCodeEntity>> code_entities_;
    std::unordered_map<std::string, std::shared_ptr<FloatingCodeEntity>> entity_map_;

    // Visual environment
    Vector3 background_color_{0.02f, 0.02f, 0.05f};
    float ambient_light_intensity_{0.3f};
    bool depth_of_field_enabled_{true};
    float field_of_view_{60.0f};
    std::string current_theme_{"dark_coding"};

    // Performance tracking
    std::chrono::steady_clock::time_point last_frame_time_;
    std::vector<float> frame_times_;
    size_t max_frame_time_history_{300}; // 5 seconds at 60fps
    mutable EnvironmentStatistics cached_stats_;

    // Search and filtering state
    std::string current_search_query_;
    std::vector<std::string> active_filters_;
    bool show_only_connected_{false};

    // Callbacks
    EntitySelectedCallback entity_selected_callback_;
    ConnectionActivatedCallback connection_activated_callback_;
    ProjectLoadedCallback project_loaded_callback_;
    ExecutionStartedCallback execution_started_callback_;

    // Helper methods
    void initialize_visual_environment();
    void load_theme_settings(const std::string& theme_name);
    void setup_default_layout();
    void update_statistics();
    void handle_entity_selection(std::shared_ptr<FloatingCodeEntity> entity);
    void handle_connection_activation(std::shared_ptr<ConnectionLine> connection);
    void apply_filters();
    void optimize_layout();

    // Project management helpers
    bool parse_project_structure(const std::string& project_path);
    void arrange_entities_in_space();
    void create_initial_connections();
    void save_project_state(const std::string& filepath);
    bool load_project_state(const std::string& filepath);

    // Input handling helpers
    void process_navigation_input();
    void process_editing_input();
    void process_debugging_input();
    void process_search_input(const std::string& text);

    // Rendering helpers
    void render_entities();
    void render_connections();
    void render_execution_visualization();
    void render_ui_overlays();
    void render_performance_overlay();

    // Audio helpers
    void play_entity_select_sound();
    void play_connection_activate_sound();
    void play_navigation_sound();
    void update_spatial_audio();
};

/**
 * @brief Factory for creating development environments
 */
class ImmersiveEnvironmentFactory {
public:
    static std::unique_ptr<ImmersiveDevelopmentEnvironment> create_environment(
        std::shared_ptr<application::SphericalSceneService> scene_service,
        int width = 1920,
        int height = 1080
    );

    static std::unique_ptr<ImmersiveDevelopmentEnvironment> create_from_project(
        const std::string& project_path,
        std::shared_ptr<application::SphericalSceneService> scene_service,
        int width = 1920,
        int height = 1080
    );
};

/**
 * @brief Configuration for the development environment
 */
struct DevelopmentEnvironmentConfig {
    // Display settings
    int window_width{1920};
    int window_height{1080};
    bool fullscreen{false};
    bool vsync{true};
    int antialiasing_samples{4};

    // Visual settings
    std::string default_theme{"dark_coding"};
    Vector3 background_color{0.02f, 0.02f, 0.05f};
    float ambient_light_intensity{0.3f};
    bool depth_of_field_enabled{true};
    float field_of_view{60.0f};

    // Entity settings
    float default_entity_radius{2.0f};
    Vector3 default_entity_color{0.2f, 0.6f, 1.0f};
    bool syntax_highlighting_enabled{true};

    // Connection settings
    float default_connection_thickness{2.0f};
    Vector3 default_connection_color{0.0f, 1.0f, 0.0f};
    bool animate_connections{true};

    // Performance settings
    bool performance_monitoring_enabled{true};
    size_t max_entities{1000};
    size_t max_connections{5000};
    float target_frame_rate{60.0f};

    // Audio settings
    bool spatial_audio_enabled{true};
    float master_volume{0.7f};
    bool sound_effects_enabled{true};

    // Interaction settings
    bool mouse_sensitivity_enabled{true};
    float mouse_sensitivity{1.0f};
    bool gesture_recognition_enabled{false};
    bool haptic_feedback_enabled{false};
};

/**
 * @brief Environment presets for different use cases
 */
class EnvironmentPresets {
public:
    static DevelopmentEnvironmentConfig create_coding_preset();
    static DevelopmentEnvironmentConfig create_debugging_preset();
    static DevelopmentEnvironmentConfig create_collaboration_preset();
    static DevelopmentEnvironmentConfig create_presentation_preset();
    static DevelopmentEnvironmentConfig create_minimal_preset();
    static DevelopmentEnvironmentConfig create_maximum_preset();
};

} // namespace presentation
} // namespace hsml
