/** @file Perfect3DGUI.h
 * @brief Perfect 3D GUI Integration Layer
 *
 * Phase 6: Perfect Architecture - Main integration layer
 * that brings together all components into a cohesive system.
 */

#pragma once

#include "src/presentation/gui/3d/architecture/ArchitectureCore.h"
#include "src/presentation/gui/3d/architecture/EntitySystem.h"
#include "src/presentation/gui/3d/architecture/SystemImplementations.h"
#include "../../../application/services/SphericalSceneService.h"
#include "../../../core/spherical_coords.h"
#include "../../../core/vector3.h"
#include <memory>
#include <string>
#include <vector>
#include <functional>

namespace hsml {
namespace gui3d {

/**
 * @brief Perfect 3D GUI - Futureproofed Integration
 *
 * This is the main integration layer that brings together:
 * - Futureproofed ECS Architecture
 * - Plugin-based System Management
 * - Advanced Entity Component System
 * - Real-time Performance Monitoring
 * - Comprehensive Memory Management
 * - Thread-safe Operations
 */
class Perfect3DGUI {
public:
    /**
     * @brief Create a Perfect 3D GUI instance with dependency injection
     *
     * @param scene_service The application service for scene management
     * @param message_bus Message bus for event communication
     * @param resource_manager Resource manager for asset management
     * @param config_manager Configuration manager for runtime settings
     * @param width Initial viewport width
     * @param height Initial viewport height
     * @return Perfect 3D GUI instance or nullptr on failure
     */
    static std::unique_ptr<Perfect3DGUI> create(
        std::shared_ptr<application::SphericalSceneService> scene_service,
        std::shared_ptr<architecture::IMessageBus> message_bus = nullptr,
        std::shared_ptr<architecture::IResourceManager> resource_manager = nullptr,
        std::shared_ptr<architecture::IConfigurationManager> config_manager = nullptr,
        int width = 1920,
        int height = 1080
    );

    /**
     * @brief Initialize the Perfect 3D GUI System (legacy static interface)
     *
     * @param scene_service The application service for scene management
     * @param width Initial viewport width
     * @param height Initial viewport height
     * @return true if initialization successful
     */
    static bool initialize(
        std::shared_ptr<application::SphericalSceneService> scene_service,
        int width = 1920,
        int height = 1080
    );

    /**
     * @brief Shutdown the 3D GUI system
     */
    static void shutdown();

    /**
     * @brief Check if the system is initialized
     */
    static bool is_initialized();

    /**
     * @brief Get the architecture manager
     */
    static architecture::ArchitectureManager& get_architecture_manager();

    // Instance methods (when using dependency injection)
    bool initialize_instance();
    void shutdown_instance();
    bool is_initialized_instance() const { return initialized_; }

    architecture::ArchitectureManager& get_architecture_manager_instance() {
        return *architecture_manager_;
    }

    /**
     * @brief Create a floating code entity
     */
    static std::shared_ptr<architecture::Entity> create_code_entity(
        const std::string& filename,
        const std::string& content,
        const SphericalCoords& position
    );

    /**
     * @brief Create a connection between entities
     */
    static std::shared_ptr<architecture::Entity> create_connection_entity(
        std::shared_ptr<architecture::Entity> from,
        std::shared_ptr<architecture::Entity> to,
        architecture::ConnectionType type
    );

    // Entity lookup functions
    static std::shared_ptr<architecture::Entity> find_entity_by_filename(const std::string& filename);
    static std::shared_ptr<architecture::Entity> find_entity_by_id(uint64_t id);

    /**
     * @brief Load a project into the 3D space
     */
    static bool load_project(const std::string& project_path);

    /**
     * @brief Save current project state
     */
    static bool save_project(const std::string& project_path = "");

    /**
     * @brief Get system statistics
     */
    static architecture::ArchitectureStatistics get_statistics();

    /**
     * @brief Enable/disable performance monitoring
     */
    static void set_performance_monitoring(bool enabled);

    /**
     * @brief Set target frame rate
     */
    static void set_target_frame_rate(double fps);

    /**
     * @brief Main update loop
     */
    static void update(double delta_time);

    /**
     * @brief Main render loop
     */
    static void render();

    /**
     * @brief Handle input events
     */
    static void handle_input(int type, int key, int action, int mods);
    static void handle_mouse(double x, double y, int button, int action);
    static void handle_scroll(double x_offset, double y_offset);

    /**
     * @brief Camera controls
     */
    static void set_camera_position(const SphericalCoords& position);
    static const SphericalCoords& get_camera_position();
    static void orbit_around_entity(uint64_t entity_id, double duration = 2.0);
    static void focus_on_entity(uint64_t entity_id);

    /**
     * @brief Entity management
     */
    static std::shared_ptr<architecture::Entity> get_entity(uint64_t id);
    static std::vector<std::shared_ptr<architecture::Entity>> get_all_entities();
    static std::vector<std::shared_ptr<architecture::Entity>> query_entities(
        const std::function<bool(const architecture::Entity&)>& predicate
    );

    /**
     * @brief System management
     */
    static std::shared_ptr<architecture::ISystem> get_system(const std::string& name);
    static bool enable_system(const std::string& name);
    static bool disable_system(const std::string& name);

    /**
     * @brief Plugin management
     */
    static bool load_plugin(const std::string& path);
    static bool unload_plugin(const std::string& name);
    static std::vector<std::string> get_loaded_plugins();

    /**
     * @brief Configuration management
     */
    static bool load_config(const std::string& path);
    static bool save_config(const std::string& path);

    /**
     * @brief Resource management
     */
    static bool load_resource(const std::string& path, const std::string& type);
    static void* get_resource(const std::string& path);

    /**
     * @brief Event system
     */
    using EntitySelectedCallback = std::function<void(std::shared_ptr<architecture::Entity>)>;
    using ConnectionCreatedCallback = std::function<void(std::shared_ptr<architecture::Entity>)>;
    using SystemEventCallback = std::function<void(const std::string&, const std::string&)>;

    static void set_entity_selected_callback(EntitySelectedCallback callback);
    static void set_connection_created_callback(ConnectionCreatedCallback callback);
    static void set_system_event_callback(SystemEventCallback callback);

private:
    // Private constructor for dependency injection
    Perfect3DGUI(
        std::shared_ptr<application::SphericalSceneService> scene_service,
        std::shared_ptr<architecture::IMessageBus> message_bus,
        std::shared_ptr<architecture::IResourceManager> resource_manager,
        std::shared_ptr<architecture::IConfigurationManager> config_manager,
        int width = 1920,
        int height = 1080
    );
    ~Perfect3DGUI();

    // Legacy singleton support
    static Perfect3DGUI* instance_;
    static std::mutex instance_mutex_;

    // Core systems
    std::shared_ptr<application::SphericalSceneService> scene_service_;
    std::shared_ptr<architecture::IMessageBus> message_bus_;
    std::shared_ptr<architecture::IResourceManager> resource_manager_;
    std::shared_ptr<architecture::IConfigurationManager> config_manager_;
    std::unique_ptr<architecture::ArchitectureManager> architecture_manager_;

    // State
    bool initialized_;
    std::string current_project_;
    int viewport_width_;
    int viewport_height_;

    // Callbacks
    EntitySelectedCallback entity_selected_callback_;
    ConnectionCreatedCallback connection_created_callback_;
    SystemEventCallback system_event_callback_;

    // Helper methods
    bool initialize_architecture();
    bool initialize_default_systems();
    bool initialize_default_plugins();
    void setup_event_handling();
    void configure_default_settings();

    // Project management
    bool parse_project_structure(const std::string& project_path);
    void arrange_project_entities();
    void create_project_connections();

    // Entity creation helpers
    std::shared_ptr<architecture::Entity> create_transform_component(
        const SphericalCoords& position
    );
    std::shared_ptr<architecture::Entity> create_visual_component(
        const Vector3& color = Vector3(1.0f, 1.0f, 1.0f)
    );
    std::shared_ptr<architecture::Entity> create_code_component(
        const std::string& filename,
        const std::string& content
    );
    std::shared_ptr<architecture::Entity> create_connection_component(
        uint64_t target_id,
        architecture::ConnectionType type
    );

    // System helpers
    void update_camera_system();
    void update_interaction_system();
    void process_pending_events();

    // Memory and performance helpers
    void optimize_memory_usage();
    void balance_system_load();
    void update_performance_metrics();
};

/**
 * @brief 3D GUI Factory - Easy initialization
 */
class Perfect3DGUIFactory {
public:
    /**
     * @brief Create and initialize a 3D GUI instance
     */
    static std::unique_ptr<Perfect3DGUI> create(
        std::shared_ptr<application::SphericalSceneService> scene_service,
        int width = 1920,
        int height = 1080
    );

    /**
     * @brief Create with project loading
     */
    static std::unique_ptr<Perfect3DGUI> create_with_project(
        const std::string& project_path,
        std::shared_ptr<application::SphericalSceneService> scene_service,
        int width = 1920,
        int height = 1080
    );

    /**
     * @brief Create with custom configuration
     */
    static std::unique_ptr<Perfect3DGUI> create_with_config(
        const std::string& config_path,
        std::shared_ptr<application::SphericalSceneService> scene_service,
        int width = 1920,
        int height = 1080
    );
};

/**
 * @brief 3D GUI Configuration
 */
struct Perfect3DGUIConfig {
    // Display settings
    int window_width{1920};
    int window_height{1080};
    double target_frame_rate{60.0};
    bool fullscreen{false};
    bool vsync{true};
    int antialiasing_samples{4};

    // System settings
    bool enable_physics{true};
    bool enable_animations{true};
    bool enable_particles{true};
    bool enable_post_processing{true};

    // Performance settings
    size_t max_entities{10000};
    size_t max_memory_mb{2048};
    size_t thread_pool_size{4};
    bool enable_performance_monitoring{true};

    // Visual settings
    std::string theme{"dark_coding"};
    Vector3 background_color{0.02f, 0.02f, 0.05f};
    float ambient_light_intensity{0.3f};
    bool depth_of_field_enabled{true};
    float field_of_view{60.0f};

    // Interaction settings
    float mouse_sensitivity{1.0f};
    float keyboard_repeat_rate{30.0f};
    bool gesture_recognition_enabled{false};
    bool haptic_feedback_enabled{false};

    // Plugin settings
    std::vector<std::string> default_plugins{
        "rendering", "physics", "interaction", "animation"
    };
    std::vector<std::string> optional_plugins{
        "networking", "audio", "video", "scripting"
    };

    // Debug settings
    bool enable_debug_overlay{false};
    bool enable_profiling{false};
    bool enable_memory_tracking{true};
    std::string log_level{"info"};
};

/**
 * @brief Global configuration instance
 */
extern Perfect3DGUIConfig global_config;

/**
 * @brief Quick start macros
 */
#define INIT_PERFECT_3D_GUI(scene_service) \
    Perfect3DGUI::initialize(scene_service)

#define CREATE_CODE_ENTITY(filename, content, position) \
    Perfect3DGUI::create_code_entity(filename, content, position)

#define CREATE_CONNECTION(from, to, type) \
    Perfect3DGUI::create_connection_entity(from, to, type)

#define UPDATE_3D_GUI(delta_time) \
    Perfect3DGUI::update(delta_time)

#define RENDER_3D_GUI() \
    Perfect3DGUI::render()

#define SHUTDOWN_3D_GUI() \
    Perfect3DGUI::shutdown()

} // namespace gui3d
} // namespace hsml
