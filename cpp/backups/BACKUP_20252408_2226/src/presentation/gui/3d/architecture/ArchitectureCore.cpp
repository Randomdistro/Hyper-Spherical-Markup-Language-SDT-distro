/** @file ArchitectureCore.h
 * @brief Perfect 3D GUI Architecture Core - Futureproofed Design
 *
 * Phase 6: Perfect Architecture Implementation
 * Core architectural framework for the 3D GUI system with futureproofing.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <future>
#include <chrono>

namespace hsml {
namespace gui3d {
namespace architecture {

// Forward declarations
class IComponent;
class ISystem;
class IEntity;
class IPlugin;
class IEvent;
class IMessageBus;
class IResourceManager;
class IConfigurationManager;

// Type definitions for the architecture
enum class ConnectionType {
    FUNCTION_CALL,
    INHERITANCE,
    COMPOSITION,
    DEPENDENCY,
    DATA_FLOW,
    EVENT_FLOW,
    ASYNC_CALL,
    TEMPLATE_USAGE,
    INTERFACE_IMPLEMENTATION,
    FRIEND_RELATIONSHIP
};

// Simple vector and coordinate types for the architecture
struct Vector3 {
    float x, y, z;
    Vector3(float x = 0.0f, float y = 0.0f, float z = 0.0f) : x(x), y(y), z(z) {}
};

struct SphericalCoords {
    double radius;
    double theta;  // azimuthal angle
    double phi;    // polar angle
    SphericalCoords(double r = 0.0, double t = 0.0, double p = 0.0) : radius(r), theta(t), phi(p) {}
};

/**
 * @brief Architecture constants and configuration
 */
namespace constants {
    const size_t MAX_ENTITIES = 10000;
    const size_t MAX_CONNECTIONS = 50000;
    const size_t MAX_COMPONENTS_PER_ENTITY = 64;
    const double TARGET_FRAME_RATE = 60.0;
    const double MAX_FRAME_TIME_MS = 16.67;
    const size_t DEFAULT_THREAD_POOL_SIZE = 4;
    const size_t MAX_MEMORY_USAGE_MB = 2048;
}

/**
 * @brief Entity Component System (ECS) Architecture
 *
 * Futureproofed ECS design that allows for:
 * - Dynamic component addition/removal
 * - Plugin-based system extensions
 * - Efficient memory management
 * - Thread-safe operations
 * - Serialization/deserialization
 */
class IEntity {
public:
    virtual ~IEntity() = default;

    virtual uint64_t get_id() const = 0;
    virtual const std::string& get_name() const = 0;
    virtual void set_name(const std::string& name) = 0;

    // Component management
    virtual bool has_component(const std::string& type) const = 0;
    virtual IComponent* get_component(const std::string& type) = 0;
    virtual const IComponent* get_component(const std::string& type) const = 0;

    virtual bool add_component(std::unique_ptr<IComponent> component) = 0;
    virtual bool remove_component(const std::string& type) = 0;
    virtual void clear_components() = 0;

    virtual std::vector<std::string> get_component_types() const = 0;
    virtual size_t get_component_count() const = 0;

    // Entity state
    virtual bool is_active() const = 0;
    virtual void set_active(bool active) = 0;

    virtual bool is_visible() const = 0;
    virtual void set_visible(bool visible) = 0;

    // Lifecycle
    virtual void initialize() = 0;
    virtual void update(double delta_time) = 0;
    virtual void render() = 0;
    virtual void shutdown() = 0;

    // Serialization
    virtual std::string serialize() const = 0;
    virtual bool deserialize(const std::string& data) = 0;

    // Event handling
    virtual void send_event(std::shared_ptr<IEvent> event) = 0;
    virtual void receive_event(std::shared_ptr<IEvent> event) = 0;
};

/**
 * @brief Component Interface - Futureproofed Component System
 *
 * Components are pure data with minimal logic, following ECS principles.
 * This design allows for maximum flexibility and extensibility.
 */
class IComponent {
public:
    virtual ~IComponent() = default;

    virtual const std::string& get_type() const = 0;
    virtual uint64_t get_entity_id() const = 0;
    virtual void set_entity_id(uint64_t id) = 0;

    virtual bool is_enabled() const = 0;
    virtual void set_enabled(bool enabled) = 0;

    // Lifecycle
    virtual void initialize() = 0;
    virtual void update(double delta_time) = 0;
    virtual void render() = 0;
    virtual void shutdown() = 0;

    // Serialization
    virtual std::string serialize() const = 0;
    virtual bool deserialize(const std::string& data) = 0;

    // Memory management
    virtual size_t get_memory_usage() const = 0;
    virtual void optimize_memory() = 0;
};

/**
 * @brief System Interface - Plugin-based Architecture
 *
 * Systems contain the logic and operate on components.
 * Plugin-based design allows for runtime extension and customization.
 */
class ISystem {
public:
    virtual ~ISystem() = default;

    virtual const std::string& get_name() const = 0;
    virtual int get_priority() const = 0; // Lower numbers = higher priority
    virtual bool is_enabled() const = 0;
    virtual void set_enabled(bool enabled) = 0;

    // System lifecycle
    virtual bool initialize() = 0;
    virtual void update(double delta_time) = 0;
    virtual void render() = 0;
    virtual void shutdown() = 0;

    // Entity management
    virtual void on_entity_created(std::shared_ptr<IEntity> entity) = 0;
    virtual void on_entity_destroyed(uint64_t entity_id) = 0;
    virtual void on_component_added(uint64_t entity_id, const std::string& component_type) = 0;
    virtual void on_component_removed(uint64_t entity_id, const std::string& component_type) = 0;

    // Query interface
    virtual std::vector<std::shared_ptr<IEntity>> query_entities(
        const std::function<bool(const IEntity&)>& predicate) const = 0;

    // Performance monitoring
    virtual double get_average_update_time() const = 0;
    virtual double get_average_render_time() const = 0;
    virtual size_t get_processed_entities_count() const = 0;
};

/**
 * @brief Plugin Interface - Futureproofing through Extensibility
 *
 * Plugin system allows for runtime loading of new functionality,
 * enabling the system to grow and adapt without core changes.
 */
class IPlugin {
public:
    virtual ~IPlugin() = default;

    virtual const std::string& get_name() const = 0;
    virtual const std::string& get_version() const = 0;
    virtual const std::string& get_description() const = 0;

    virtual bool is_compatible() const = 0;
    virtual std::vector<std::string> get_dependencies() const = 0;

    // Plugin lifecycle
    virtual bool load() = 0;
    virtual bool unload() = 0;
    virtual bool is_loaded() const = 0;

    // System integration
    virtual std::vector<std::shared_ptr<ISystem>> get_systems() const = 0;
    virtual std::vector<std::string> get_provided_components() const = 0;

    // Resource management
    virtual size_t get_memory_usage() const = 0;
    virtual void optimize_resources() = 0;
};

/**
 * @brief Message Bus - Decoupled Communication
 *
 * Event-driven architecture for loose coupling between systems and components.
 * Supports both synchronous and asynchronous message delivery.
 */
class IMessageBus {
public:
    virtual ~IMessageBus() = default;

    // Event publishing
    virtual void publish(std::shared_ptr<IEvent> event) = 0;
    virtual void publish_async(std::shared_ptr<IEvent> event) = 0;

    // Event subscription
    using EventHandler = std::function<void(std::shared_ptr<IEvent>)>;
    virtual uint64_t subscribe(const std::string& event_type, EventHandler handler) = 0;
    virtual bool unsubscribe(uint64_t subscription_id) = 0;

    // Event filtering
    virtual std::vector<std::shared_ptr<IEvent>> get_events(
        const std::string& event_type,
        std::chrono::milliseconds timeout = std::chrono::milliseconds(0)) = 0;

    // Performance monitoring
    virtual size_t get_published_event_count() const = 0;
    virtual size_t get_subscriber_count() const = 0;
    virtual double get_average_delivery_time() const = 0;
};

/**
 * @brief Event Base Class
 */
class IEvent {
public:
    virtual ~IEvent() = default;

    virtual const std::string& get_type() const = 0;
    virtual uint64_t get_source_id() const = 0;
    virtual std::chrono::steady_clock::time_point get_timestamp() const = 0;
    virtual uint64_t get_sequence_number() const = 0;

    virtual std::string serialize() const = 0;
    virtual bool deserialize(const std::string& data) = 0;
};

/**
 * @brief Resource Manager - Futureproofed Asset Management
 */
class IResourceManager {
public:
    virtual ~IResourceManager() = default;

    // Resource loading
    virtual bool load_resource(const std::string& path, const std::string& type) = 0;
    virtual bool unload_resource(const std::string& path) = 0;
    virtual bool is_resource_loaded(const std::string& path) const = 0;

    // Resource access
    virtual void* get_resource(const std::string& path) const = 0;
    virtual size_t get_resource_size(const std::string& path) const = 0;
    virtual std::string get_resource_type(const std::string& path) const = 0;

    // Resource streaming (for large assets)
    virtual bool start_streaming(const std::string& path) = 0;
    virtual bool stop_streaming(const std::string& path) = 0;
    virtual float get_streaming_progress(const std::string& path) const = 0;

    // Memory management
    virtual size_t get_total_memory_usage() const = 0;
    virtual size_t get_cache_size() const = 0;
    virtual void set_memory_limit(size_t bytes) = 0;
    virtual void optimize_memory() = 0;
};

/**
 * @brief Configuration Manager - Runtime Configuration
 */
class IConfigurationManager {
public:
    virtual ~IConfigurationManager() = default;

    // Configuration loading/saving
    virtual bool load_config(const std::string& path) = 0;
    virtual bool save_config(const std::string& path) const = 0;

    // Value access
    virtual bool get_bool(const std::string& key, bool default_value = false) const = 0;
    virtual int get_int(const std::string& key, int default_value = 0) const = 0;
    virtual double get_double(const std::string& key, double default_value = 0.0) const = 0;
    virtual std::string get_string(const std::string& key, const std::string& default_value = "") const = 0;

    virtual bool set_bool(const std::string& key, bool value) = 0;
    virtual bool set_int(const std::string& key, int value) = 0;
    virtual bool set_double(const std::string& key, double value) = 0;
    virtual bool set_string(const std::string& key, const std::string& value) = 0;

    // Section management
    virtual std::vector<std::string> get_sections() const = 0;
    virtual std::vector<std::string> get_keys(const std::string& section) const = 0;

    // Runtime configuration updates
    virtual void reload_config() = 0;
    virtual bool is_config_modified() const = 0;
};

/**
 * @brief Architecture Statistics - Performance Monitoring
 */
struct ArchitectureStatistics {
    // Entity statistics
    size_t total_entities{0};
    size_t active_entities{0};
    size_t visible_entities{0};

    // Component statistics
    size_t total_components{0};
    size_t enabled_components{0};
    std::unordered_map<std::string, size_t> components_by_type;

    // System statistics
    size_t total_systems{0};
    size_t enabled_systems{0};
    double average_system_update_time{0.0};
    double average_system_render_time{0.0};

    // Memory statistics
    size_t total_memory_usage{0};
    size_t resource_memory_usage{0};
    size_t entity_memory_usage{0};
    size_t component_memory_usage{0};

    // Performance statistics
    double frame_rate{0.0};
    double average_frame_time{0.0};
    double max_frame_time{0.0};
    double min_frame_time{0.0};

    // Event statistics
    size_t events_published{0};
    size_t events_delivered{0};
    size_t events_pending{0};

    // Threading statistics
    size_t active_threads{0};
    size_t pending_tasks{0};
    double average_task_completion_time{0.0};
};

/**
 * @brief Architecture Manager - Central Coordination
 */
class ArchitectureManager {
public:
    static ArchitectureManager& instance();

    // Initialization
    bool initialize();
    void shutdown();
    bool is_initialized() const { return initialized_; }

    // Entity management
    std::shared_ptr<IEntity> create_entity(const std::string& name = "");
    bool destroy_entity(uint64_t entity_id);
    std::shared_ptr<IEntity> get_entity(uint64_t entity_id) const;
    std::vector<std::shared_ptr<IEntity>> get_all_entities() const;
    std::vector<std::shared_ptr<IEntity>> query_entities(
        const std::function<bool(const IEntity&)>& predicate) const;

    // System management
    bool register_system(std::shared_ptr<ISystem> system);
    bool unregister_system(const std::string& name);
    std::shared_ptr<ISystem> get_system(const std::string& name) const;

    // Plugin management
    bool load_plugin(const std::string& path);
    bool unload_plugin(const std::string& name);
    std::shared_ptr<IPlugin> get_plugin(const std::string& name) const;
    std::vector<std::shared_ptr<IPlugin>> get_loaded_plugins() const;

    // Core services
    std::shared_ptr<IMessageBus> get_message_bus() const { return message_bus_; }
    std::shared_ptr<IResourceManager> get_resource_manager() const { return resource_manager_; }
    std::shared_ptr<IConfigurationManager> get_config_manager() const { return config_manager_; }

    // Main loop
    void update(double delta_time);
    void render();

    // Statistics and monitoring
    const ArchitectureStatistics& get_statistics() const { return statistics_; }
    void reset_statistics();

    // Configuration
    void set_target_frame_rate(double fps);
    void set_memory_limit(size_t mb);
    void set_thread_pool_size(size_t count);

private:
    ArchitectureManager();
    ~ArchitectureManager();

    bool initialized_{false};
    ArchitectureStatistics statistics_;
    std::chrono::steady_clock::time_point last_update_time_;

    // Core systems
    std::shared_ptr<IMessageBus> message_bus_;
    std::shared_ptr<IResourceManager> resource_manager_;
    std::shared_ptr<IConfigurationManager> config_manager_;

    // Entity and system management
    std::unordered_map<uint64_t, std::shared_ptr<IEntity>> entities_;
    std::unordered_map<std::string, std::shared_ptr<ISystem>> systems_;
    std::unordered_map<std::string, std::shared_ptr<IPlugin>> plugins_;

    // Thread safety
    mutable std::shared_mutex entities_mutex_;
    mutable std::shared_mutex systems_mutex_;
    mutable std::shared_mutex plugins_mutex_;

    // Entity ID generation
    std::atomic<uint64_t> next_entity_id_{1};

    // Helper methods
    void update_statistics(double delta_time);
    void process_events();
    void optimize_memory();
    void validate_system_dependencies();
};

} // namespace architecture
} // namespace gui3d
} // namespace hsml
