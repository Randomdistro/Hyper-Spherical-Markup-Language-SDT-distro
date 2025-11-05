/** @file SystemFactory.h
 * @brief System Factory Implementation - Dynamic System Creation
 */

#pragma once

#include "src/presentation/gui/3d/architecture/SystemImplementations.h"
#include <memory>
#include <string>
#include <unordered_map>
#include <functional>
#include <mutex>

namespace hsml {
namespace gui3d {
namespace architecture {

/**
 * @brief System Factory - Dynamic System Creation
 */
class SystemFactory {
public:
    using SystemCreator = std::function<std::shared_ptr<ISystem>()>;

    static std::shared_ptr<ISystem> create_system(const std::string& system_type);

    // Concrete system creation methods
    static std::shared_ptr<RenderingSystem> create_rendering_system();
    static std::shared_ptr<PhysicsSystem> create_physics_system();
    static std::shared_ptr<InteractionSystem> create_interaction_system();
    static std::shared_ptr<AnimationSystem> create_animation_system();

    // Plugin system creation
    static std::shared_ptr<ISystem> create_plugin_system(
        const std::string& plugin_name,
        const std::string& system_name
    );

    // Factory registration for extensibility
    static bool register_system_creator(const std::string& type, SystemCreator creator);
    static bool unregister_system_creator(const std::string& type);

    // Factory introspection
    static std::vector<std::string> get_available_system_types();
    static bool is_system_type_available(const std::string& type);

    // System configuration
    static std::shared_ptr<ISystem> create_system_with_config(
        const std::string& system_type,
        const std::unordered_map<std::string, ConfigValue>& config
    );

private:
    static std::unordered_map<std::string, SystemCreator> system_creators_;
    static std::mutex creators_mutex_;

    // Built-in system creators
    static void register_builtin_creators();
    static bool builtin_creators_registered_;
};

/**
 * @brief System Configuration Helper
 */
class SystemConfigurator {
public:
    static void configure_rendering_system(
        std::shared_ptr<RenderingSystem> system,
        const std::unordered_map<std::string, ConfigValue>& config
    );

    static void configure_physics_system(
        std::shared_ptr<PhysicsSystem> system,
        const std::unordered_map<std::string, ConfigValue>& config
    );

    static void configure_interaction_system(
        std::shared_ptr<InteractionSystem> system,
        const std::unordered_map<std::string, ConfigValue>& config
    );

    static void configure_animation_system(
        std::shared_ptr<AnimationSystem> system,
        const std::unordered_map<std::string, ConfigValue>& config
    );

private:
    static bool get_bool_config(
        const std::unordered_map<std::string, ConfigValue>& config,
        const std::string& key,
        bool default_value = false
    );

    static int get_int_config(
        const std::unordered_map<std::string, ConfigValue>& config,
        const std::string& key,
        int default_value = 0
    );

    static double get_double_config(
        const std::unordered_map<std::string, ConfigValue>& config,
        const std::string& key,
        double default_value = 0.0
    );

    static std::string get_string_config(
        const std::unordered_map<std::string, ConfigValue>& config,
        const std::string& key,
        const std::string& default_value = ""
    );
};

/**
 * @brief System Dependency Resolver
 */
class SystemDependencyResolver {
public:
    static std::vector<std::shared_ptr<ISystem>> resolve_dependencies(
        const std::vector<std::string>& system_types
    );

    static std::vector<std::string> get_system_dependencies(const std::string& system_type);
    static bool has_circular_dependency(const std::vector<std::string>& system_types);

    static std::vector<std::string> topological_sort(const std::vector<std::string>& system_types);

private:
    static const std::unordered_map<std::string, std::vector<std::string>> system_dependencies_;
};

/**
 * @brief System Pool - Object Pool for System Reuse
 */
class SystemPool {
public:
    static SystemPool& instance();

    std::shared_ptr<ISystem> acquire(const std::string& system_type);
    void release(std::shared_ptr<ISystem> system);

    size_t get_pool_size(const std::string& system_type) const;
    void set_pool_limit(const std::string& system_type, size_t limit);

    void clear();
    void optimize();

private:
    SystemPool() = default;
    ~SystemPool() = default;

    struct PooledSystem {
        std::shared_ptr<ISystem> system;
        std::chrono::steady_clock::time_point last_used;
        bool in_use;
    };

    std::unordered_map<std::string, std::vector<PooledSystem>> pools_;
    std::unordered_map<std::string, size_t> pool_limits_;
    mutable std::shared_mutex pools_mutex_;
};

} // namespace architecture
} // namespace gui3d
} // namespace hsml
