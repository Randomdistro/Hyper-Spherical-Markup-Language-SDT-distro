/** @file ArchitectureManager.cpp
 * @brief Futureproofed Architecture Manager Implementation
 */

#include "src/presentation/gui/3d/architecture/ArchitectureCore.h"
#include "src/presentation/gui/3d/architecture/SystemImplementations.h"
#include <algorithm>
#include <sstream>

namespace hsml {
namespace gui3d {
namespace architecture {

ArchitectureManager::ArchitectureManager()
    : initialized_(false)
    , message_bus_(std::make_shared<MessageBus>())
    , resource_manager_(std::make_shared<ResourceManager>())
    , config_manager_(std::make_shared<ConfigurationManager>())
{
}

ArchitectureManager::ArchitectureManager(
    std::shared_ptr<IMessageBus> message_bus,
    std::shared_ptr<IResourceManager> resource_manager,
    std::shared_ptr<IConfigurationManager> config_manager)
    : initialized_(false)
    , message_bus_(message_bus)
    , resource_manager_(resource_manager)
    , config_manager_(config_manager)
{
}

ArchitectureManager::~ArchitectureManager() {
    if (initialized_) {
        shutdown();
    }
}

bool ArchitectureManager::initialize() {
    if (initialized_) {
        return true;
    }

    try {
        // Initialize core services
        if (!message_bus_->initialize()) {
            return false;
        }

        if (!resource_manager_->initialize()) {
            return false;
        }

        if (!config_manager_->initialize()) {
            return false;
        }

        // Initialize default systems
        auto rendering_system = SystemFactory::create_rendering_system();
        if (!register_system(rendering_system)) {
            return false;
        }

        auto physics_system = SystemFactory::create_physics_system();
        if (!register_system(physics_system)) {
            return false;
        }

        auto interaction_system = SystemFactory::create_interaction_system();
        if (!register_system(interaction_system)) {
            return false;
        }

        auto animation_system = SystemFactory::create_animation_system();
        if (!register_system(animation_system)) {
            return false;
        }

        initialized_ = true;
        last_update_time_ = std::chrono::steady_clock::now();

        return true;

    } catch (const std::exception& e) {
        // Log error
        return false;
    }
}

void ArchitectureManager::shutdown() {
    if (!initialized_) {
        return;
    }

    // Shutdown all systems in reverse priority order
    auto systems = get_systems_by_priority();
    std::reverse(systems.begin(), systems.end());

    for (auto system : systems) {
        system->shutdown();
    }

    systems_.clear();

    // Shutdown core services
    if (message_bus_) {
        message_bus_->shutdown();
    }

    if (resource_manager_) {
        resource_manager_->shutdown();
    }

    if (config_manager_) {
        config_manager_->shutdown();
    }

    // Clear entities
    entities_.clear();

    initialized_ = false;
}

bool ArchitectureManager::is_initialized() const {
    return initialized_;
}

std::shared_ptr<IEntity> ArchitectureManager::create_entity(const std::string& name) {
    auto entity = std::make_shared<Entity>(name);
    uint64_t id = entity->get_id();

    {
        std::unique_lock<std::shared_mutex> lock(entities_mutex_);
        entities_[id] = entity;
    }

    // Notify systems of new entity
    for (const auto& [name, system] : systems_) {
        system->on_entity_created(entity);
    }

    return entity;
}

bool ArchitectureManager::destroy_entity(uint64_t entity_id) {
    std::shared_ptr<IEntity> entity;

    {
        std::shared_lock<std::shared_mutex> lock(entities_mutex_);
        auto it = entities_.find(entity_id);
        if (it == entities_.end()) {
            return false;
        }
        entity = it->second;
    }

    // Notify systems before destruction
    for (const auto& [name, system] : systems_) {
        system->on_entity_destroyed(entity_id);
    }

    {
        std::unique_lock<std::shared_mutex> lock(entities_mutex_);
        entities_.erase(entity_id);
    }

    return true;
}

std::shared_ptr<IEntity> ArchitectureManager::get_entity(uint64_t entity_id) const {
    std::shared_lock<std::shared_mutex> lock(entities_mutex_);
    auto it = entities_.find(entity_id);
    return (it != entities_.end()) ? it->second : nullptr;
}

std::vector<std::shared_ptr<IEntity>> ArchitectureManager::get_all_entities() const {
    std::shared_lock<std::shared_mutex> lock(entities_mutex_);
    std::vector<std::shared_ptr<IEntity>> result;
    result.reserve(entities_.size());

    for (const auto& [id, entity] : entities_) {
        result.push_back(entity);
    }

    return result;
}

std::vector<std::shared_ptr<IEntity>> ArchitectureManager::query_entities(
    const std::function<bool(const IEntity&)>& predicate) const {

    std::shared_lock<std::shared_mutex> lock(entities_mutex_);
    std::vector<std::shared_ptr<IEntity>> result;

    for (const auto& [id, entity] : entities_) {
        if (predicate(*entity)) {
            result.push_back(entity);
        }
    }

    return result;
}

bool ArchitectureManager::register_system(std::shared_ptr<ISystem> system) {
    if (!system) {
        return false;
    }

    const std::string& name = system->get_name();

    {
        std::unique_lock<std::shared_mutex> lock(systems_mutex_);
        if (systems_.find(name) != systems_.end()) {
            return false; // System already registered
        }
        systems_[name] = system;
    }

    // Initialize the system
    if (!system->initialize()) {
        systems_.erase(name);
        return false;
    }

    return true;
}

bool ArchitectureManager::unregister_system(const std::string& name) {
    std::shared_ptr<ISystem> system;

    {
        std::unique_lock<std::shared_mutex> lock(systems_mutex_);
        auto it = systems_.find(name);
        if (it == systems_.end()) {
            return false;
        }
        system = it->second;
        systems_.erase(it);
    }

    // Shutdown the system
    if (system) {
        system->shutdown();
    }

    return true;
}

std::shared_ptr<ISystem> ArchitectureManager::get_system(const std::string& name) const {
    std::shared_lock<std::shared_mutex> lock(systems_mutex_);
    auto it = systems_.find(name);
    return (it != systems_.end()) ? it->second : nullptr;
}

bool ArchitectureManager::load_plugin(const std::string& path) {
    // TODO: Implement plugin loading
    return false;
}

bool ArchitectureManager::unload_plugin(const std::string& name) {
    // TODO: Implement plugin unloading
    return false;
}

std::shared_ptr<IPlugin> ArchitectureManager::get_plugin(const std::string& name) const {
    std::shared_lock<std::shared_mutex> lock(plugins_mutex_);
    auto it = plugins_.find(name);
    return (it != plugins_.end()) ? it->second : nullptr;
}

std::vector<std::shared_ptr<IPlugin>> ArchitectureManager::get_loaded_plugins() const {
    std::shared_lock<std::shared_mutex> lock(plugins_mutex_);
    std::vector<std::shared_ptr<IPlugin>> result;
    result.reserve(plugins_.size());

    for (const auto& [name, plugin] : plugins_) {
        result.push_back(plugin);
    }

    return result;
}

void ArchitectureManager::update(double delta_time) {
    if (!initialized_) {
        return;
    }

    auto current_time = std::chrono::steady_clock::now();
    double actual_delta_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
        current_time - last_update_time_).count() / 1e9;
    last_update_time_ = current_time;

    // Update statistics
    update_statistics(actual_delta_time);

    // Update all entities
    {
        std::shared_lock<std::shared_mutex> lock(entities_mutex_);
        for (const auto& [id, entity] : entities_) {
            if (entity->is_active()) {
                entity->update(actual_delta_time);
            }
        }
    }

    // Update all systems in priority order
    auto systems = get_systems_by_priority();
    for (auto system : systems) {
        if (system->is_enabled()) {
            system->update(actual_delta_time);
        }
    }

    // Process pending events
    process_events();

    // Optimize memory usage periodically
    static double time_since_optimization = 0.0;
    time_since_optimization += actual_delta_time;
    if (time_since_optimization > 60.0) { // Every minute
        optimize_memory();
        time_since_optimization = 0.0;
    }
}

void ArchitectureManager::render() {
    if (!initialized_) {
        return;
    }

    // Render all entities
    {
        std::shared_lock<std::shared_mutex> lock(entities_mutex_);
        for (const auto& [id, entity] : entities_) {
            if (entity->is_active() && entity->is_visible()) {
                entity->render();
            }
        }
    }

    // Render all systems
    auto systems = get_systems_by_priority();
    for (auto system : systems) {
        if (system->is_enabled()) {
            system->render();
        }
    }
}

std::vector<std::shared_ptr<ISystem>> ArchitectureManager::get_systems_by_priority() const {
    std::shared_lock<std::shared_mutex> lock(systems_mutex_);
    std::vector<std::shared_ptr<ISystem>> result;
    result.reserve(systems_.size());

    for (const auto& [name, system] : systems_) {
        result.push_back(system);
    }

    // Sort by priority (lower numbers = higher priority)
    std::sort(result.begin(), result.end(),
        [](const std::shared_ptr<ISystem>& a, const std::shared_ptr<ISystem>& b) {
            return a->get_priority() < b->get_priority();
        });

    return result;
}

std::vector<std::shared_ptr<ISystem>> ArchitectureManager::get_enabled_systems() const {
    std::shared_lock<std::shared_mutex> lock(systems_mutex_);
    std::vector<std::shared_ptr<ISystem>> result;

    for (const auto& [name, system] : systems_) {
        if (system->is_enabled()) {
            result.push_back(system);
        }
    }

    return result;
}

void ArchitectureManager::set_target_frame_rate(double fps) {
    // Pass to rendering system
    auto rendering_system = get_system(RenderingSystem::SYSTEM_NAME);
    if (rendering_system) {
        // Cast and set frame rate
    }
}

void ArchitectureManager::set_memory_limit(size_t mb) {
    // Configure memory limits
}

void ArchitectureManager::set_thread_pool_size(size_t count) {
    // Configure threading
}

void ArchitectureManager::update_statistics(double delta_time) {
    // Update entity statistics
    {
        std::shared_lock<std::shared_mutex> lock(entities_mutex_);
        statistics_.total_entities = entities_.size();
        statistics_.active_entities = 0;
        statistics_.visible_entities = 0;

        for (const auto& [id, entity] : entities_) {
            if (entity->is_active()) {
                statistics_.active_entities++;
            }
            if (entity->is_visible()) {
                statistics_.visible_entities++;
            }
        }
    }

    // Update system statistics
    {
        std::shared_lock<std::shared_mutex> lock(systems_mutex_);
        statistics_.total_systems = systems_.size();
        statistics_.enabled_systems = 0;
        statistics_.average_system_update_time = 0.0;
        statistics_.average_system_render_time = 0.0;

        for (const auto& [name, system] : systems_) {
            if (system->is_enabled()) {
                statistics_.enabled_systems++;
                statistics_.average_system_update_time += system->get_average_update_time();
                statistics_.average_system_render_time += system->get_average_render_time();
            }
        }

        if (statistics_.enabled_systems > 0) {
            statistics_.average_system_update_time /= statistics_.enabled_systems;
            statistics_.average_system_render_time /= statistics_.enabled_systems;
        }
    }

    // Update performance statistics
    statistics_.frame_rate = 1.0 / delta_time;
    statistics_.average_frame_time = delta_time;
    statistics_.frame_time = delta_time;

    // Update memory statistics
    statistics_.total_memory_usage = 0; // TODO: Calculate actual memory usage
    statistics_.resource_memory_usage = 0;
    statistics_.entity_memory_usage = 0;
    statistics_.component_memory_usage = 0;
}

void ArchitectureManager::process_events() {
    // Process message bus events
    // Handle system events
    // Process entity events
}

void ArchitectureManager::optimize_memory() {
    // Optimize entity memory usage
    {
        std::shared_lock<std::shared_mutex> lock(entities_mutex_);
        for (const auto& [id, entity] : entities_) {
            // Call optimize_memory on entities that support it
        }
    }

    // Optimize system memory usage
    {
        std::shared_lock<std::shared_mutex> lock(systems_mutex_);
        for (const auto& [name, system] : systems_) {
            // Systems can optimize their own memory usage
        }
    }
}

} // namespace architecture
} // namespace gui3d
} // namespace hsml
