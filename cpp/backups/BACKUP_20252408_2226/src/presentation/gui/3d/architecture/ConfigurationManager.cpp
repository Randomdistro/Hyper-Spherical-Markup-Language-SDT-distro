/** @file ConfigurationManager.h
 * @brief Configuration Manager Implementation - Runtime Configuration System
 */

#pragma once

#include "src/presentation/gui/3d/architecture/ArchitectureCore.h"
#include <unordered_map>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <string>
#include <variant>
#include <vector>
#include <chrono>

namespace hsml {
namespace gui3d {
namespace architecture {

/**
 * @brief Configuration Value Type
 */
using ConfigValue = std::variant<bool, int, double, std::string>;

/**
 * @brief Configuration Section
 */
class ConfigurationSection {
public:
    ConfigurationSection(const std::string& name);
    ~ConfigurationSection() = default;

    const std::string& get_name() const { return name_; }

    bool has_value(const std::string& key) const;
    ConfigValue get_value(const std::string& key) const;
    void set_value(const std::string& key, const ConfigValue& value);

    std::vector<std::string> get_keys() const;
    void clear();

    std::string serialize() const;
    bool deserialize(const std::string& data);

private:
    std::string name_;
    std::unordered_map<std::string, ConfigValue> values_;
    mutable std::shared_mutex mutex_;
};

/**
 * @brief Configuration Manager Implementation
 */
class ConfigurationManager : public IConfigurationManager {
public:
    ConfigurationManager();
    ~ConfigurationManager() override;

    bool initialize() override;
    void shutdown() override;

    // Configuration loading/saving
    bool load_config(const std::string& path) override;
    bool save_config(const std::string& path) const override;

    // Value access with type safety
    bool get_bool(const std::string& key, bool default_value = false) const override;
    int get_int(const std::string& key, int default_value = 0) const override;
    double get_double(const std::string& key, double default_value = 0.0) const override;
    std::string get_string(const std::string& key, const std::string& default_value = "") const override;

    bool set_bool(const std::string& key, bool value) override;
    bool set_int(const std::string& key, int value) override;
    bool set_double(const std::string& key, double value) override;
    bool set_string(const std::string& key, const std::string& value) override;

    // Section management
    std::vector<std::string> get_sections() const override;
    std::vector<std::string> get_keys(const std::string& section) const override;

    // Runtime configuration updates
    void reload_config() override;
    bool is_config_modified() const override;

    // Advanced features
    bool create_section(const std::string& section_name);
    bool delete_section(const std::string& section_name);
    bool has_section(const std::string& section_name) const;

    // Configuration validation
    bool validate_config() const;
    std::vector<std::string> get_validation_errors() const;

    // Configuration templates
    bool load_template(const std::string& template_name);
    bool save_template(const std::string& template_name) const;
    std::vector<std::string> get_available_templates() const;

    // Configuration events
    using ConfigChangeCallback = std::function<void(const std::string&, const ConfigValue&)>;
    uint64_t add_change_listener(const std::string& key, ConfigChangeCallback callback);
    bool remove_change_listener(uint64_t listener_id);

private:
    struct ConfigListener {
        uint64_t id;
        std::string key;
        ConfigChangeCallback callback;
        bool active;
    };

    // Configuration data
    std::unordered_map<std::string, std::shared_ptr<ConfigurationSection>> sections_;
    std::string current_config_path_;
    std::chrono::steady_clock::time_point last_load_time_;
    std::chrono::steady_clock::time_point last_save_time_;

    // Change listeners
    std::unordered_map<uint64_t, ConfigListener> listeners_;
    std::atomic<uint64_t> next_listener_id_;

    // Thread safety
    mutable std::shared_mutex config_mutex_;
    mutable std::shared_mutex listeners_mutex_;

    // State management
    std::atomic<bool> initialized_;
    std::atomic<bool> modified_;

    // Default configuration
    void load_defaults();
    void create_default_sections();

    // Helper methods
    std::pair<std::string, std::string> parse_key(const std::string& key) const;
    std::string make_key(const std::string& section, const std::string& name) const;

    template<typename T>
    T get_value(const std::string& key, T default_value) const;

    template<typename T>
    bool set_value(const std::string& key, T value);

    void notify_listeners(const std::string& key, const ConfigValue& value);

    // File I/O helpers
    bool load_from_file(const std::string& path);
    bool save_to_file(const std::string& path) const;
    bool parse_ini_file(const std::string& content);
    std::string generate_ini_content() const;

    // Validation helpers
    bool validate_section(const std::string& section_name) const;
    bool validate_key(const std::string& section, const std::string& key) const;
};

/**
 * @brief Configuration Templates
 */
class ConfigurationTemplates {
public:
    static const std::string DEFAULT_TEMPLATE;
    static const std::string DEBUG_TEMPLATE;
    static const std::string PERFORMANCE_TEMPLATE;
    static const std::string VR_TEMPLATE;

    static std::shared_ptr<ConfigurationManager> create_from_template(
        const std::string& template_name
    );

private:
    static void apply_default_template(ConfigurationManager& config);
    static void apply_debug_template(ConfigurationManager& config);
    static void apply_performance_template(ConfigurationManager& config);
    static void apply_vr_template(ConfigurationManager& config);
};

/**
 * @brief Configuration Validator
 */
class ConfigurationValidator {
public:
    static bool validate(const ConfigurationManager& config);
    static std::vector<std::string> get_validation_errors(const ConfigurationManager& config);

    static bool is_valid_key(const std::string& key);
    static bool is_valid_section(const std::string& section);
    static bool is_valid_value_type(const ConfigValue& value);

private:
    static bool validate_display_settings(const ConfigurationManager& config);
    static bool validate_performance_settings(const ConfigurationManager& config);
    static bool validate_system_settings(const ConfigurationManager& config);
};

} // namespace architecture
} // namespace gui3d
} // namespace hsml
