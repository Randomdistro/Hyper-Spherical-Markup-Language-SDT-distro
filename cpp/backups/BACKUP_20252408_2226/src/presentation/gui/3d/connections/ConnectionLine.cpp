/** @file ConnectionLine.h
 * @brief 3D GUI Connection lines between code entities
 *
 * Phase 5: Immersive 3D Development Environment
 * Visual representation of code relationships and dependencies.
 */

#pragma once

#include "hsml/domain/entities/SphericalEntity.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include <memory>
#include <string>
#include <chrono>
#include <vector>

namespace hsml {
namespace presentation {

// Forward declarations
class FloatingCodeEntity;

/**
 * @brief Types of connections between code entities
 */
enum class ConnectionType {
    FUNCTION_CALL,      // Function/method call
    INHERITANCE,        // Class inheritance
    COMPOSITION,        // Object composition
    DEPENDENCY,         // Import/include dependency
    DATA_FLOW,          // Data flow between functions
    EVENT_FLOW,         // Event propagation
    ASYNC_CALL,         // Asynchronous call
    TEMPLATE_USAGE,     // Template/generic usage
    INTERFACE_IMPLEMENTATION, // Interface implementation
    FRIEND_RELATIONSHIP       // Friend relationship
};

/**
 * @brief Represents a visual connection line between code entities
 *
 * Connection lines visualize relationships between different parts of the codebase,
 * showing function calls, dependencies, and data flow in real-time.
 */
class ConnectionLine : public domain::SphericalEntity {
public:
    /**
     * @brief Create a connection line
     *
     * @param from_entity Source entity
     * @param to_entity Target entity
     * @param connection_type Type of relationship
     * @param color Line color
     * @param thickness Line thickness
     */
    ConnectionLine(
        std::shared_ptr<FloatingCodeEntity> from_entity,
        std::shared_ptr<FloatingCodeEntity> to_entity,
        ConnectionType connection_type,
        const Vector3& color = Vector3(0.0f, 1.0f, 0.0f), // Default green
        float thickness = 2.0f
    );

    // Connection information
    ConnectionType get_connection_type() const { return connection_type_; }
    const std::string& get_description() const { return description_; }
    void set_description(const std::string& description) { description_ = description; }

    // Visual properties
    void set_color(const Vector3& color) { color_ = color; }
    const Vector3& get_color() const { return color_; }
    void set_thickness(float thickness) { thickness_ = thickness; }
    float get_thickness() const { return thickness_; }

    // Animation and effects
    void set_pulse_frequency(float hz) { pulse_frequency_ = hz; }
    float get_pulse_frequency() const { return pulse_frequency_; }
    void set_pulse_amplitude(float amplitude) { pulse_amplitude_ = amplitude; }
    float get_pulse_amplitude() const { return pulse_amplitude_; }

    // Data flow visualization
    void set_data_flow_direction(bool bidirectional) { bidirectional_ = bidirectional; }
    bool is_bidirectional() const { return bidirectional_; }
    void set_data_flow_rate(float rate) { data_flow_rate_ = rate; }
    float get_data_flow_rate() const { return data_flow_rate_; }

    // Connection strength and activity
    void set_activity_level(float level) { activity_level_ = level; }
    float get_activity_level() const { return activity_level_; }
    void set_max_activity_level(float level) { max_activity_level_ = level; }
    float get_max_activity_level() const { return max_activity_level_; }

    // Error and warning states
    void set_error_state(bool has_error) { has_error_ = has_error; }
    bool has_error() const { return has_error_; }
    void set_warning_state(bool has_warning) { has_warning_ = has_warning; }
    bool has_warning() const { return has_warning_; }

    // Connection lifecycle
    void activate();
    void deactivate();
    bool is_active() const { return is_active_; }

    // Entity references
    std::shared_ptr<FloatingCodeEntity> get_from_entity() const { return from_entity_.lock(); }
    std::shared_ptr<FloatingCodeEntity> get_to_entity() const { return to_entity_.lock(); }

    // Update and rendering
    void update(float delta_time) override;
    void render();

    // Event callbacks
    using ConnectionActivatedCallback = std::function<void(ConnectionLine*)>;
    using ConnectionDeactivatedCallback = std::function<void(ConnectionLine*)>;
    using DataFlowCallback = std::function<void(ConnectionLine*, float)>;

    void set_activated_callback(ConnectionActivatedCallback callback) {
        activated_callback_ = callback;
    }
    void set_deactivated_callback(ConnectionDeactivatedCallback callback) {
        deactivated_callback_ = callback;
    }
    void set_data_flow_callback(DataFlowCallback callback) {
        data_flow_callback_ = callback;
    }

private:
    // Connection endpoints
    std::weak_ptr<FloatingCodeEntity> from_entity_;
    std::weak_ptr<FloatingCodeEntity> to_entity_;

    // Connection properties
    ConnectionType connection_type_;
    std::string description_;
    Vector3 color_;
    float thickness_;
    bool is_active_{true};

    // Animation properties
    float pulse_frequency_{0.0f};
    float pulse_amplitude_{0.0f};
    float animation_time_{0.0f};
    float current_pulse_value_{0.0f};

    // Data flow properties
    bool bidirectional_{false};
    float data_flow_rate_{0.0f};
    float flow_animation_offset_{0.0f};

    // Activity monitoring
    float activity_level_{0.0f};
    float max_activity_level_{1.0f};
    std::chrono::steady_clock::time_point last_activity_time_;

    // Error states
    bool has_error_{false};
    bool has_warning_{false};

    // Visual effects
    struct Particle {
        Vector3 position;
        Vector3 velocity;
        float life;
        float max_life;
    };
    std::vector<Particle> particles_;
    size_t max_particles_{50};

    // Callbacks
    ConnectionActivatedCallback activated_callback_;
    ConnectionDeactivatedCallback deactivated_callback_;
    DataFlowCallback data_flow_callback_;

    // Helper methods
    void update_position();
    void update_animation(float delta_time);
    void update_particles(float delta_time);
    void generate_particles();
    Vector3 calculate_midpoint() const;
    Vector3 calculate_bezier_point(float t) const;
    void render_line();
    void render_particles();
    void render_data_flow();
    void render_error_effects();
    void render_activity_effects();

    // Connection type specific properties
    void apply_connection_type_styling();
    Vector3 get_connection_type_color() const;
    float get_connection_type_thickness() const;
};

/**
 * @brief Manager for connection lines
 */
class ConnectionLineManager {
public:
    static ConnectionLineManager& instance();

    // Connection creation
    std::shared_ptr<ConnectionLine> create_connection(
        std::shared_ptr<FloatingCodeEntity> from,
        std::shared_ptr<FloatingCodeEntity> to,
        ConnectionType type,
        const std::string& description = ""
    );

    // Connection management
    void remove_connection(std::shared_ptr<ConnectionLine> connection);
    void clear_all_connections();
    const std::vector<std::shared_ptr<ConnectionLine>>& get_all_connections() const;

    // Query connections
    std::vector<std::shared_ptr<ConnectionLine>> get_connections_from_entity(
        std::shared_ptr<FloatingCodeEntity> entity
    ) const;
    std::vector<std::shared_ptr<ConnectionLine>> get_connections_to_entity(
        std::shared_ptr<FloatingCodeEntity> entity
    ) const;
    std::vector<std::shared_ptr<ConnectionLine>> get_connections_of_type(
        ConnectionType type
    ) const;

    // Update all connections
    void update_all(float delta_time);
    void render_all();

    // Connection analysis
    struct ConnectionStatistics {
        size_t total_connections;
        size_t active_connections;
        size_t error_connections;
        size_t warning_connections;
        std::unordered_map<ConnectionType, size_t> connections_by_type;
        float average_activity_level;
        float max_activity_level;
    };

    ConnectionStatistics get_statistics() const;

private:
    ConnectionLineManager() = default;
    std::vector<std::shared_ptr<ConnectionLine>> connections_;

    mutable std::mutex mutex_;
};

/**
 * @brief Automatic connection analyzer
 */
class ConnectionAnalyzer {
public:
    static std::vector<std::shared_ptr<ConnectionLine>> analyze_code_dependencies(
        const std::vector<std::shared_ptr<FloatingCodeEntity>>& entities
    );

    static std::vector<std::shared_ptr<ConnectionLine>> analyze_function_calls(
        const std::vector<std::shared_ptr<FloatingCodeEntity>>& entities
    );

    static std::vector<std::shared_ptr<ConnectionLine>> analyze_inheritance_hierarchy(
        const std::vector<std::shared_ptr<FloatingCodeEntity>>& entities
    );

private:
    static ConnectionType determine_connection_type(
        const std::string& from_line,
        const std::string& to_entity_name
    );
    static bool is_function_call(const std::string& line, const std::string& function_name);
    static bool is_inheritance(const std::string& line, const std::string& base_class);
    static bool is_include_dependency(const std::string& line, const std::string& include_name);
};

} // namespace presentation
} // namespace hsml
