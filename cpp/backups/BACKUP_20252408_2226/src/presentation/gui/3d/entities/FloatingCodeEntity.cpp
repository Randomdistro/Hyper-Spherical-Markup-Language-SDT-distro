/** @file FloatingCodeEntity.h
 * @brief 3D GUI Entity representing floating code files
 *
 * Phase 5: Immersive 3D Development Environment
 * Core component for the floating code visualization system.
 */

#pragma once

#include "hsml/domain/entities/SphericalEntity.h"
#include "hsml/domain/entities/SphereShape.h"
#include "hsml/domain/entities/LinearMovementStrategy.h"
#include "hsml/core/spherical_coords.h"
#include "hsml/core/vector3.h"
#include <string>
#include <vector>
#include <memory>
#include <chrono>
#include <functional>

namespace hsml {
namespace presentation {

// Forward declarations
class ConnectionLine;
class CodeExecutionVisualizer;

/**
 * @brief Represents a code file floating in 3D space
 *
 * This entity visualizes source code files as floating 3D objects
 * with interactive capabilities for code navigation and editing.
 */
class FloatingCodeEntity : public domain::SphericalEntity {
public:
    /**
     * @brief Create a floating code entity
     *
     * @param filename The source file name
     * @param content The file content
     * @param position Spherical position in 3D space
     * @param radius Visual size of the entity
     */
    FloatingCodeEntity(
        const std::string& filename,
        const std::string& content,
        const SphericalCoords& position,
        double radius = 2.0
    );

    // Entity identification
    const std::string& get_filename() const { return filename_; }
    const std::string& get_file_extension() const { return extension_; }
    const std::string& get_language() const { return language_; }
    size_t get_line_count() const { return lines_.size(); }

    // Content management
    void set_content(const std::string& content);
    const std::string& get_content() const { return content_; }
    const std::vector<std::string>& get_lines() const { return lines_; }
    std::string get_line(size_t line_number) const;

    // Syntax highlighting
    void set_syntax_highlighting(bool enabled) { syntax_highlighting_ = enabled; }
    bool is_syntax_highlighting_enabled() const { return syntax_highlighting_; }

    // Visual properties
    void set_base_color(const Vector3& color) { base_color_ = color; }
    const Vector3& get_base_color() const { return base_color_; }
    void set_glow_intensity(float intensity) { glow_intensity_ = intensity; }
    float get_glow_intensity() const { return glow_intensity_; }

    // Interactive features
    void set_selected(bool selected);
    bool is_selected() const { return is_selected_; }
    void highlight_line(size_t line_number, const Vector3& highlight_color = Vector3(1.0f, 1.0f, 0.0f));
    void clear_highlights();

    // Connection management
    void add_connection(std::shared_ptr<ConnectionLine> connection);
    void remove_connection(std::shared_ptr<ConnectionLine> connection);
    const std::vector<std::shared_ptr<ConnectionLine>>& get_connections() const { return connections_; }

    // Function and class detection
    std::vector<std::pair<std::string, size_t>> find_functions() const;
    std::vector<std::pair<std::string, size_t>> find_classes() const;
    std::vector<std::pair<std::string, size_t>> find_includes() const;

    // Navigation
    void scroll_to_line(size_t line_number);
    void focus_on_entity();
    void zoom_to_entity(float duration = 1.0f);

    // Animation and interaction
    void update(float delta_time) override;
    void on_mouse_enter();
    void on_mouse_exit();
    void on_click();
    void on_double_click();

    // Event callbacks
    using ContentChangedCallback = std::function<void(const std::string&)>;
    using LineHighlightedCallback = std::function<void(size_t)>;
    using EntitySelectedCallback = std::function<void(FloatingCodeEntity*)>;

    void set_content_changed_callback(ContentChangedCallback callback) {
        content_changed_callback_ = callback;
    }
    void set_line_highlighted_callback(LineHighlightedCallback callback) {
        line_highlighted_callback_ = callback;
    }
    void set_selected_callback(EntitySelectedCallback callback) {
        selected_callback_ = callback;
    }

private:
    // File information
    std::string filename_;
    std::string extension_;
    std::string language_;
    std::string content_;
    std::vector<std::string> lines_;

    // Visual state
    Vector3 base_color_{0.2f, 0.6f, 1.0f}; // Default blue
    float glow_intensity_{0.0f};
    bool is_selected_{false};
    bool syntax_highlighting_{true};

    // Interactive state
    bool is_hovered_{false};
    float hover_time_{0.0f};
    float selection_time_{0.0f};

    // Line highlighting
    struct LineHighlight {
        size_t line_number;
        Vector3 color;
        float intensity;
        std::chrono::steady_clock::time_point start_time;
    };
    std::vector<LineHighlight> highlights_;

    // Connections to other entities
    std::vector<std::shared_ptr<ConnectionLine>> connections_;

    // Animation state
    Vector3 target_scale_{1.0f, 1.0f, 1.0f};
    Vector3 current_scale_{1.0f, 1.0f, 1.0f};
    float scale_animation_time_{0.0f};

    // Callbacks
    ContentChangedCallback content_changed_callback_;
    LineHighlightedCallback line_highlighted_callback_;
    EntitySelectedCallback selected_callback_;

    // Helper methods
    void parse_content();
    void detect_language();
    void update_visual_state(float delta_time);
    void animate_scale(const Vector3& target_scale, float duration);
    Vector3 get_language_color() const;
    std::string extract_function_name(const std::string& line) const;
    std::string extract_class_name(const std::string& line) const;
    std::string extract_include_name(const std::string& line) const;
};

/**
 * @brief Factory for creating floating code entities
 */
class FloatingCodeEntityFactory {
public:
    static std::shared_ptr<FloatingCodeEntity> create_from_file(
        const std::string& filepath,
        const SphericalCoords& position
    );

    static std::shared_ptr<FloatingCodeEntity> create_from_content(
        const std::string& filename,
        const std::string& content,
        const SphericalCoords& position
    );

    static std::vector<std::shared_ptr<FloatingCodeEntity>> create_from_project(
        const std::string& project_root,
        const SphericalCoords& center_position,
        double radius = 50.0
    );

private:
    static std::string read_file_content(const std::string& filepath);
    static SphericalCoords calculate_file_position(
        size_t file_index,
        size_t total_files,
        const SphericalCoords& center,
        double radius
    );
};

} // namespace presentation
} // namespace hsml
