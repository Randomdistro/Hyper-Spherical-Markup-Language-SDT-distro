/**
 * Spatial Document Parser - C++20 Implementation
 * Revolutionary HSML document parser with multi-paradigm design
 * Integrates with existing HSML foundation for maximum performance
 */

#pragma once

#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <unordered_map>
#include <expected>
#include <concepts>
#include <functional>
#include <atomic>
#include <chrono>
#include <future>

// Existing HSML foundation
#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"
#include "hsml/parsers/hsml_parser.h"

namespace p0rt3r::parser {

// Forward declarations
class SpatialElement;
class SpatialStylesheet;
class SpatialScript;

/**
 * // [The OOP Architect]: Well-structured document representation
 * Revolutionary 3D document structure for spatial web pages
 */
class SpatialDocument {
public:
    // Document metadata
    struct DocumentMetadata {
        std::string document_uri;
        std::string title;
        std::string language = "en";
        std::string encoding = "utf-8";
        
        // Spatial properties
        hsml::core::SphericalCoords document_origin{650.0, 0.0, 0.0};
        double document_radius = 2000.0;
        
        // Performance hints
        uint32_t estimated_element_count = 0;
        uint32_t estimated_memory_usage_mb = 0;
        std::chrono::milliseconds estimated_load_time{0};
        
        // Dependencies
        std::vector<std::string> stylesheets;
        std::vector<std::string> scripts;
        std::vector<std::string> resources;
        
        DocumentMetadata(std::string_view uri) : document_uri(uri) {}
    };

private:
    DocumentMetadata metadata_;
    std::vector<std::unique_ptr<SpatialElement>> root_elements_;
    std::vector<std::unique_ptr<SpatialStylesheet>> stylesheets_;
    std::vector<std::unique_ptr<SpatialScript>> inline_scripts_;
    
    // Performance optimization
    mutable std::atomic<bool> is_parsed_{false};
    mutable std::atomic<bool> is_ready_{false};

public:
    explicit SpatialDocument(std::string_view uri) : metadata_(uri) {}
    
    // Document metadata access
    const DocumentMetadata& metadata() const noexcept { return metadata_; }
    const std::string& document_uri() const noexcept { return metadata_.document_uri; }
    const hsml::core::SphericalCoords& document_origin() const noexcept { return metadata_.document_origin; }
    
    // Element management
    void add_root_element(std::unique_ptr<SpatialElement> element);
    const std::vector<std::unique_ptr<SpatialElement>>& root_elements() const noexcept { return root_elements_; }
    
    // Dependencies
    const std::vector<std::string>& dependencies() const noexcept { return metadata_.resources; }
    const std::vector<std::unique_ptr<SpatialScript>>& inline_scripts() const noexcept { return inline_scripts_; }
    const std::vector<std::unique_ptr<SpatialStylesheet>>& stylesheets() const noexcept { return stylesheets_; }
    
    // Document state
    bool is_parsed() const noexcept { return is_parsed_.load(); }
    bool is_ready() const noexcept { return is_ready_.load(); }
    
    void mark_parsed() const noexcept { is_parsed_.store(true); }
    void mark_ready() const noexcept { is_ready_.store(true); }
};

/**
 * // [The Functional Purist]: Immutable spatial element representation
 * Represents a single spatial element in 3D space
 */
class SpatialElement {
public:
    enum class ElementType {
        BUBBLE,           // 3D bubble element
        PRESENCE,         // Physics-enabled presence
        PORTAL,           // Navigation portal
        FIELD,            // Force field
        MEDIA,            // 3D media element
        INTERFACE,        // 2D interface projected in 3D
        CUSTOM            // Custom element type
    };
    
    struct ElementProperties {
        std::string element_id;
        std::string element_name;
        ElementType type = ElementType::BUBBLE;
        
        // Spatial properties
        hsml::core::SphericalCoords position{0, 0, 0};
        double radius = 50.0;
        double mass = 1.0;
        
        // Visual properties
        struct Color { double r, g, b, a; } color{1.0, 1.0, 1.0, 1.0};
        double opacity = 1.0;
        bool visible = true;
        
        // Interaction properties
        bool interactive = true;
        double interaction_radius = 100.0;
        
        // Attributes and content
        std::unordered_map<std::string, std::string> attributes;
        std::string text_content;
        std::string html_content;
        
        ElementProperties(std::string_view id, ElementType elem_type) 
            : element_id(id), type(elem_type) {}
    };

private:
    ElementProperties properties_;
    std::vector<std::unique_ptr<SpatialElement>> children_;
    SpatialElement* parent_ = nullptr;

public:
    explicit SpatialElement(const ElementProperties& props) : properties_(props) {}
    
    // Element properties
    const ElementProperties& properties() const noexcept { return properties_; }
    const std::string& element_id() const noexcept { return properties_.element_id; }
    ElementType element_type() const noexcept { return properties_.type; }
    
    // Spatial properties
    const hsml::core::SphericalCoords& position() const noexcept { return properties_.position; }
    double radius() const noexcept { return properties_.radius; }
    
    // Hierarchy management
    void add_child(std::unique_ptr<SpatialElement> child);
    const std::vector<std::unique_ptr<SpatialElement>>& children() const noexcept { return children_; }
    SpatialElement* parent() const noexcept { return parent_; }
    
    // Content access
    const std::string& text_content() const noexcept { return properties_.text_content; }
    const std::unordered_map<std::string, std::string>& attributes() const noexcept { return properties_.attributes; }
};

/**
 * // [The Modern Hipster]: Cutting-edge spatial stylesheet support
 * CSSS (Cascading Spatial Style Sheets) parser and processor
 */
class SpatialStylesheet {
public:
    struct StyleRule {
        std::string selector;
        std::unordered_map<std::string, std::string> properties;
        
        // Spatial-specific properties
        std::optional<hsml::core::SphericalCoords> position;
        std::optional<double> radius;
        std::optional<double> mass;
        std::optional<double> orbital_velocity;
        
        StyleRule(std::string_view sel) : selector(sel) {}
    };

private:
    std::string stylesheet_source_;
    std::vector<StyleRule> rules_;
    bool is_parsed_ = false;

public:
    explicit SpatialStylesheet(std::string_view source) : stylesheet_source_(source) {}
    
    // Parse CSSS content
    [[nodiscard]] auto parse() -> std::expected<void, std::string>;
    
    // Rule access
    const std::vector<StyleRule>& rules() const noexcept { return rules_; }
    [[nodiscard]] std::vector<const StyleRule*> rules_for_element(const SpatialElement& element) const;
    
    // Stylesheet state
    bool is_parsed() const noexcept { return is_parsed_; }
    const std::string& source() const noexcept { return stylesheet_source_; }
};

/**
 * // [The Security Paranoid]: Safely contained spatial scripts
 * JavaScript integration with spatial APIs
 */
class SpatialScript {
public:
    enum class ScriptType {
        JAVASCRIPT,
        SHAPESCRIPT,    // HSML's domain-specific language
        WEBASSEMBLY,
        GLSL_COMPUTE    // GPU compute shaders
    };

private:
    std::string script_source_;
    ScriptType script_type_;
    std::string script_id_;
    bool is_compiled_ = false;

public:
    SpatialScript(std::string_view source, ScriptType type, std::string_view id = "")
        : script_source_(source), script_type_(type), script_id_(id) {}
    
    // Script properties
    const std::string& source() const noexcept { return script_source_; }
    ScriptType type() const noexcept { return script_type_; }
    const std::string& script_id() const noexcept { return script_id_; }
    
    // Compilation
    [[nodiscard]] auto compile() -> std::expected<void, std::string>;
    bool is_compiled() const noexcept { return is_compiled_; }
};

/**
 * // [The Performance Demon]: SIMD-optimized parsing engine
 * High-performance spatial document parser with multiple paradigms
 */
template<typename CharType = char>
requires std::same_as<CharType, char> || std::same_as<CharType, wchar_t>
class SpatialDocumentParser {
public:
    struct ParserConfig {
        bool enable_javascript = true;
        bool enable_css_parsing = true;
        bool enable_validation = true;
        bool enable_optimization = true;
        
        // Performance settings
        uint32_t max_parse_threads = std::thread::hardware_concurrency();
        size_t max_document_size_mb = 100;
        std::chrono::milliseconds parse_timeout{30000};
        
        // Spatial settings
        double default_viewer_distance = 650.0;
        double default_document_radius = 2000.0;
        hsml::core::SphericalCoords default_origin{650.0, 0.0, 0.0};
        
        ParserConfig() = default;
    };

private:    
    ParserConfig config_;
    std::unique_ptr<hsml::parsers::HSMLParser> base_parser_;
    
    // Performance metrics
    mutable std::atomic<uint64_t> documents_parsed_{0};
    mutable std::atomic<uint64_t> elements_parsed_{0};
    mutable std::atomic<double> average_parse_time_ms_{0.0};

    // // [The Functional Purist]: Pure parsing functions
    [[nodiscard]] auto parse_spatial_attributes(std::string_view attr_string) const 
        -> std::expected<std::unordered_map<std::string, std::string>, std::string>;
    
    [[nodiscard]] auto parse_spherical_coordinate(std::string_view coord_string) const
        -> std::expected<hsml::core::SphericalCoords, std::string>;
    
    [[nodiscard]] auto parse_solid_angle(std::string_view angle_string) const
        -> std::expected<hsml::core::SolidAngle, std::string>;

    // // [The Template Metaprogrammer]: Compile-time optimizations
    template<typename ElementCallback>
    requires std::invocable<ElementCallback, const SpatialElement&>
    void traverse_elements_compile_time(const SpatialDocument& doc, ElementCallback&& callback) const;

public:
    explicit SpatialDocumentParser(const ParserConfig& config = {}) : config_(config) {
        base_parser_ = std::make_unique<hsml::parsers::HSMLParser>();
    }
    
    // Main parsing interface
    [[nodiscard]] auto parse_document(std::string_view hsml_content, std::string_view document_uri)
        -> std::expected<std::unique_ptr<SpatialDocument>, std::string>;
    
    [[nodiscard]] auto parse_document_async(std::string_view hsml_content, std::string_view document_uri)
        -> std::future<std::expected<std::unique_ptr<SpatialDocument>, std::string>>;
    
    // Element parsing
    [[nodiscard]] auto parse_spatial_element(std::string_view element_markup)
        -> std::expected<std::unique_ptr<SpatialElement>, std::string>;
    
    // Stylesheet parsing
    [[nodiscard]] auto parse_spatial_stylesheet(std::string_view csss_content)
        -> std::expected<std::unique_ptr<SpatialStylesheet>, std::string>;
    
    // Script parsing
    [[nodiscard]] auto parse_spatial_script(std::string_view script_content, SpatialScript::ScriptType type)
        -> std::expected<std::unique_ptr<SpatialScript>, std::string>;
    
    // Configuration
    void update_config(const ParserConfig& new_config) { config_ = new_config; }
    const ParserConfig& config() const noexcept { return config_; }
    
    // Performance metrics
    struct ParserMetrics {
        uint64_t documents_parsed;
        uint64_t elements_parsed;
        double average_parse_time_ms;
        size_t memory_usage_mb;
        uint32_t active_parsing_threads;
    };
    
    [[nodiscard]] ParserMetrics get_metrics() const;
    
    // Validation
    [[nodiscard]] auto validate_document(const SpatialDocument& document) const
        -> std::vector<std::string>;
    
    [[nodiscard]] auto validate_hsml_syntax(std::string_view hsml_content) const
        -> std::expected<void, std::vector<std::string>>;
};

// // [The Hacktivist]: Quick utility functions that shouldn't work but do
namespace parser_utils {
    // Single-line parsers for common patterns
    [[nodiscard]] inline auto quick_parse_position(std::string_view pos_str) -> hsml::core::SphericalCoords {
        // Regex wizardry: "spherical(r, theta, phi)" -> coordinates
        static const std::regex coord_pattern(R"(spherical\s*\(\s*([0-9.]+)\s*,\s*([0-9.]+)\s*,\s*([0-9.]+)\s*\))");
        std::smatch matches;
        std::string str(pos_str);
        
        return std::regex_match(str, matches, coord_pattern) 
            ? hsml::core::SphericalCoords{std::stod(matches[1]), std::stod(matches[2]), std::stod(matches[3])}
            : hsml::core::SphericalCoords{650.0, 0.0, 0.0}; // Fallback
    }
    
    [[nodiscard]] inline auto extract_element_id(std::string_view markup) -> std::string {
        // Quick and dirty ID extraction
        auto id_pos = markup.find("id=\"");
        if (id_pos == std::string_view::npos) return "";
        
        id_pos += 4; // Skip 'id="'
        auto end_pos = markup.find("\"", id_pos);
        return end_pos != std::string_view::npos 
            ? std::string(markup.substr(id_pos, end_pos - id_pos))
            : "";
    }
    
    [[nodiscard]] inline auto count_spatial_elements(std::string_view hsml) -> size_t {
        // Count occurrences of spatial elements
        size_t count = 0;
        const std::vector<std::string_view> spatial_tags = {"<bubble", "<presence", "<portal", "<field"};
        
        for (auto tag : spatial_tags) {
            size_t pos = 0;
            while ((pos = hsml.find(tag, pos)) != std::string_view::npos) {
                ++count;
                pos += tag.length();
            }
        }
        
        return count;
    }
}

// Factory functions
[[nodiscard]] auto create_spatial_document_parser(const SpatialDocumentParser<>::ParserConfig& config = {})
    -> std::unique_ptr<SpatialDocumentParser<>>;

} // namespace p0rt3r::parser

// Include implementation details
#include "spatial_document_parser_impl.h"