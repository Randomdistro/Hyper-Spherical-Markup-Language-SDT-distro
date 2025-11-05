#pragma once

#include "hsml/core/bubble.h"
#include "hsml/core/presence.h"
#include "hsml_semantic_analyzer.h"
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <functional>

namespace hsml {
namespace parsers {

using core::Bubble;
using core::BubblePtr;
using core::Presence;
using core::PresencePtr;
using core::Material;
using core::MaterialState;
using core::StateTensor;
using core::SphericalCoords;
using core::Vector3;

enum class TokenType {
    OPEN_TAG,           // <
    CLOSE_TAG,          // >
    OPEN_CLOSE_TAG,     // </
    SELF_CLOSING,       // />
    IDENTIFIER,         // bubble, presence, etc.
    ATTRIBUTE_NAME,     // r, theta, phi, etc.
    ATTRIBUTE_VALUE,    // "1000", "0.5", etc.
    TEXT,               // Content text
    COMMENT,            // <!-- -->
    WHITESPACE,
    END_OF_FILE,
    INVALID
};

struct Token {
    TokenType type;
    std::string value;
    size_t line;
    size_t column;
    
    Token(TokenType t = TokenType::INVALID, const std::string& v = "", size_t l = 0, size_t c = 0)
        : type(t), value(v), line(l), column(c) {}
};

struct ParseError {
    std::string message;
    size_t line;
    size_t column;
    
    ParseError(const std::string& msg, size_t l, size_t c)
        : message(msg), line(l), column(c) {}
    
    std::string to_string() const {
        return "Parse Error at line " + std::to_string(line) + 
               ", column " + std::to_string(column) + ": " + message;
    }
};

class HSMLLexer {
public:
    explicit HSMLLexer(const std::string& input);
    
    Token next_token();
    Token peek_token();
    bool has_more_tokens() const;
    
    size_t current_line() const { return current_line_; }
    size_t current_column() const { return current_column_; }

private:
    std::string input_;
    size_t position_;
    size_t current_line_;
    size_t current_column_;
    Token peeked_token_;
    bool has_peeked_;
    
    char current_char() const;
    char peek_char(size_t offset = 1) const;
    void advance();
    void skip_whitespace();
    Token read_string_literal();
    Token read_identifier();
    Token read_number();
    Token read_comment();
    bool is_alpha(char c) const;
    bool is_digit(char c) const;
    bool is_alnum(char c) const;
};

struct AttributeMap {
    std::unordered_map<std::string, std::string> values;
    
    bool has(const std::string& name) const {
        return values.find(name) != values.end();
    }
    
    std::string get(const std::string& name, const std::string& default_value = "") const {
        auto it = values.find(name);
        return it != values.end() ? it->second : default_value;
    }
    
    double get_double(const std::string& name, double default_value = 0.0) const {
        auto it = values.find(name);
        if (it != values.end()) {
            try {
                return std::stod(it->second);
            } catch (...) {
                return default_value;
            }
        }
        return default_value;
    }
    
    bool get_bool(const std::string& name, bool default_value = false) const {
        auto it = values.find(name);
        if (it != values.end()) {
            const std::string& value = it->second;
            return value == "true" || value == "1" || value == "yes";
        }
        return default_value;
    }
};

struct HSMLElement : public std::enable_shared_from_this<HSMLElement> {
    std::string tag_name;
    AttributeMap attributes;
    std::string text_content;
    std::vector<std::shared_ptr<HSMLElement>> children;
    std::weak_ptr<HSMLElement> parent;
    
    HSMLElement(const std::string& name) : tag_name(name) {}
    
    void add_child(std::shared_ptr<HSMLElement> child) {
        children.push_back(child);
        child->parent = shared_from_this();
    }
};

class HSMLParser {
public:
    HSMLParser();
    ~HSMLParser() = default;
    
    // Main parsing methods with semantic analysis
    BubblePtr parse(const std::string& hsml_content);
    BubblePtr parse_file(const std::string& filename);
    
    // Parse to element tree first (for validation/debugging)
    std::shared_ptr<HSMLElement> parse_to_elements(const std::string& hsml_content);
    
    // Convert element tree to bubble hierarchy with semantic validation
    BubblePtr elements_to_bubbles(std::shared_ptr<HSMLElement> root_element);
    
    // Semantic analysis integration
    HSMLSemanticAnalyzer& semantic_analyzer() { return semantic_analyzer_; }
    const HSMLSemanticAnalyzer& semantic_analyzer() const { return semantic_analyzer_; }
    
    // Validate parsed content semantically
    bool validate_semantics(std::shared_ptr<HSMLElement> root_element);
    // Semantic results accessors (no ValidationResult type)
    const std::vector<SemanticError>& semantic_errors() const { return semantic_analyzer_.get_errors(); }
    bool has_semantic_errors() const { return semantic_analyzer_.has_errors(); }
    
    // Error handling (includes both parse and semantic errors)
    const std::vector<ParseError>& errors() const { return errors_; }
    bool has_errors() const { return !errors_.empty() || semantic_analyzer_.has_errors(); }
    void clear_errors() { 
        errors_.clear(); 
        semantic_analyzer_.clear_errors();
    }
    
    // Get all errors (parse + semantic)
    std::vector<ParseError> get_all_errors() const;
    
    // Configuration
    void set_strict_mode(bool strict) { strict_mode_ = strict; }
    bool is_strict_mode() const { return strict_mode_; }
    
    void set_default_viewer_distance(double distance) { default_viewer_distance_ = distance; }
    double default_viewer_distance() const { return default_viewer_distance_; }
    
    // Semantic analysis configuration
    void enable_semantic_validation(bool enable) { semantic_validation_enabled_ = enable; }
    bool is_semantic_validation_enabled() const { return semantic_validation_enabled_; }

private:
    HSMLLexer lexer_;
    std::vector<ParseError> errors_;
    bool strict_mode_;
    double default_viewer_distance_;
    HSMLSemanticAnalyzer semantic_analyzer_;
    bool semantic_validation_enabled_;
    
    // Parser state
    Token current_token_;
    bool initialized_;
    
    // Parsing methods
    void init_parsing(const std::string& content);
    Token consume_token(TokenType expected_type);
    bool match_token(TokenType type);
    void advance_token();
    
    std::shared_ptr<HSMLElement> parse_element();
    AttributeMap parse_attributes();
    std::string parse_text_content();
    
    // Conversion methods
    BubblePtr convert_bubble_element(std::shared_ptr<HSMLElement> element, BubblePtr parent = nullptr);
    PresencePtr convert_presence_element(std::shared_ptr<HSMLElement> element);
    Material convert_material_element(std::shared_ptr<HSMLElement> element);
    StateTensor convert_state_element(std::shared_ptr<HSMLElement> element);
    
    SphericalCoords parse_coordinates(const AttributeMap& attributes);
    Vector3 parse_vector3(const AttributeMap& attributes, const std::string& prefix = "");
    
    // Validation methods
    void validate_element(std::shared_ptr<HSMLElement> element);
    void validate_bubble_attributes(const AttributeMap& attributes);
    void validate_presence_attributes(const AttributeMap& attributes);
    
    // Error reporting
    void report_error(const std::string& message);
    void report_error(const std::string& message, size_t line, size_t column);
    
    // Helper methods
    MaterialState parse_material_state(const std::string& state_string);
    std::string material_state_to_string(MaterialState state);
    
    bool is_valid_tag_name(const std::string& name) const;
    bool is_valid_attribute_name(const std::string& name) const;
};

// Factory methods for creating parsers with different configurations
namespace factory {

std::unique_ptr<HSMLParser> create_strict_parser(double viewer_distance = 650.0);
std::unique_ptr<HSMLParser> create_lenient_parser(double viewer_distance = 650.0);
std::unique_ptr<HSMLParser> create_debug_parser(double viewer_distance = 650.0);

} // namespace factory

// Utility functions
namespace utils {

std::string bubble_to_hsml(BubblePtr bubble, int indent_level = 0);
std::string presence_to_hsml(PresencePtr presence, int indent_level = 0);
bool validate_hsml_syntax(const std::string& hsml_content);
std::vector<std::string> extract_hsml_errors(const std::string& hsml_content);

} // namespace utils

} // namespace parsers
} // namespace hsml