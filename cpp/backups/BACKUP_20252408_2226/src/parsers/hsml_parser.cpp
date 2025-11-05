#include "hsml/parsers/hsml_parser.h"
#include <fstream>
#include <sstream>
#include <cctype>
#include <algorithm>

namespace hsml {
namespace parsers {

// HSMLLexer Implementation
HSMLLexer::HSMLLexer(const std::string& input) 
    : input_(input), position_(0), current_line_(1), current_column_(1), has_peeked_(false) {}

Token HSMLLexer::next_token() {
    if (has_peeked_) {
        has_peeked_ = false;
        return peeked_token_;
    }
    
    skip_whitespace();
    
    if (position_ >= input_.length()) {
        return Token(TokenType::END_OF_FILE, "", current_line_, current_column_);
    }
    
    char c = current_char();
    size_t start_line = current_line_;
    size_t start_column = current_column_;
    
    // Handle different token types
    if (c == '<') {
        advance();
        if (current_char() == '/') {
            advance();
            return Token(TokenType::OPEN_CLOSE_TAG, "</", start_line, start_column);
        } else if (current_char() == '!') {
            return read_comment();
        } else {
            return Token(TokenType::OPEN_TAG, "<", start_line, start_column);
        }
    }
    
    if (c == '>') {
        advance();
        return Token(TokenType::CLOSE_TAG, ">", start_line, start_column);
    }
    
    if (c == '/' && peek_char() == '>') {
        advance();
        advance();
        return Token(TokenType::SELF_CLOSING, "/>", start_line, start_column);
    }
    
    if (c == '"' || c == '\'') {
        return read_string_literal();
    }
    
    if (is_alpha(c) || c == '_') {
        return read_identifier();
    }
    
    if (is_digit(c) || c == '.' || c == '-' || c == '+') {
        return read_number();
    }
    
    // Single character token
    advance();
    return Token(TokenType::INVALID, std::string(1, c), start_line, start_column);
}

Token HSMLLexer::peek_token() {
    if (!has_peeked_) {
        peeked_token_ = next_token();
        has_peeked_ = true;
    }
    return peeked_token_;
}

bool HSMLLexer::has_more_tokens() const {
    return position_ < input_.length();
}

char HSMLLexer::current_char() const {
    return position_ < input_.length() ? input_[position_] : '\0';
}

char HSMLLexer::peek_char(size_t offset) const {
    size_t peek_pos = position_ + offset;
    return peek_pos < input_.length() ? input_[peek_pos] : '\0';
}

void HSMLLexer::advance() {
    if (position_ < input_.length()) {
        if (input_[position_] == '\n') {
            current_line_++;
            current_column_ = 1;
        } else {
            current_column_++;
        }
        position_++;
    }
}

void HSMLLexer::skip_whitespace() {
    while (position_ < input_.length() && std::isspace(current_char())) {
        advance();
    }
}

Token HSMLLexer::read_string_literal() {
    char quote_char = current_char();
    size_t start_line = current_line_;
    size_t start_column = current_column_;
    
    advance(); // Skip opening quote
    
    std::string value;
    while (position_ < input_.length() && current_char() != quote_char) {
        if (current_char() == '\\') {
            advance();
            if (position_ < input_.length()) {
                char escaped = current_char();
                switch (escaped) {
                    case 'n': value += '\n'; break;
                    case 't': value += '\t'; break;
                    case 'r': value += '\r'; break;
                    case '\\': value += '\\'; break;
                    case '"': value += '"'; break;
                    case '\'': value += '\''; break;
                    default: value += escaped; break;
                }
                advance();
            }
        } else {
            value += current_char();
            advance();
        }
    }
    
    if (current_char() == quote_char) {
        advance(); // Skip closing quote
    }
    
    return Token(TokenType::ATTRIBUTE_VALUE, value, start_line, start_column);
}

Token HSMLLexer::read_identifier() {
    size_t start_line = current_line_;
    size_t start_column = current_column_;
    
    std::string value;
    while (position_ < input_.length() && (is_alnum(current_char()) || current_char() == '_' || current_char() == '-')) {
        value += current_char();
        advance();
    }
    
    return Token(TokenType::IDENTIFIER, value, start_line, start_column);
}

Token HSMLLexer::read_number() {
    size_t start_line = current_line_;
    size_t start_column = current_column_;
    
    std::string value;
    bool has_dot = false;
    
    if (current_char() == '-' || current_char() == '+') {
        value += current_char();
        advance();
    }
    
    while (position_ < input_.length() && 
           (is_digit(current_char()) || (current_char() == '.' && !has_dot))) {
        if (current_char() == '.') {
            has_dot = true;
        }
        value += current_char();
        advance();
    }
    
    return Token(TokenType::ATTRIBUTE_VALUE, value, start_line, start_column);
}

Token HSMLLexer::read_comment() {
    size_t start_line = current_line_;
    size_t start_column = current_column_;
    
    std::string value;
    advance(); // Skip '!'
    
    if (current_char() == '-' && peek_char() == '-') {
        advance();
        advance();
        
        // Read until -->
        while (position_ + 2 < input_.length()) {
            if (current_char() == '-' && peek_char() == '-' && peek_char(2) == '>') {
                advance();
                advance();
                advance();
                break;
            }
            value += current_char();
            advance();
        }
    }
    
    return Token(TokenType::COMMENT, value, start_line, start_column);
}

bool HSMLLexer::is_alpha(char c) const {
    return std::isalpha(c);
}

bool HSMLLexer::is_digit(char c) const {
    return std::isdigit(c);
}

bool HSMLLexer::is_alnum(char c) const {
    return std::isalnum(c);
}

// HSMLParser Implementation
HSMLParser::HSMLParser() 
    : lexer_(""), strict_mode_(false), default_viewer_distance_(650.0), initialized_(false) {}

BubblePtr HSMLParser::parse(const std::string& hsml_content) {
    clear_errors();
    
    try {
        auto root_element = parse_to_elements(hsml_content);
        if (!root_element) {
            report_error("Failed to parse HSML content");
            return nullptr;
        }
        
        return elements_to_bubbles(root_element);
    } catch (const std::exception& e) {
        report_error(std::string("Parse exception: ") + e.what());
        return nullptr;
    }
}

BubblePtr HSMLParser::parse_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        report_error("Cannot open file: " + filename);
        return nullptr;
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    return parse(buffer.str());
}

std::shared_ptr<HSMLElement> HSMLParser::parse_to_elements(const std::string& hsml_content) {
    init_parsing(hsml_content);
    
    if (!initialized_) {
        return nullptr;
    }
    
    // Parse root element
    auto root = parse_element();
    
    if (current_token_.type != TokenType::END_OF_FILE) {
        report_error("Unexpected content after root element");
    }
    
    return root;
}

BubblePtr HSMLParser::elements_to_bubbles(std::shared_ptr<HSMLElement> root_element) {
    if (!root_element) {
        return nullptr;
    }
    
    return convert_bubble_element(root_element);
}

void HSMLParser::init_parsing(const std::string& content) {
    lexer_ = HSMLLexer(content);
    advance_token();
    initialized_ = true;
}

Token HSMLParser::consume_token(TokenType expected_type) {
    Token token = current_token_;
    if (token.type != expected_type) {
        report_error("Expected " + std::to_string(static_cast<int>(expected_type)) + 
                    " but got " + std::to_string(static_cast<int>(token.type)));
    }
    advance_token();
    return token;
}

bool HSMLParser::match_token(TokenType type) {
    return current_token_.type == type;
}

void HSMLParser::advance_token() {
    current_token_ = lexer_.next_token();
    
    // Skip whitespace and comments
    while (current_token_.type == TokenType::WHITESPACE || current_token_.type == TokenType::COMMENT) {
        current_token_ = lexer_.next_token();
    }
}

std::shared_ptr<HSMLElement> HSMLParser::parse_element() {
    if (!match_token(TokenType::OPEN_TAG)) {
        report_error("Expected opening tag");
        return nullptr;
    }
    
    consume_token(TokenType::OPEN_TAG);
    
    if (!match_token(TokenType::IDENTIFIER)) {
        report_error("Expected element name");
        return nullptr;
    }
    
    Token name_token = consume_token(TokenType::IDENTIFIER);
    auto element = std::make_shared<HSMLElement>(name_token.value);
    
    // Parse attributes
    element->attributes = parse_attributes();
    
    // Check for self-closing tag
    if (match_token(TokenType::SELF_CLOSING)) {
        consume_token(TokenType::SELF_CLOSING);
        return element;
    }
    
    consume_token(TokenType::CLOSE_TAG);
    
    // Parse content (children or text)
    while (!match_token(TokenType::OPEN_CLOSE_TAG) && !match_token(TokenType::END_OF_FILE)) {
        if (match_token(TokenType::OPEN_TAG)) {
            auto child = parse_element();
            if (child) {
                element->add_child(child);
            }
        } else if (match_token(TokenType::IDENTIFIER) || match_token(TokenType::ATTRIBUTE_VALUE)) {
            // Text content
            element->text_content += current_token_.value;
            advance_token();
        } else {
            advance_token(); // Skip unknown tokens
        }
    }
    
    // Parse closing tag
    if (match_token(TokenType::OPEN_CLOSE_TAG)) {
        consume_token(TokenType::OPEN_CLOSE_TAG);
        Token closing_name = consume_token(TokenType::IDENTIFIER);
        
        if (closing_name.value != name_token.value) {
            report_error("Mismatched closing tag: expected " + name_token.value + 
                        " but got " + closing_name.value);
        }
        
        consume_token(TokenType::CLOSE_TAG);
    }
    
    if (strict_mode_) {
        validate_element(element);
    }
    
    return element;
}

AttributeMap HSMLParser::parse_attributes() {
    AttributeMap attributes;
    
    while (match_token(TokenType::IDENTIFIER)) {
        Token name_token = consume_token(TokenType::IDENTIFIER);
        
        if (current_token_.value == "=") {
            advance_token(); // Skip '='
            
            if (match_token(TokenType::ATTRIBUTE_VALUE)) {
                Token value_token = consume_token(TokenType::ATTRIBUTE_VALUE);
                attributes.values[name_token.value] = value_token.value;
            } else {
                report_error("Expected attribute value");
            }
        } else {
            // Boolean attribute (present = true)
            attributes.values[name_token.value] = "true";
        }
    }
    
    return attributes;
}

BubblePtr HSMLParser::convert_bubble_element(std::shared_ptr<HSMLElement> element, BubblePtr parent) {
    if (element->tag_name != "bubble") {
        report_error("Expected 'bubble' element, got '" + element->tag_name + "'");
        return nullptr;
    }
    
    // Parse coordinates
    SphericalCoords coords = parse_coordinates(element->attributes);
    
    // Create bubble
    BubblePtr bubble = std::make_shared<Bubble>(coords, parent);
    
    // Set attributes
    for (const auto& attr : element->attributes.values) {
        bubble->set_attribute(attr.first, attr.second);
    }
    
    // Process children
    for (const auto& child : element->children) {
        if (child->tag_name == "bubble") {
            auto child_bubble = convert_bubble_element(child, bubble);
            if (child_bubble) {
                bubble->add_child(child_bubble);
            }
        } else if (child->tag_name == "presence") {
            auto presence = convert_presence_element(child);
            if (presence) {
                bubble->set_presence(presence);
                presence->set_bubble(bubble);
            }
        }
    }
    
    return bubble;
}

PresencePtr HSMLParser::convert_presence_element(std::shared_ptr<HSMLElement> element) {
    if (element->tag_name != "presence") {
        report_error("Expected 'presence' element, got '" + element->tag_name + "'");
        return nullptr;
    }
    
    // Create basic presence
    auto presence = std::make_shared<Presence>();
    
    // Set density
    double density = element->attributes.get_double("density", 1.0);
    presence->set_density(density);
    
    // Set material if specified
    std::string material_name = element->attributes.get("material", "default");
    std::string state_str = element->attributes.get("state", "solid");
    MaterialState state = parse_material_state(state_str);
    
    Material material(material_name, state, density);
    presence->set_material(material);
    
    // Set visibility
    presence->set_visible(element->attributes.get_bool("visible", true));
    
    return presence;
}

SphericalCoords HSMLParser::parse_coordinates(const AttributeMap& attributes) {
    double r = attributes.get_double("r", 0.0);
    double theta = attributes.get_double("theta", 0.0);
    double phi = attributes.get_double("phi", 0.0);
    
    // Convert degrees to radians if specified
    std::string angle_unit = attributes.get("angle_unit", "radians");
    if (angle_unit == "degrees") {
        theta = theta * SphericalCoords::PI / 180.0;
        phi = phi * SphericalCoords::PI / 180.0;
    }
    
    return SphericalCoords(r, theta, phi);
}

MaterialState HSMLParser::parse_material_state(const std::string& state_string) {
    std::string lower_state = state_string;
    std::transform(lower_state.begin(), lower_state.end(), lower_state.begin(), ::tolower);
    
    if (lower_state == "solid") return MaterialState::SOLID;
    if (lower_state == "liquid") return MaterialState::LIQUID;
    if (lower_state == "gas") return MaterialState::GAS;
    if (lower_state == "plasma") return MaterialState::PLASMA;
    
    return MaterialState::UNKNOWN;
}

void HSMLParser::validate_element(std::shared_ptr<HSMLElement> element) {
    if (element->tag_name == "bubble") {
        validate_bubble_attributes(element->attributes);
    } else if (element->tag_name == "presence") {
        validate_presence_attributes(element->attributes);
    }
}

void HSMLParser::validate_bubble_attributes(const AttributeMap& attributes) {
    // Check required attributes
    if (!attributes.has("r")) {
        report_error("Bubble element missing required 'r' attribute");
    }
    
    // Validate ranges
    double r = attributes.get_double("r", 0.0);
    if (r < 0.0) {
        report_error("Bubble radius 'r' must be non-negative");
    }
    
    double theta = attributes.get_double("theta", 0.0);
    if (theta < 0.0 || theta > SphericalCoords::PI) {
        report_error("Bubble theta must be between 0 and π");
    }
}

void HSMLParser::validate_presence_attributes(const AttributeMap& attributes) {
    double density = attributes.get_double("density", 1.0);
    if (density <= 0.0) {
        report_error("Presence density must be positive");
    }
}

void HSMLParser::report_error(const std::string& message) {
    errors_.emplace_back(message, lexer_.current_line(), lexer_.current_column());
}

void HSMLParser::report_error(const std::string& message, size_t line, size_t column) {
    errors_.emplace_back(message, line, column);
}

// Factory methods
namespace factory {

std::unique_ptr<HSMLParser> create_strict_parser(double viewer_distance) {
    auto parser = std::make_unique<HSMLParser>();
    parser->set_strict_mode(true);
    parser->set_default_viewer_distance(viewer_distance);
    return parser;
}

std::unique_ptr<HSMLParser> create_lenient_parser(double viewer_distance) {
    auto parser = std::make_unique<HSMLParser>();
    parser->set_strict_mode(false);
    parser->set_default_viewer_distance(viewer_distance);
    return parser;
}

std::unique_ptr<HSMLParser> create_debug_parser(double viewer_distance) {
    return create_strict_parser(viewer_distance);
}

} // namespace factory

// Utility functions
namespace utils {

std::string bubble_to_hsml(BubblePtr bubble, int indent_level) {
    if (!bubble) return "";
    
    std::string indent(indent_level * 2, ' ');
    std::stringstream ss;
    
    const auto& coords = bubble->coordinates();
    ss << indent << "<bubble r=\"" << coords.radius() 
       << "\" theta=\"" << coords.theta() 
       << "\" phi=\"" << coords.phi() << "\">\n";
    
    // Add presence if exists
    if (auto presence = bubble->presence()) {
        ss << presence_to_hsml(presence, indent_level + 1);
    }
    
    // Add children
    for (const auto& child : bubble->children()) {
        ss << bubble_to_hsml(child, indent_level + 1);
    }
    
    ss << indent << "</bubble>\n";
    return ss.str();
}

std::string presence_to_hsml(PresencePtr presence, int indent_level) {
    if (!presence) return "";
    
    std::string indent(indent_level * 2, ' ');
    std::stringstream ss;
    
    ss << indent << "<presence density=\"" << presence->density() << "\"";
    
    if (!presence->material().name().empty() && presence->material().name() != "default") {
        ss << " material=\"" << presence->material().name() << "\"";
    }
    
    ss << " visible=\"" << (presence->is_visible() ? "true" : "false") << "\"";
    ss << "/>\n";
    
    return ss.str();
}

bool validate_hsml_syntax(const std::string& hsml_content) {
    HSMLParser parser;
    auto result = parser.parse(hsml_content);
    return result != nullptr && !parser.has_errors();
}

std::vector<std::string> extract_hsml_errors(const std::string& hsml_content) {
    HSMLParser parser;
    parser.parse(hsml_content);
    
    std::vector<std::string> error_strings;
    for (const auto& error : parser.errors()) {
        error_strings.push_back(error.to_string());
    }
    
    return error_strings;
}

} // namespace utils

} // namespace parsers
} // namespace hsml