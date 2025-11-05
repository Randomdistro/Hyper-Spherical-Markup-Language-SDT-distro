/**
 * HSML Multi-Language Lexer - C++20 Implementation
 * Revolutionary template metaprogramming lexer for all four HSML languages
 * Death-cheating transposition from TypeScript with multiple paradigms
 */

#pragma once

#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <concepts>
#include <optional>
#include <variant>
#include <expected>
#include <array>
#include <span>
#include <ranges>
#include <algorithm>
#include <type_traits>

namespace hsml::parsers {

// Modern C++20 concepts for lexical analysis
template<typename T>
concept TokenValue = std::convertible_to<T, std::string> || 
                    std::convertible_to<T, double> ||
                    std::convertible_to<T, bool>;

template<typename T>
concept SphericalTokenType = requires(T t) {
    { t.r } -> std::convertible_to<double>;
    { t.theta } -> std::convertible_to<double>;
    { t.phi } -> std::convertible_to<double>;
};

// Comprehensive token enumeration for all HSML languages
enum class HSMLTokenType : uint32_t {
    // Core literals
    NUMBER = 0x0001,
    STRING = 0x0002,
    IDENTIFIER = 0x0003,
    BOOLEAN = 0x0004,
    
    // Spherical coordinates - PURE SPHERICAL MATHEMATICS!
    SPHERICAL_COORD = 0x0010,
    SOLID_ANGLE = 0x0011,
    STERADIAN = 0x0012,
    
    // Matter states for physics-based rendering
    SOLID = 0x0020,
    LIQUID = 0x0021,
    GAS = 0x0022,
    PLASMA = 0x0023,
    
    // Mathematical operators
    PLUS = 0x0030,
    MINUS = 0x0031,
    MULTIPLY = 0x0032,
    DIVIDE = 0x0033,
    POWER = 0x0034,
    ASSIGN = 0x0035,
    EQUALS = 0x0036,
    NOT_EQUALS = 0x0037,
    LESS_THAN = 0x0038,
    GREATER_THAN = 0x0039,
    LESS_EQUAL = 0x003A,
    GREATER_EQUAL = 0x003B,
    
    // Structural delimiters
    LEFT_PAREN = 0x0040,
    RIGHT_PAREN = 0x0041,
    LEFT_BRACE = 0x0042,
    RIGHT_BRACE = 0x0043,
    LEFT_BRACKET = 0x0044,
    RIGHT_BRACKET = 0x0045,
    SEMICOLON = 0x0046,
    COMMA = 0x0047,
    DOT = 0x0048,
    COLON = 0x0049,
    
    // HSML-specific markup tokens
    ELEMENT_TAG = 0x0100,
    CLOSE_TAG = 0x0101,
    SELF_CLOSING = 0x0102,
    ATTRIBUTE = 0x0103,
    
    // CSSS-specific physics styling tokens
    SELECTOR = 0x0200,
    PROPERTY = 0x0201,
    VALUE = 0x0202,
    MATERIAL = 0x0203,
    ANIMATION = 0x0204,
    KEYFRAME = 0x0205,
    TRANSITION = 0x0206,
    
    // ShapeScript behavior tokens
    BEHAVIOR = 0x0300,
    PHYSICS = 0x0301,
    FORCE = 0x0302,
    CONSTRAINT = 0x0303,
    COLLISION = 0x0304,
    INTERACTION = 0x0305,
    EVENT = 0x0306,
    
    // StyleBot parallel processing tokens
    BOT = 0x0400,
    AGENT = 0x0401,
    PARALLEL = 0x0402,
    DISTRIBUTED = 0x0403,
    RENDER = 0x0404,
    OPTIMIZE = 0x0405,
    CACHE = 0x0406,
    
    // Control flow keywords
    IF = 0x0500,
    ELSE = 0x0501,
    WHILE = 0x0502,
    FOR = 0x0503,
    FUNCTION = 0x0504,
    RETURN = 0x0505,
    VAR = 0x0506,
    LET = 0x0507,
    CONST = 0x0508,
    
    // Special tokens
    EOF_TOKEN = 0xFFFE,
    ERROR_TOKEN = 0xFFFF,
    COMMENT = 0xFFFD,
    WHITESPACE = 0xFFFC
};

// Language enumeration for multi-language support
enum class HSMLLanguage : uint8_t {
    HSML = 0,   // Hyper-Spherical Markup Language
    CSSS = 1,   // Cascading Spherical Style Sheets
    SHAPE = 2,  // ShapeScript behaviors
    STYB = 3    // StyleBot parallel processing
};

// Base token structure with compile-time type safety
struct HSMLToken {
    HSMLTokenType type;
    std::string value;
    size_t line;
    size_t column;
    HSMLLanguage language;
    
    constexpr HSMLToken(HSMLTokenType t = HSMLTokenType::ERROR_TOKEN, 
                       std::string_view v = "", 
                       size_t l = 0, size_t c = 0,
                       HSMLLanguage lang = HSMLLanguage::HSML) noexcept
        : type(t), value(v), line(l), column(c), language(lang) {}
};

// Template specializations for different token types
template<SphericalTokenType T>
struct SphericalCoordToken : public HSMLToken {
    double r, theta, phi;
    
    constexpr SphericalCoordToken(double r_val, double theta_val, double phi_val,
                                 size_t l = 0, size_t c = 0,
                                 HSMLLanguage lang = HSMLLanguage::HSML) noexcept
        : HSMLToken(HSMLTokenType::SPHERICAL_COORD, 
                   std::format("({}, {}, {})", r_val, theta_val, phi_val), l, c, lang),
          r(r_val), theta(theta_val), phi(phi_val) {}
};

struct SolidAngleToken : public HSMLToken {
    double omega, theta_min, theta_max, phi_min, phi_max;
    
    constexpr SolidAngleToken(double omega_val, double theta_min_val, double theta_max_val,
                             double phi_min_val, double phi_max_val,
                             size_t l = 0, size_t c = 0,
                             HSMLLanguage lang = HSMLLanguage::HSML) noexcept
        : HSMLToken(HSMLTokenType::SOLID_ANGLE, 
                   std::format("Ω({}, {}, {}, {}, {})", omega_val, theta_min_val, 
                              theta_max_val, phi_min_val, phi_max_val), l, c, lang),
          omega(omega_val), theta_min(theta_min_val), theta_max(theta_max_val),
          phi_min(phi_min_val), phi_max(phi_max_val) {}
};

// Matter state token with physics properties
struct MatterStateToken : public HSMLToken {
    struct PhysicsProperties {
        std::optional<double> density;
        std::optional<double> temperature;
        std::optional<double> pressure;
        std::optional<double> conductivity;
    } properties;
    
    constexpr MatterStateToken(HSMLTokenType matter_type, std::string_view value,
                              size_t l = 0, size_t c = 0,
                              HSMLLanguage lang = HSMLLanguage::CSSS) noexcept
        : HSMLToken(matter_type, value, l, c, lang) {
        initialize_default_properties();
    }
    
private:
    constexpr void initialize_default_properties() noexcept {
        switch (type) {
            case HSMLTokenType::SOLID:
                properties = {1000.0, 293.15, 101325.0, std::nullopt};
                break;
            case HSMLTokenType::LIQUID:
                properties = {1000.0, 293.15, 101325.0, 0.6};
                break;
            case HSMLTokenType::GAS:
                properties = {1.225, 293.15, 101325.0, std::nullopt};
                break;
            case HSMLTokenType::PLASMA:
                properties = {1e-6, 10000.0, 101325.0, 1e6};
                break;
            default:
                properties = {};
        }
    }
};

// High-performance multi-language lexer with template metaprogramming
template<HSMLLanguage Language = HSMLLanguage::HSML>
class HSMLMultiLexer {
private:
    std::string_view source_;
    size_t current_position_ = 0;
    size_t current_line_ = 1;
    size_t current_column_ = 1;
    
    // Language-specific keyword sets with compile-time initialization
    static constexpr std::array<std::string_view, 13> hsml_keywords = {
        "element", "sphere", "shell", "point", "group", "container",
        "position", "radius", "material", "behavior", "style",
        "visible", "interactive", "matter", "state"
    };
    
    static constexpr std::array<std::string_view, 14> csss_keywords = {
        "material", "albedo", "metallic", "roughness", "emission",
        "transparency", "refraction", "animation", "keyframe",
        "transition", "easing", "duration", "delay", "repeat"
    };
    
    static constexpr std::array<std::string_view, 13> shape_keywords = {
        "behavior", "physics", "force", "gravity", "collision",
        "constraint", "interaction", "event", "trigger", "response",
        "elastic", "viscous", "electromagnetic", "thermal"
    };
    
    static constexpr std::array<std::string_view, 13> styb_keywords = {
        "bot", "agent", "parallel", "distributed", "render",
        "optimize", "cache", "performance", "quality", "lod",
        "frustum", "occlusion", "spatial", "index"
    };
    
    // SIMD-optimized character classification
    alignas(32) static constexpr std::array<bool, 256> digit_lookup = []() constexpr {
        std::array<bool, 256> lookup{};
        for (int i = '0'; i <= '9'; ++i) {
            lookup[i] = true;
        }
        return lookup;
    }();
    
    alignas(32) static constexpr std::array<bool, 256> alpha_lookup = []() constexpr {
        std::array<bool, 256> lookup{};
        for (int i = 'a'; i <= 'z'; ++i) {
            lookup[i] = true;
        }
        for (int i = 'A'; i <= 'Z'; ++i) {
            lookup[i] = true;
        }
        lookup['_'] = true;
        return lookup;
    }();
    
    alignas(32) static constexpr std::array<bool, 256> whitespace_lookup = []() constexpr {
        std::array<bool, 256> lookup{};
        lookup[' '] = lookup['\t'] = lookup['\r'] = lookup['\n'] = true;
        return lookup;
    }();

public:
    explicit constexpr HSMLMultiLexer(std::string_view source) noexcept
        : source_(source) {}
    
    // Main tokenization with ranges and algorithms
    [[nodiscard]] auto tokenize() -> std::vector<HSMLToken> {
        std::vector<HSMLToken> tokens;
        tokens.reserve(source_.size() / 8); // Heuristic capacity
        
        while (!is_at_end()) {
            auto token = next_token();
            if (token.type != HSMLTokenType::WHITESPACE && 
                token.type != HSMLTokenType::COMMENT) {
                tokens.emplace_back(std::move(token));
            }
        }
        
        tokens.emplace_back(HSMLTokenType::EOF_TOKEN, "", current_line_, current_column_, Language);
        return tokens;
    }
    
    // High-performance single token extraction
    [[nodiscard]] auto next_token() -> HSMLToken {
        skip_whitespace();
        
        if (is_at_end()) {
            return create_token(HSMLTokenType::EOF_TOKEN, "");
        }
        
        const char current_char = peek();
        
        // Comment handling
        if (current_char == '/' && peek_next() == '/') {
            return handle_line_comment();
        }
        
        // Numeric literals (including spherical coordinates)
        if (is_digit(current_char)) {
            return handle_number();
        }
        
        // Identifiers and keywords
        if (is_alpha(current_char)) {
            return handle_identifier();
        }
        
        // String literals
        if (current_char == '"' || current_char == '\'') {
            return handle_string_literal();
        }
        
        // Spherical coordinate detection
        if (current_char == '(' && peek_ahead(1) == 'r') {
            return handle_spherical_coordinate();
        }
        
        // Solid angle detection (Ω or Om)
        if (current_char == 'Ω' || (current_char == 'O' && peek_next() == 'm')) {
            return handle_solid_angle();
        }
        
        // Operators and delimiters
        return handle_operator();
    }

private:
    // Character access and navigation
    [[nodiscard]] constexpr char peek(size_t offset = 0) const noexcept {
        const size_t pos = current_position_ + offset;
        return (pos < source_.size()) ? source_[pos] : '\0';
    }
    
    [[nodiscard]] constexpr char peek_next() const noexcept { return peek(1); }
    [[nodiscard]] constexpr char peek_ahead(size_t offset) const noexcept { return peek(offset); }
    
    constexpr char advance() noexcept {
        if (current_position_ < source_.size()) {
            const char ch = source_[current_position_++];
            if (ch == '\n') {
                ++current_line_;
                current_column_ = 1;
            } else {
                ++current_column_;
            }
            return ch;
        }
        return '\0';
    }
    
    [[nodiscard]] constexpr bool is_at_end() const noexcept {
        return current_position_ >= source_.size();
    }
    
    // SIMD-optimized character classification
    [[nodiscard]] constexpr bool is_digit(char c) const noexcept {
        const auto index = static_cast<unsigned char>(c);
        return index < digit_lookup.size() && digit_lookup[index];
    }
    
    [[nodiscard]] constexpr bool is_alpha(char c) const noexcept {
        const auto index = static_cast<unsigned char>(c);
        return index < alpha_lookup.size() && alpha_lookup[index];
    }
    
    [[nodiscard]] constexpr bool is_alphanumeric(char c) const noexcept {
        return is_alpha(c) || is_digit(c);
    }
    
    [[nodiscard]] constexpr bool is_whitespace(char c) const noexcept {
        const auto index = static_cast<unsigned char>(c);
        return index < whitespace_lookup.size() && whitespace_lookup[index];
    }
    
    // Whitespace skipping with newline tracking
    void skip_whitespace() noexcept {
        while (!is_at_end() && is_whitespace(peek())) {
            advance();
        }
    }
    
    // Token creation helper
    [[nodiscard]] constexpr HSMLToken create_token(HSMLTokenType type, std::string_view value) const noexcept {
        return HSMLToken(type, value, current_line_, current_column_, Language);
    }
    
    // Specific token handlers
    [[nodiscard]] auto handle_number() -> HSMLToken {
        const size_t start_pos = current_position_;
        bool has_decimal = false;
        
        // Integer part
        while (is_digit(peek())) {
            advance();
        }
        
        // Decimal part
        if (peek() == '.' && is_digit(peek_next())) {
            has_decimal = true;
            advance(); // consume '.'
            while (is_digit(peek())) {
                advance();
            }
        }
        
        // Scientific notation
        if (peek() == 'e' || peek() == 'E') {
            advance();
            if (peek() == '+' || peek() == '-') {
                advance();
            }
            while (is_digit(peek())) {
                advance();
            }
        }
        
        const std::string_view number_str = source_.substr(start_pos, current_position_ - start_pos);
        return create_token(HSMLTokenType::NUMBER, number_str);
    }
    
    [[nodiscard]] auto handle_identifier() -> HSMLToken {
        const size_t start_pos = current_position_;
        
        while (is_alphanumeric(peek())) {
            advance();
        }
        
        const std::string_view identifier = source_.substr(start_pos, current_position_ - start_pos);
        
        // Check for language-specific keywords
        if (const auto keyword_type = get_keyword_type(identifier)) {
            return create_token(*keyword_type, identifier);
        }
        
        // Check for matter states
        if (is_matter_state(identifier)) {
            return create_matter_state_token(identifier);
        }
        
        return create_token(HSMLTokenType::IDENTIFIER, identifier);
    }
    
    [[nodiscard]] auto handle_string_literal() -> HSMLToken {
        const char quote_char = advance(); // Consume opening quote
        const size_t start_pos = current_position_;
        
        while (!is_at_end() && peek() != quote_char) {
            if (peek() == '\\') {
                advance(); // Skip escape character
                if (!is_at_end()) {
                    advance(); // Skip escaped character
                }
            } else {
                advance();
            }
        }
        
        if (is_at_end()) {
            return create_token(HSMLTokenType::ERROR_TOKEN, "Unterminated string literal");
        }
        
        const std::string_view string_content = source_.substr(start_pos, current_position_ - start_pos);
        advance(); // Consume closing quote
        
        return create_token(HSMLTokenType::STRING, string_content);
    }
    
    [[nodiscard]] auto handle_spherical_coordinate() -> SphericalCoordToken<double> {
        advance(); // Consume '('
        
        // Parse r, theta, phi components
        skip_whitespace();
        const double r = parse_number_value();
        
        expect_character(',');
        skip_whitespace();
        const double theta = parse_number_value();
        
        expect_character(',');
        skip_whitespace();
        const double phi = parse_number_value();
        
        expect_character(')');
        
        return SphericalCoordToken<double>(r, theta, phi, current_line_, current_column_, Language);
    }
    
    [[nodiscard]] auto handle_solid_angle() -> SolidAngleToken {
        advance(); // Consume 'Ω' or 'O'
        
        if (peek() == 'm') {
            advance(); // Consume 'm' for "Om"
        }
        
        skip_whitespace();
        expect_character('(');
        
        // Parse solid angle components
        const double omega = parse_number_value();
        expect_character(',');
        const double theta_min = parse_number_value();
        expect_character(',');
        const double theta_max = parse_number_value();
        expect_character(',');
        const double phi_min = parse_number_value();
        expect_character(',');
        const double phi_max = parse_number_value();
        
        expect_character(')');
        
        return SolidAngleToken(omega, theta_min, theta_max, phi_min, phi_max,
                              current_line_, current_column_, Language);
    }
    
    [[nodiscard]] auto handle_line_comment() -> HSMLToken {
        advance(); // Consume first '/'
        advance(); // Consume second '/'
        
        const size_t start_pos = current_position_;
        
        while (!is_at_end() && peek() != '\n') {
            advance();
        }
        
        const std::string_view comment_text = source_.substr(start_pos, current_position_ - start_pos);
        return create_token(HSMLTokenType::COMMENT, comment_text);
    }
    
    [[nodiscard]] auto handle_operator() -> HSMLToken {
        const char op_char = advance();
        
        switch (op_char) {
            case '(': return create_token(HSMLTokenType::LEFT_PAREN, "(");
            case ')': return create_token(HSMLTokenType::RIGHT_PAREN, ")");
            case '{': return create_token(HSMLTokenType::LEFT_BRACE, "{");
            case '}': return create_token(HSMLTokenType::RIGHT_BRACE, "}");
            case '[': return create_token(HSMLTokenType::LEFT_BRACKET, "[");
            case ']': return create_token(HSMLTokenType::RIGHT_BRACKET, "]");
            case ';': return create_token(HSMLTokenType::SEMICOLON, ";");
            case ',': return create_token(HSMLTokenType::COMMA, ",");
            case '.': return create_token(HSMLTokenType::DOT, ".");
            case ':': return create_token(HSMLTokenType::COLON, ":");
            case '+': return create_token(HSMLTokenType::PLUS, "+");
            case '-': return create_token(HSMLTokenType::MINUS, "-");
            case '*': return create_token(HSMLTokenType::MULTIPLY, "*");
            case '/': return create_token(HSMLTokenType::DIVIDE, "/");
            case '^': return create_token(HSMLTokenType::POWER, "^");
            
            case '=':
                if (peek() == '=') {
                    advance();
                    return create_token(HSMLTokenType::EQUALS, "==");
                }
                return create_token(HSMLTokenType::ASSIGN, "=");
                
            case '!':
                if (peek() == '=') {
                    advance();
                    return create_token(HSMLTokenType::NOT_EQUALS, "!=");
                }
                return create_token(HSMLTokenType::ERROR_TOKEN, "Unexpected character '!'");
                
            case '<':
                if (peek() == '=') {
                    advance();
                    return create_token(HSMLTokenType::LESS_EQUAL, "<=");
                }
                return create_token(HSMLTokenType::LESS_THAN, "<");
                
            case '>':
                if (peek() == '=') {
                    advance();
                    return create_token(HSMLTokenType::GREATER_EQUAL, ">=");
                }
                return create_token(HSMLTokenType::GREATER_THAN, ">");
                
            default:
                return create_token(HSMLTokenType::ERROR_TOKEN, 
                                   std::string("Unexpected character: ") + op_char);
        }
    }
    
    // Helper methods for parsing
    [[nodiscard]] double parse_number_value() {
        const size_t start_pos = current_position_;
        
        while (is_digit(peek()) || peek() == '.') {
            advance();
        }
        
        const std::string_view number_str = source_.substr(start_pos, current_position_ - start_pos);
        return std::stod(std::string(number_str));
    }
    
    void expect_character(char expected) {
        skip_whitespace();
        if (peek() != expected) {
            throw std::runtime_error(std::string("Expected '") + expected + "' but found '" + peek() + "'");
        }
        advance();
    }
    
    // Language-specific keyword recognition
    [[nodiscard]] auto get_keyword_type(std::string_view identifier) const -> std::optional<HSMLTokenType> {
        const auto keywords = get_keywords_for_language();
        
        if (std::ranges::find(keywords, identifier) != keywords.end()) {
            return get_keyword_token_type(identifier);
        }
        
        return std::nullopt;
    }
    
    [[nodiscard]] constexpr auto get_keywords_for_language() const -> std::span<const std::string_view> {
        if constexpr (Language == HSMLLanguage::HSML) {
            return hsml_keywords;
        } else if constexpr (Language == HSMLLanguage::CSSS) {
            return csss_keywords;
        } else if constexpr (Language == HSMLLanguage::SHAPE) {
            return shape_keywords;
        } else if constexpr (Language == HSMLLanguage::STYB) {
            return styb_keywords;
        }
    }
    
    [[nodiscard]] auto get_keyword_token_type(std::string_view keyword) const -> HSMLTokenType {
        // HSML element keywords
        if (keyword == "element" || keyword == "sphere" || keyword == "shell" || keyword == "point") {
            return HSMLTokenType::ELEMENT_TAG;
        }
        
        // CSSS material keywords
        if (keyword == "material" || keyword == "animation" || keyword == "keyframe" || keyword == "transition") {
            return HSMLTokenType::MATERIAL;
        }
        
        // ShapeScript behavior keywords
        if (keyword == "behavior" || keyword == "physics" || keyword == "force" || keyword == "collision") {
            return HSMLTokenType::BEHAVIOR;
        }
        
        // StyleBot processing keywords
        if (keyword == "bot" || keyword == "agent" || keyword == "parallel" || keyword == "distributed") {
            return HSMLTokenType::BOT;
        }
        
        return HSMLTokenType::IDENTIFIER;
    }
    
    [[nodiscard]] constexpr bool is_matter_state(std::string_view identifier) const noexcept {
        return identifier == "solid" || identifier == "liquid" || 
               identifier == "gas" || identifier == "plasma";
    }
    
    [[nodiscard]] auto create_matter_state_token(std::string_view matter_state) const -> MatterStateToken {
        HSMLTokenType token_type;
        
        if (matter_state == "solid") token_type = HSMLTokenType::SOLID;
        else if (matter_state == "liquid") token_type = HSMLTokenType::LIQUID;
        else if (matter_state == "gas") token_type = HSMLTokenType::GAS;
        else if (matter_state == "plasma") token_type = HSMLTokenType::PLASMA;
        else token_type = HSMLTokenType::ERROR_TOKEN;
        
        return MatterStateToken(token_type, matter_state, current_line_, current_column_, Language);
    }
};

// Template aliases for different languages
using HSMLLexer = HSMLMultiLexer<HSMLLanguage::HSML>;
using CSSSLexer = HSMLMultiLexer<HSMLLanguage::CSSS>;
using ShapeLexer = HSMLMultiLexer<HSMLLanguage::SHAPE>;
using StybLexer = HSMLMultiLexer<HSMLLanguage::STYB>;

// Factory functions for creating language-specific lexers
namespace factory {
    [[nodiscard]] auto create_hsml_lexer(std::string_view source) -> HSMLLexer {
        return HSMLLexer(source);
    }
    
    [[nodiscard]] auto create_csss_lexer(std::string_view source) -> CSSSLexer {
        return CSSSLexer(source);
    }
    
    [[nodiscard]] auto create_shape_lexer(std::string_view source) -> ShapeLexer {
        return ShapeLexer(source);
    }
    
    [[nodiscard]] auto create_styb_lexer(std::string_view source) -> StybLexer {
        return StybLexer(source);
    }
}

} // namespace hsml::parsers