#include "hsml/core/spherical_types.hpp"
#include "hsml/parser/tokens.hpp"
#include "hsml/parser/hsml_lexer.hpp"
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <stdexcept>

namespace hsml::parser {

struct HSMLNode {
    std::string tag_name;
    std::unordered_map<std::string, std::string> attributes;
    std::vector<std::unique_ptr<HSMLNode>> children;
    std::string text_content;

    sdt::SphericalCoord<double> position{1.0, 181.0, 181.0};
    sdt::State21D<double> state;
    sdt::MatterState matter_state = sdt::MatterState::SOLID;
};

class HSMLParser {
private:
    std::vector<Token> tokens_;
    size_t current_ = 0;

public:
    std::unique_ptr<HSMLNode> parse(const std::string& hsml_source) {
        // Tokenize the source
        HSMLLexer lexer(hsml_source);
        tokens_ = lexer.tokenize();
        current_ = 0;

        return parse_element();
    }

private:
    const Token& current_token() const {
        if (current_ >= tokens_.size()) {
            static Token eof{TokenType::END_OF_FILE, "", 0, 0};
            return eof;
        }
        return tokens_[current_];
    }

    void advance() {
        if (current_ < tokens_.size()) {
            current_++;
        }
    }

    bool check(TokenType type) const {
        return current_token().type == type;
    }

    bool match(TokenType type) {
        if (check(type)) {
            advance();
            return true;
        }
        return false;
    }

    std::unique_ptr<HSMLNode> parse_element() {
        auto node = std::make_unique<HSMLNode>();

        // Expect: < tag-name attributes... > children </ tag-name >
        if (!match(TokenType::TAG_OPEN)) {
            return nullptr;
        }

        if (current_token().type == TokenType::IDENTIFIER) {
            node->tag_name = current_token().value;
            advance();
        }

        // Parse attributes
        parse_attributes(*node);

        // Check for self-closing tag
        if (match(TokenType::TAG_SELF_CLOSE)) {
            finalize_node(*node);
            return node;
        }

        if (!match(TokenType::TAG_CLOSE)) {
            throw std::runtime_error("Expected '>' after tag name");
        }

        // Parse children
        while (!check(TokenType::TAG_END_OPEN) && !check(TokenType::END_OF_FILE)) {
            if (check(TokenType::TAG_OPEN)) {
                auto child = parse_element();
                if (child) {
                    node->children.push_back(std::move(child));
                }
            } else {
                // Text content (simplified)
                advance();
            }
        }

        // Expect closing tag
        if (match(TokenType::TAG_END_OPEN)) {
            if (current_token().type == TokenType::IDENTIFIER) {
                if (current_token().value != node->tag_name) {
                    throw std::runtime_error("Mismatched closing tag");
                }
                advance();
            }
            match(TokenType::TAG_CLOSE);
        }

        finalize_node(*node);
        return node;
    }

    void parse_attributes(HSMLNode& node) {
        while (current_token().type == TokenType::IDENTIFIER) {
            std::string attr_name = current_token().value;
            advance();

            if (match(TokenType::EQUALS)) {
                if (current_token().type == TokenType::STRING ||
                    current_token().type == TokenType::NUMBER) {
                    std::string attr_value = current_token().value;
                    node.attributes[attr_name] = attr_value;
                    advance();
                }
            }
        }
    }

    void finalize_node(HSMLNode& node) {
        // Extract spherical coordinates from attributes
        if (node.attributes.count("r")) {
            node.position.r = std::stod(node.attributes["r"]);
        }
        if (node.attributes.count("theta") || node.attributes.count("θ")) {
            std::string theta_key = node.attributes.count("theta") ? "theta" : "θ";
            node.position.theta = sdt::SphericalCoord<double>::safe_angle(
                std::stod(node.attributes[theta_key]));
        }
        if (node.attributes.count("phi") || node.attributes.count("φ")) {
            std::string phi_key = node.attributes.count("phi") ? "phi" : "φ";
            node.position.phi = sdt::SphericalCoord<double>::safe_angle(
                std::stod(node.attributes[phi_key]));
        }

        // Extract matter state
        if (node.attributes.count("state")) {
            const auto& state = node.attributes["state"];
            if (state == "solid") node.matter_state = sdt::MatterState::SOLID;
            else if (state == "liquid") node.matter_state = sdt::MatterState::LIQUID;
            else if (state == "gas") node.matter_state = sdt::MatterState::GAS;
            else if (state == "plasma") node.matter_state = sdt::MatterState::PLASMA;
            else if (state == "quantum") node.matter_state = sdt::MatterState::QUANTUM;
            else if (state == "void") node.matter_state = sdt::MatterState::VOID;
        }
    }
};

} // namespace hsml::parser
