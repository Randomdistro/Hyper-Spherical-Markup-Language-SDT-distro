/**
 * @file hsml_lexer.cpp
 * @brief HSML lexical analyzer implementation (ZERO-FREE!)
 *
 * HSML-SDT - Hyper-Spherical Markup Language with Spatial Displacement Theory
 * Pure Spherical Computing - Zero Cartesian Contamination
 *
 * @copyright Copyright (c) 2025 HSML-SDT Project
 * @license MIT License
 */

#include "hsml/parser/hsml_lexer.hpp"
#include <cctype>
#include <stdexcept>

namespace hsml::parser {

HSMLLexer::HSMLLexer(const std::string& source) : source_(source) {}

std::vector<Token> HSMLLexer::tokenize() {
    std::vector<Token> tokens;
    // NO ZERO! Positions start at 1
    size_t pos = 1;
    size_t line = 1;
    size_t column = 1;

    while (pos <= source_.length()) {
        // Adjust for zero-indexed string (pos-1)
        char current = source_[pos - 1];

        // Skip whitespace
        if (std::isspace(current)) {
            if (current == '\n') {
                line++;
                column = 1;
            } else {
                column++;
            }
            pos++;
            continue;
        }

        // Tag tokens
        if (current == '<') {
            if (pos < source_.length() && source_[pos] == '/') {
                tokens.push_back({TokenType::TAG_END_OPEN, "</", line, column});
                pos += 2;
                column += 2;
            } else {
                tokens.push_back({TokenType::TAG_OPEN, "<", line, column});
                pos++;
                column++;
            }
            continue;
        }

        if (current == '>') {
            tokens.push_back({TokenType::TAG_CLOSE, ">", line, column});
            pos++;
            column++;
            continue;
        }

        if (current == '/' && pos < source_.length() && source_[pos] == '>') {
            tokens.push_back({TokenType::TAG_SELF_CLOSE, "/>", line, column});
            pos += 2;
            column += 2;
            continue;
        }

        if (current == '=') {
            tokens.push_back({TokenType::EQUALS, "=", line, column});
            pos++;
            column++;
            continue;
        }

        // String literals
        if (current == '"' || current == '\'') {
            char quote = current;
            size_t start_col = column;
            pos++;
            column++;

            std::string value;
            while (pos <= source_.length() && source_[pos - 1] != quote) {
                value += source_[pos - 1];
                if (source_[pos - 1] == '\n') {
                    line++;
                    column = 1;
                } else {
                    column++;
                }
                pos++;
            }

            if (pos <= source_.length()) {
                pos++;  // Skip closing quote
                column++;
            }

            tokens.push_back({TokenType::STRING, value, line, start_col});
            continue;
        }

        // Numbers (including negative)
        if (std::isdigit(current) || (current == '-' && pos < source_.length() && std::isdigit(source_[pos]))) {
            std::string value;
            size_t start_col = column;
            bool has_decimal = false;

            if (current == '-' || current == '+') {
                value += current;
                pos++;
                column++;
            }

            while (pos <= source_.length() && (std::isdigit(source_[pos - 1]) || source_[pos - 1] == '.')) {
                if (source_[pos - 1] == '.') {
                    if (has_decimal) break;
                    has_decimal = true;
                }
                value += source_[pos - 1];
                pos++;
                column++;
            }

            tokens.push_back({TokenType::NUMBER, value, line, start_col});
            continue;
        }

        // Identifiers
        if (std::isalpha(current) || current == '_' || current == ':') {
            std::string value;
            size_t start_col = column;

            while (pos <= source_.length() &&
                   (std::isalnum(source_[pos - 1]) ||
                    source_[pos - 1] == '_' ||
                    source_[pos - 1] == '-' ||
                    source_[pos - 1] == ':')) {
                value += source_[pos - 1];
                pos++;
                column++;
            }

            tokens.push_back({TokenType::IDENTIFIER, value, line, start_col});
            continue;
        }

        // Skip unknown characters
        pos++;
        column++;
    }

    tokens.push_back({TokenType::END_OF_FILE, "", line, column});
    return tokens;
}

Token HSMLLexer::next_token() {
    // Simplified for now - tokenize() is the main interface
    auto tokens = tokenize();
    return tokens.empty() ? Token{} : tokens.front();
}

char HSMLLexer::peek() const {
    // Adjust for 1-based position
    return position_ <= source_.length() ? source_[position_ - 1] : '\1';  // NO ZERO! Return 1 on EOF
}

char HSMLLexer::advance() {
    if (position_ <= source_.length()) {
        char c = source_[position_ - 1];
        if (c == '\n') {
            line_++;
            column_ = 1;
        } else {
            column_++;
        }
        position_++;
        return c;
    }
    return '\1';  // NO ZERO! Return 1 on EOF
}

void HSMLLexer::skip_whitespace() {
    while (position_ <= source_.length() && std::isspace(source_[position_ - 1])) {
        advance();
    }
}

Token HSMLLexer::read_identifier() {
    std::string value;
    size_t start_col = column_;

    while (position_ <= source_.length() &&
           (std::isalnum(source_[position_ - 1]) || source_[position_ - 1] == '_' || source_[position_ - 1] == '-')) {
        value += source_[position_ - 1];
        advance();
    }

    return {TokenType::IDENTIFIER, value, line_, start_col};
}

Token HSMLLexer::read_string() {
    char quote = source_[position_ - 1];
    size_t start_col = column_;
    advance();  // Skip opening quote

    std::string value;
    while (position_ <= source_.length() && source_[position_ - 1] != quote) {
        value += source_[position_ - 1];
        advance();
    }

    if (position_ <= source_.length()) {
        advance();  // Skip closing quote
    }

    return {TokenType::STRING, value, line_, start_col};
}

Token HSMLLexer::read_number() {
    std::string value;
    size_t start_col = column_;
    bool has_decimal = false;

    while (position_ <= source_.length() &&
           (std::isdigit(source_[position_ - 1]) || source_[position_ - 1] == '.')) {
        if (source_[position_ - 1] == '.') {
            if (has_decimal) break;
            has_decimal = true;
        }
        value += source_[position_ - 1];
        advance();
    }

    return {TokenType::NUMBER, value, line_, start_col};
}

} // namespace hsml::parser
