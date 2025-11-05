/**
 * @file tokens.hpp
 * @brief HSML token definitions for lexical analysis
 *
 * HSML-SDT - Hyper-Spherical Markup Language with Spatial Displacement Theory
 * Pure Spherical Computing - Zero Cartesian Contamination
 *
 * @copyright Copyright (c) 2025 HSML-SDT Project
 * @license MIT License
 */

#pragma once

#include <string>
#include <cstddef>

namespace hsml::parser {

/**
 * @brief Token types for HSML markup language
 */
enum class TokenType {
    TAG_OPEN,       ///< < - Opening tag bracket
    TAG_CLOSE,      ///< > - Closing tag bracket
    TAG_SELF_CLOSE, ///< /> - Self-closing tag
    TAG_END_OPEN,   ///< </ - End tag opening
    EQUALS,         ///< = - Attribute assignment
    STRING,         ///< "value" - String literal
    IDENTIFIER,     ///< element-name - Identifier
    NUMBER,         ///< 123.45 - Numeric literal
    WHITESPACE,     ///< Whitespace (usually skipped)
    END_OF_FILE     ///< End of input
};

/**
 * @brief Token structure containing type, value, and location
 */
struct Token {
    TokenType type;        ///< Type of token
    std::string value;     ///< Lexeme/value of token
    size_t line;           ///< Line number in source
    size_t column;         ///< Column number in source

    /**
     * @brief Construct a token
     */
    Token(TokenType t = TokenType::END_OF_FILE,
          std::string v = "",
          size_t l = 0,
          size_t c = 0)
        : type(t), value(std::move(v)), line(l), column(c) {}
};

/**
 * @brief Convert token type to string for debugging
 * @param type Token type to convert
 * @return String representation of token type
 */
inline const char* token_type_to_string(TokenType type) {
    switch (type) {
        case TokenType::TAG_OPEN: return "TAG_OPEN";
        case TokenType::TAG_CLOSE: return "TAG_CLOSE";
        case TokenType::TAG_SELF_CLOSE: return "TAG_SELF_CLOSE";
        case TokenType::TAG_END_OPEN: return "TAG_END_OPEN";
        case TokenType::EQUALS: return "EQUALS";
        case TokenType::STRING: return "STRING";
        case TokenType::IDENTIFIER: return "IDENTIFIER";
        case TokenType::NUMBER: return "NUMBER";
        case TokenType::WHITESPACE: return "WHITESPACE";
        case TokenType::END_OF_FILE: return "END_OF_FILE";
        default: return "UNKNOWN";
    }
}

} // namespace hsml::parser
