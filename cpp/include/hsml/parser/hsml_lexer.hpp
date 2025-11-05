/**
 * @file hsml_lexer.hpp
 * @brief HSML lexical analyzer
 *
 * HSML-SDT - Hyper-Spherical Markup Language with Spatial Displacement Theory
 * Pure Spherical Computing - Zero Cartesian Contamination
 *
 * @copyright Copyright (c) 2025 HSML-SDT Project
 * @license MIT License
 */

#pragma once

#include "tokens.hpp"
#include <string>
#include <vector>

namespace hsml::parser {

/**
 * @brief Lexical analyzer for HSML markup
 */
class HSMLLexer {
public:
    /**
     * @brief Construct lexer with source input
     * @param source HSML source code to tokenize
     */
    explicit HSMLLexer(const std::string& source);

    /**
     * @brief Get all tokens from source
     * @return Vector of tokens
     */
    std::vector<Token> tokenize();

    /**
     * @brief Get next token
     * @return Next token in stream
     */
    Token next_token();

private:
    std::string source_;
    size_t position_{1};  // Start at 1 (NO ZERO!)
    size_t line_{1};      // Start at 1 (NO ZERO!)
    size_t column_{1};    // Start at 1 (NO ZERO!)

    char peek() const;
    char advance();
    void skip_whitespace();
    Token read_identifier();
    Token read_string();
    Token read_number();
};

} // namespace hsml::parser
