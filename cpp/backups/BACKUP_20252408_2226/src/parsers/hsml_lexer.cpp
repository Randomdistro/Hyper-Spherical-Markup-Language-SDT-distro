/**
 * HSML Multi-Language Lexer Implementation - C++20
 * Death-cheating transposition with SIMD optimizations and template instantiations
 */

#include "hsml/parsers/hsml_lexer.h"
#include <immintrin.h>
#include <execution>
#include <numeric>
#include <bit>

namespace hsml::parsers {

// Template instantiations for all supported languages
template class HSMLMultiLexer<HSMLLanguage::HSML>;
template class HSMLMultiLexer<HSMLLanguage::CSSS>;
template class HSMLMultiLexer<HSMLLanguage::SHAPE>;
template class HSMLMultiLexer<HSMLLanguage::STYB>;

// SIMD-optimized character classification using AVX2
namespace simd_optimized {
    
    // Vectorized whitespace detection
    [[nodiscard]] bool has_whitespace_simd(const char* text, size_t length) noexcept {
        if (length < 32) {
            // Fallback to scalar for small strings
            return std::any_of(text, text + length, [](char c) {
                return c == ' ' || c == '\t' || c == '\r' || c == '\n';
            });
        }
        
#ifdef __AVX2__
        const __m256i whitespace_chars = _mm256_set_epi8(
            '\n', '\r', '\t', ' ', '\n', '\r', '\t', ' ',
            '\n', '\r', '\t', ' ', '\n', '\r', '\t', ' ',
            '\n', '\r', '\t', ' ', '\n', '\r', '\t', ' ',
            '\n', '\r', '\t', ' ', '\n', '\r', '\t', ' '
        );
        
        for (size_t i = 0; i < length; i += 32) {
            const __m256i chunk = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(text + i)
            );
            
            const __m256i cmp_result = _mm256_cmpeq_epi8(chunk, whitespace_chars);
            
            if (_mm256_movemask_epi8(cmp_result) != 0) {
                return true;
            }
        }
        
        return false;
#else
        // Fallback to scalar implementation
        return std::any_of(text, text + length, [](char c) {
            return c == ' ' || c == '\t' || c == '\r' || c == '\n';
        });
#endif
    }
    
    // Vectorized digit detection
    [[nodiscard]] size_t count_digits_simd(const char* text, size_t length) noexcept {
#ifdef __AVX2__
        if (length < 32) {
            return std::count_if(text, text + length, [](char c) {
                return c >= '0' && c <= '9';
            });
        }
        
        const __m256i digit_min = _mm256_set1_epi8('0');
        const __m256i digit_max = _mm256_set1_epi8('9');
        size_t count = 0;
        
        for (size_t i = 0; i < length; i += 32) {
            const __m256i chunk = _mm256_loadu_si256(
                reinterpret_cast<const __m256i*>(text + i)
            );
            
            const __m256i ge_min = _mm256_cmpgt_epi8(chunk, _mm256_sub_epi8(digit_min, _mm256_set1_epi8(1)));
            const __m256i le_max = _mm256_cmpgt_epi8(_mm256_add_epi8(digit_max, _mm256_set1_epi8(1)), chunk);
            const __m256i is_digit = _mm256_and_si256(ge_min, le_max);
            
            count += std::popcount(static_cast<uint32_t>(_mm256_movemask_epi8(is_digit)));
        }
        
        return count;
#else
        return std::count_if(text, text + length, [](char c) {
            return c >= '0' && c <= '9';
        });
#endif
    }
}

// Optimized string parsing utilities
namespace parsing_utils {
    
    // Fast floating-point parser with scientific notation support
    [[nodiscard]] auto parse_double_fast(std::string_view str) -> std::expected<double, std::string> {
        if (str.empty()) {
            return std::unexpected("Empty string cannot be parsed as double");
        }
        
        double result = 0.0;
        bool negative = false;
        size_t pos = 0;
        
        // Handle sign
        if (str[pos] == '-') {
            negative = true;
            ++pos;
        } else if (str[pos] == '+') {
            ++pos;
        }
        
        // Parse integer part
        while (pos < str.size() && std::isdigit(str[pos])) {
            result = result * 10.0 + (str[pos] - '0');
            ++pos;
        }
        
        // Parse fractional part
        if (pos < str.size() && str[pos] == '.') {
            ++pos;
            double fractional = 0.0;
            double divisor = 10.0;
            
            while (pos < str.size() && std::isdigit(str[pos])) {
                fractional += (str[pos] - '0') / divisor;
                divisor *= 10.0;
                ++pos;
            }
            
            result += fractional;
        }
        
        // Parse scientific notation
        if (pos < str.size() && (str[pos] == 'e' || str[pos] == 'E')) {
            ++pos;
            bool exp_negative = false;
            int exponent = 0;
            
            if (pos < str.size() && str[pos] == '-') {
                exp_negative = true;
                ++pos;
            } else if (pos < str.size() && str[pos] == '+') {
                ++pos;
            }
            
            while (pos < str.size() && std::isdigit(str[pos])) {
                exponent = exponent * 10 + (str[pos] - '0');
                ++pos;
            }
            
            if (exp_negative) {
                exponent = -exponent;
            }
            
            result *= std::pow(10.0, exponent);
        }
        
        return negative ? -result : result;
    }
    
    // Unicode-aware spherical coordinate symbol detection
    [[nodiscard]] bool is_spherical_symbol(char32_t codepoint) noexcept {
        switch (codepoint) {
            case U'Ω':  // Greek capital omega
            case U'ω':  // Greek small omega
            case U'θ':  // Greek small theta
            case U'Θ':  // Greek capital theta
            case U'φ':  // Greek small phi
            case U'Φ':  // Greek capital phi
            case U'ρ':  // Greek small rho
            case U'Ρ':  // Greek capital rho
                return true;
            default:
                return false;
        }
    }
}

// Performance-optimized token buffer management
class TokenBuffer {
private:
    static constexpr size_t INITIAL_CAPACITY = 1024;
    static constexpr size_t GROWTH_FACTOR = 2;
    
    std::vector<HSMLToken> tokens_;
    size_t current_capacity_;
    
public:
    TokenBuffer() : current_capacity_(INITIAL_CAPACITY) {
        tokens_.reserve(current_capacity_);
    }
    
    void push_token(HSMLToken&& token) {
        if (tokens_.size() >= current_capacity_) {
            current_capacity_ *= GROWTH_FACTOR;
            tokens_.reserve(current_capacity_);
        }
        
        tokens_.emplace_back(std::move(token));
    }
    
    [[nodiscard]] std::vector<HSMLToken> extract_tokens() && {
        return std::move(tokens_);
    }
    
    [[nodiscard]] size_t size() const noexcept { return tokens_.size(); }
    [[nodiscard]] bool empty() const noexcept { return tokens_.empty(); }
};

// Advanced lexical analysis metrics
struct LexicalMetrics {
    size_t total_tokens = 0;
    size_t unique_identifiers = 0;
    size_t spherical_coordinates = 0;
    size_t solid_angles = 0;
    size_t matter_state_tokens = 0;
    double parse_time_microseconds = 0.0;
    
    void reset() noexcept {
        *this = LexicalMetrics{};
    }
    
    [[nodiscard]] double tokens_per_second() const noexcept {
        return parse_time_microseconds > 0.0 
            ? (total_tokens * 1'000'000.0 / parse_time_microseconds)
            : 0.0;
    }
};

// Global lexical analysis performance tracker
static thread_local LexicalMetrics current_metrics;

// Enhanced error reporting with source context
class LexicalError {
private:
    std::string message_;
    size_t line_;
    size_t column_;
    size_t position_;
    std::string_view source_context_;
    
public:
    LexicalError(std::string message, size_t line, size_t column, 
                size_t position, std::string_view context)
        : message_(std::move(message)), line_(line), column_(column),
          position_(position), source_context_(context) {}
    
    [[nodiscard]] std::string format_error() const {
        return std::format("Lexical Error at line {}, column {}: {}\n"
                          "Context: '{}'\n"
                          "Position: {}", 
                          line_, column_, message_, source_context_, position_);
    }
    
    [[nodiscard]] const std::string& message() const noexcept { return message_; }
    [[nodiscard]] size_t line() const noexcept { return line_; }
    [[nodiscard]] size_t column() const noexcept { return column_; }
};

// Multi-threaded lexical analysis for large files
namespace parallel_lexing {
    
    // Chunk-based parallel tokenization
    [[nodiscard]] auto tokenize_parallel(std::string_view source, size_t chunk_size = 64'000) 
        -> std::vector<HSMLToken> {
        
        if (source.size() < chunk_size * 2) {
            // For small inputs, use single-threaded approach
            HSMLLexer lexer(source);
            return lexer.tokenize();
        }
        
        const size_t num_chunks = (source.size() + chunk_size - 1) / chunk_size;
        std::vector<std::future<std::vector<HSMLToken>>> futures;
        futures.reserve(num_chunks);
        
        // Launch parallel tokenization tasks
        for (size_t i = 0; i < num_chunks; ++i) {
            const size_t start_pos = i * chunk_size;
            const size_t end_pos = std::min(start_pos + chunk_size, source.size());
            const std::string_view chunk = source.substr(start_pos, end_pos - start_pos);
            
            futures.emplace_back(std::async(std::launch::async, [chunk]() {
                HSMLLexer lexer(chunk);
                return lexer.tokenize();
            }));
        }
        
        // Collect and merge results
        std::vector<HSMLToken> merged_tokens;
        merged_tokens.reserve(source.size() / 8); // Heuristic
        
        for (auto& future : futures) {
            auto chunk_tokens = future.get();
            merged_tokens.insert(merged_tokens.end(),
                               std::make_move_iterator(chunk_tokens.begin()),
                               std::make_move_iterator(chunk_tokens.end()));
        }
        
        return merged_tokens;
    }
}

// Spherical mathematics validation for coordinate tokens
namespace spherical_validation {
    
    [[nodiscard]] constexpr bool is_valid_spherical_coord(double r, double theta, double phi) noexcept {
        return r >= 0.0 &&                    // Radius must be non-negative
               theta >= 0.0 && theta <= M_PI && // Polar angle [0, π]
               phi >= 0.0 && phi < 2.0 * M_PI;  // Azimuthal angle [0, 2π)
    }
    
    [[nodiscard]] constexpr bool is_valid_solid_angle(double omega) noexcept {
        return omega >= 0.0 && omega <= 4.0 * M_PI; // [0, 4π] steradians
    }
    
    // Normalize spherical coordinates to canonical ranges
    [[nodiscard]] constexpr auto normalize_spherical_coord(double r, double theta, double phi) 
        -> std::tuple<double, double, double> {
        
        // Normalize theta to [0, π]
        double norm_theta = std::fmod(std::abs(theta), 2.0 * M_PI);
        if (norm_theta > M_PI) {
            norm_theta = 2.0 * M_PI - norm_theta;
        }
        
        // Normalize phi to [0, 2π)
        double norm_phi = std::fmod(phi, 2.0 * M_PI);
        if (norm_phi < 0.0) {
            norm_phi += 2.0 * M_PI;
        }
        
        return {std::abs(r), norm_theta, norm_phi};
    }
}

// Language-specific token validation
namespace language_validation {
    
    [[nodiscard]] bool is_valid_hsml_structure(const std::vector<HSMLToken>& tokens) {
        size_t brace_depth = 0;
        size_t paren_depth = 0;
        size_t bracket_depth = 0;
        
        for (const auto& token : tokens) {
            switch (token.type) {
                case HSMLTokenType::LEFT_BRACE: ++brace_depth; break;
                case HSMLTokenType::RIGHT_BRACE: 
                    if (brace_depth == 0) return false;
                    --brace_depth; 
                    break;
                case HSMLTokenType::LEFT_PAREN: ++paren_depth; break;
                case HSMLTokenType::RIGHT_PAREN:
                    if (paren_depth == 0) return false;
                    --paren_depth;
                    break;
                case HSMLTokenType::LEFT_BRACKET: ++bracket_depth; break;
                case HSMLTokenType::RIGHT_BRACKET:
                    if (bracket_depth == 0) return false;
                    --bracket_depth;
                    break;
                default: break;
            }
        }
        
        return brace_depth == 0 && paren_depth == 0 && bracket_depth == 0;
    }
    
    [[nodiscard]] bool validate_matter_state_transitions(const std::vector<HSMLToken>& tokens) {
        // Validate physics-based matter state transitions
        std::unordered_set<HSMLTokenType> seen_states;
        
        for (const auto& token : tokens) {
            if (token.type == HSMLTokenType::SOLID ||
                token.type == HSMLTokenType::LIQUID ||
                token.type == HSMLTokenType::GAS ||
                token.type == HSMLTokenType::PLASMA) {
                
                seen_states.insert(token.type);
            }
        }
        
        // Validate that state transitions follow physical laws
        // (This is a simplified validation - real physics would be more complex)
        if (seen_states.contains(HSMLTokenType::SOLID) && 
            seen_states.contains(HSMLTokenType::PLASMA)) {
            // Direct solid-to-plasma transition requires special handling
            return std::ranges::any_of(tokens, [](const HSMLToken& t) {
                return t.value.find("sublimation") != std::string::npos ||
                       t.value.find("ionization") != std::string::npos;
            });
        }
        
        return true;
    }
}

// Export C-style interface for FFI compatibility
extern "C" {
    
    struct HSMLTokenC {
        uint32_t type;
        const char* value;
        size_t value_length;
        size_t line;
        size_t column;
        uint8_t language;
    };
    
    struct HSMLLexerC;
    
    HSMLLexerC* hsml_lexer_create(const char* source, size_t source_length, uint8_t language) {
        try {
            const std::string_view source_view(source, source_length);
            auto* lexer = new HSMLLexer(source_view);
            return reinterpret_cast<HSMLLexerC*>(lexer);
        } catch (...) {
            return nullptr;
        }
    }
    
    void hsml_lexer_destroy(HSMLLexerC* lexer) {
        delete reinterpret_cast<HSMLLexer*>(lexer);
    }
    
    size_t hsml_lexer_tokenize(HSMLLexerC* lexer, HSMLTokenC* output_tokens, size_t max_tokens) {
        try {
            auto* cpp_lexer = reinterpret_cast<HSMLLexer*>(lexer);
            auto tokens = cpp_lexer->tokenize();
            
            const size_t token_count = std::min(tokens.size(), max_tokens);
            
            for (size_t i = 0; i < token_count; ++i) {
                output_tokens[i].type = static_cast<uint32_t>(tokens[i].type);
                output_tokens[i].value = tokens[i].value.c_str();
                output_tokens[i].value_length = tokens[i].value.length();
                output_tokens[i].line = tokens[i].line;
                output_tokens[i].column = tokens[i].column;
                output_tokens[i].language = static_cast<uint8_t>(tokens[i].language);
            }
            
            return token_count;
        } catch (...) {
            return 0;
        }
    }
}

} // namespace hsml::parsers