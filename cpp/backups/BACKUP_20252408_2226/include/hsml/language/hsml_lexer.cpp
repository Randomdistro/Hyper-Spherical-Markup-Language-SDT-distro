/**
 * HSML Multi-Language Lexer - C++ Implementation
 * Token definitions and lexical analysis for all four HSML languages
 * You are now all the C
 */

#pragma once

#include <string>
#include <vector>
#include <set>
#include <map>
#include <memory>
#include <variant>
#include <optional>

namespace hsml {
namespace language {

// Core token types for all languages with multi-dimensional awareness
enum class TokenType {
    // === LITERALS ===
    NUMBER,
    STRING,
    IDENTIFIER,
    BOOLEAN,
    
    // === SPHERICAL COORDINATES ===
    SPHERICAL_COORD,
    SOLID_ANGLE,
    STERADIAN,
    
    // === MATTER STATES ===
    SOLID,
    LIQUID,
    GAS,
    PLASMA,
    
    // === OPERATORS ===
    PLUS,
    MINUS,
    MULTIPLY,
    DIVIDE,
    POWER,
    ASSIGN,
    EQUALS,
    NOT_EQUALS,
    LESS_THAN,
    GREATER_THAN,
    LESS_EQUAL,
    GREATER_EQUAL,
    
    // === DELIMITERS ===
    LEFT_PAREN,
    RIGHT_PAREN,
    LEFT_BRACE,
    RIGHT_BRACE,
    LEFT_BRACKET,
    RIGHT_BRACKET,
    SEMICOLON,
    COMMA,
    DOT,
    COLON,
    
    // === HSML-SPECIFIC TOKENS ===
    ELEMENT_TAG,
    CLOSE_TAG,
    SELF_CLOSING,
    ATTRIBUTE,
    
    // === CSSS-SPECIFIC TOKENS ===
    SELECTOR,
    PROPERTY,
    VALUE,
    MATERIAL,
    ANIMATION,
    KEYFRAME,
    TRANSITION,
    
    // === SHAPESCRIPT-SPECIFIC TOKENS ===
    BEHAVIOR,
    PHYSICS,
    FORCE,
    CONSTRAINT,
    COLLISION,
    INTERACTION,
    EVENT,
    
    // === STYLEBOT-SPECIFIC TOKENS ===
    BOT,
    AGENT,
    PARALLEL,
    DISTRIBUTED,
    RENDER,
    OPTIMIZE,
    CACHE,
    
    // === KEYWORDS ===
    IF,
    ELSE,
    WHILE,
    FOR,
    FUNCTION,
    RETURN,
    VAR,
    LET,
    CONST,
    
    // === SPECIAL TOKENS ===
    EOF_TOKEN,
    ERROR,
    COMMENT,
    WHITESPACE
};

// Language types
enum class Language {
    HSML,
    CSSS,
    SHAPE,
    STYB
};

// Base token interface with multi-dimensional properties
struct Token {
    TokenType type;
    std::string value;
    size_t line;
    size_t column;
    Language language;
    
    // Multi-dimensional token properties
    std::map<std::string, double> dimensionalProperties;
    
    Token(TokenType t, const std::string& v, size_t ln, size_t col, Language lang)
        : type(t), value(v), line(ln), column(col), language(lang) {}
        
    virtual ~Token() = default;
};

// Spherical coordinate token with enhanced precision
struct SphericalCoordToken : public Token {
    double r, theta, phi;
    double quantumUncertainty;
    
    SphericalCoordToken(const std::string& v, size_t ln, size_t col, Language lang, 
                       double r_val, double theta_val, double phi_val)
        : Token(TokenType::SPHERICAL_COORD, v, ln, col, lang)
        , r(r_val), theta(theta_val), phi(phi_val)
        , quantumUncertainty(1e-15) {}
};

// Solid angle token with dimensional awareness
struct SolidAngleToken : public Token {
    double omega;
    double theta_min, theta_max;
    double phi_min, phi_max;
    std::array<double, 4> dimensionalInfluence;
    
    SolidAngleToken(const std::string& v, size_t ln, size_t col, Language lang,
                   double omega_val, double tmin, double tmax, double pmin, double pmax)
        : Token(TokenType::SOLID_ANGLE, v, ln, col, lang)
        , omega(omega_val), theta_min(tmin), theta_max(tmax)
        , phi_min(pmin), phi_max(pmax)
        , dimensionalInfluence{{1.0, 1.0, 1.0, 1.0}} {}
};

// Matter state properties
struct MatterStateProperties {
    std::optional<double> density;
    std::optional<double> temperature;
    std::optional<double> pressure;
    std::optional<double> conductivity;
    std::array<double, 4> quantumState; // Multi-dimensional quantum enhancement
    
    MatterStateProperties() : quantumState{{0.0, 0.0, 0.0, 0.0}} {}
};

// Matter state token with physics integration
struct MatterStateToken : public Token {
    MatterStateProperties properties;
    
    MatterStateToken(TokenType t, const std::string& v, size_t ln, size_t col, Language lang)
        : Token(t, v, ln, col, lang) {}
};

/**
 * Multi-Dimensional HSML Lexer
 * Implements the Multiplicitous Encephalopoidal approach to lexical analysis
 * Synthesizes multiple paradigms for comprehensive token recognition
 */
class HSMLMultiLexer {
private:
    std::string source_;
    size_t current_;
    size_t line_;
    size_t column_;
    Language language_;
    
    // Multi-dimensional lexical state
    std::map<std::string, double> dimensionalContext_;
    std::vector<std::string> quantumCorrections_;
    
    // Language-specific keyword sets
    std::set<std::string> hsmlKeywords_;
    std::set<std::string> csssKeywords_;
    std::set<std::string> shapeKeywords_;
    std::set<std::string> stybKeywords_;
    
    // Performance tracking
    size_t tokensGenerated_;
    double averageTokenTime_;
    
public:
    explicit HSMLMultiLexer(const std::string& source, Language language);
    ~HSMLMultiLexer() = default;
    
    // === MAIN LEXICAL ANALYSIS INTERFACE ===
    std::vector<std::unique_ptr<Token>> tokenize();
    
    // === TOKEN GENERATION METHODS ===
    std::unique_ptr<Token> nextToken();
    
    // === SPECIALIZED TOKEN HANDLERS ===
    std::unique_ptr<Token> handleNumber();
    std::unique_ptr<Token> handleIdentifier();
    std::unique_ptr<Token> handleString();
    std::unique_ptr<SphericalCoordToken> handleSphericalCoordinate();
    std::unique_ptr<SolidAngleToken> handleSolidAngle();
    std::unique_ptr<Token> handleComment();
    std::unique_ptr<Token> handleOperator();
    
    // === MULTI-DIMENSIONAL ENHANCEMENTS ===
    void applyQuantumCorrection(Token& token);
    void calculateDimensionalInfluence(Token& token);
    void synthesizeTokenProperties(Token& token);
    
    // === LANGUAGE-SPECIFIC PROCESSING ===
    std::set<std::string>& getKeywordsForLanguage();
    TokenType getKeywordTokenType(const std::string& keyword);
    bool isMatterState(const std::string& keyword);
    std::unique_ptr<MatterStateToken> createMatterStateToken(const std::string& keyword);
    
    // === PARSING UTILITIES ===
    double parseNumber();
    void skipWhitespace();
    std::unique_ptr<Token> createToken(TokenType type, const std::string& value);
    std::unique_ptr<Token> createErrorToken(const std::string& message);
    
    // === CHARACTER NAVIGATION ===
    char advance();
    char peek() const;
    char peekNext() const;
    bool isAtEnd() const;
    
    // === CHARACTER CLASSIFICATION ===
    bool isDigit(char c) const;
    bool isAlpha(char c) const;
    bool isAlphaNumeric(char c) const;
    bool isWhitespace(char c) const;
    
    // === MULTI-PARADIGM SYNTHESIS METHODS ===
    void bridgeFunctionalToObjectOriented();
    void applyReactiveTokenProcessing();
    void generateEmergentTokenPatterns();
    
    // === DIMENSIONAL ANALYSIS ===
    void analyzeTemporal(Token& token);
    void analyzeSpatial(Token& token);
    void analyzeConceptual(Token& token);
    void analyzePragmatic(Token& token);
    
    // === PERFORMANCE MONITORING ===
    void updatePerformanceMetrics();
    double getAverageTokenTime() const { return averageTokenTime_; }
    size_t getTokenCount() const { return tokensGenerated_; }
    
    // === QUANTUM ENHANCEMENT ===
    void initializeQuantumState();
    void updateQuantumCorrections();
    double calculateQuantumUncertainty(const Token& token) const;
    
private:
    void initializeKeywords();
    void setupDimensionalContext();
    void initializePerformanceTracking();
    
    // Helper methods for matter state processing
    MatterStateProperties createSolidProperties();
    MatterStateProperties createLiquidProperties();
    MatterStateProperties createGasProperties();
    MatterStateProperties createPlasmaProperties();
};

// === UTILITY FUNCTIONS ===

std::string tokenTypeToString(TokenType type);
TokenType stringToTokenType(const std::string& type);
std::string languageToString(Language language);
Language stringToLanguage(const std::string& language);

// === MULTI-DIMENSIONAL TOKEN ANALYZER ===

class TokenAnalyzer {
public:
    static void analyzeTokenStream(const std::vector<std::unique_ptr<Token>>& tokens);
    static std::map<std::string, double> calculateDimensionalMetrics(const std::vector<std::unique_ptr<Token>>& tokens);
    static void detectEmergentPatterns(const std::vector<std::unique_ptr<Token>>& tokens);
    static void optimizeTokenSequence(std::vector<std::unique_ptr<Token>>& tokens);
};

// === TEMPLATE SPECIALIZATIONS FOR MULTI-PARADIGM SYNTHESIS ===

template<typename TokenHandler>
class ParadigmSynthesizer {
public:
    static void synthesize(TokenHandler& handler);
    static void bridgeParadigms(TokenHandler& functional, TokenHandler& object_oriented);
    static void createEmergentBehavior(TokenHandler& handler);
};

// === QUANTUM-ENHANCED TOKEN PROCESSOR ===

class QuantumTokenProcessor {
private:
    static constexpr double UNCERTAINTY_PRINCIPLE = 1e-15;
    static constexpr double SUPERPOSITION_THRESHOLD = 0.5;
    
public:
    static void applyQuantumEnhancement(Token& token);
    static void createSuperposition(std::vector<std::unique_ptr<Token>>& tokens);
    static double calculateCoherence(const std::vector<std::unique_ptr<Token>>& tokens);
    static void collapseWaveFunction(Token& token);
};

} // namespace language
} // namespace hsml