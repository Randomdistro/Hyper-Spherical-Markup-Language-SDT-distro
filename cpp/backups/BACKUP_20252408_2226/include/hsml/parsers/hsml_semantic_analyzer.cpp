#pragma once

#include "hsml/core/spherical_coords.h"
#include "hsml/core/solid_angle.h"
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <functional>

namespace hsml {
namespace parsers {

using namespace core;

// Forward declarations to avoid cyclic include with hsml_parser.h
struct HSMLElement;
struct AttributeMap;

/**
 * @brief Semantic error types for multi-language validation
 */
enum class SemanticErrorType {
    ERROR,
    WARNING
};

/**
 * @brief Language types supported by the semantic analyzer
 */
enum class AnalysisLanguage {
    HSML,
    CSSS,
    SHAPE,
    STYB
};

/**
 * @brief Semantic error information
 */
struct SemanticError {
    SemanticErrorType type;
    std::string message;
    std::shared_ptr<HSMLElement> node;
    size_t line;
    size_t column;
    AnalysisLanguage language;
    
    SemanticError(SemanticErrorType t, const std::string& msg, 
                  std::shared_ptr<HSMLElement> n, size_t l, size_t c, 
                  AnalysisLanguage lang)
        : type(t), message(msg), node(n), line(l), column(c), language(lang) {}
        
    std::string to_string() const {
        std::string type_str = (type == SemanticErrorType::ERROR) ? "Error" : "Warning";
        std::string lang_str;
        switch (language) {
            case AnalysisLanguage::HSML: lang_str = "HSML"; break;
            case AnalysisLanguage::CSSS: lang_str = "CSSS"; break;
            case AnalysisLanguage::SHAPE: lang_str = "SHAPE"; break;
            case AnalysisLanguage::STYB: lang_str = "STYB"; break;
        }
        return "[" + lang_str + " " + type_str + "] Line " + std::to_string(line) + 
               ":" + std::to_string(column) + " - " + message;
    }
};

/**
 * @brief Variable types supported in semantic analysis
 */
enum class VariableType {
    NUMBER,
    STRING,
    BOOLEAN,
    SPHERICAL_COORD,
    SOLID_ANGLE,
    MATTER_STATE,
    UNKNOWN
};

/**
 * @brief Variable symbol information
 */
struct VariableSymbol {
    std::string name;
    VariableType type;
    std::string scope;
    bool is_constant;
    std::string value; // Stored as string, converted when needed
    
    VariableSymbol(const std::string& n, VariableType t, const std::string& s, 
                   bool const_flag = false, const std::string& v = "")
        : name(n), type(t), scope(s), is_constant(const_flag), value(v) {}
};

/**
 * @brief Function parameter information
 */
struct ParameterSymbol {
    std::string name;
    VariableType type;
    bool required;
    std::string default_value;
    
    ParameterSymbol(const std::string& n, VariableType t, bool req = true, 
                    const std::string& def = "")
        : name(n), type(t), required(req), default_value(def) {}
};

/**
 * @brief Function symbol information
 */
struct FunctionSymbol {
    std::string name;
    std::vector<ParameterSymbol> parameters;
    VariableType return_type;
    std::string scope;
    
    FunctionSymbol(const std::string& n, VariableType ret_type, const std::string& s)
        : name(n), return_type(ret_type), scope(s) {}
};

/**
 * @brief Element symbol for HSML elements
 */
struct ElementSymbol {
    std::string id;
    std::string tag_name;
    std::unordered_map<std::string, std::string> attributes;
    std::vector<std::string> children;
    SphericalCoords position;
    std::string material;
    std::string behavior;
    
    ElementSymbol(const std::string& element_id, const std::string& tag)
        : id(element_id), tag_name(tag), position(0, 0, 0) {}
};

/**
 * @brief Material symbol for CSSS materials
 */
struct MaterialSymbol {
    std::string name;
    std::unordered_map<std::string, std::string> properties;
    std::string matter_state; // solid, liquid, gas, plasma
    std::unordered_map<std::string, double> physics_properties;
    
    MaterialSymbol(const std::string& n, const std::string& state)
        : name(n), matter_state(state) {}
};

/**
 * @brief Physics symbol for ShapeScript physics
 */
struct PhysicsSymbol {
    std::string type; // force, constraint, collision
    std::unordered_map<std::string, std::string> parameters;
    std::string matter_state;
    
    PhysicsSymbol(const std::string& t, const std::string& state)
        : type(t), matter_state(state) {}
};

/**
 * @brief Event symbol for ShapeScript events
 */
struct EventSymbol {
    std::string trigger;
    std::string response;
    std::vector<std::string> conditions;
    
    EventSymbol(const std::string& trig, const std::string& resp)
        : trigger(trig), response(resp) {}
};

/**
 * @brief Constraint symbol for ShapeScript constraints
 */
struct ConstraintSymbol {
    std::string type; // spherical_surface, radial_range, angular_cone
    std::unordered_map<std::string, std::string> parameters;
    
    ConstraintSymbol(const std::string& t) : type(t) {}
};

/**
 * @brief Behavior symbol for ShapeScript behaviors
 */
struct BehaviorSymbol {
    std::string name;
    std::vector<PhysicsSymbol> physics;
    std::vector<EventSymbol> events;
    std::vector<ConstraintSymbol> constraints;
    
    BehaviorSymbol(const std::string& n) : name(n) {}
};

/**
 * @brief Bot symbol for StyleBot agents
 */
struct BotSymbol {
    std::string name;
    std::string agent_type;
    std::unordered_map<std::string, std::string> configuration;
    std::vector<std::string> parallel_tasks;
    
    BotSymbol(const std::string& n, const std::string& type)
        : name(n), agent_type(type) {}
};

/**
 * @brief Complete symbol table for multi-language analysis
 */
struct SymbolTable {
    std::unordered_map<std::string, VariableSymbol> variables;
    std::unordered_map<std::string, FunctionSymbol> functions;
    std::unordered_map<std::string, ElementSymbol> elements;
    std::unordered_map<std::string, MaterialSymbol> materials;
    std::unordered_map<std::string, BehaviorSymbol> behaviors;
    std::unordered_map<std::string, BotSymbol> bots;
    
    void clear() {
        variables.clear();
        functions.clear();
        elements.clear();
        materials.clear();
        behaviors.clear();
        bots.clear();
    }
};

/**
 * @brief Multi-language semantic analyzer for HSML ecosystem
 * 
 * Provides comprehensive semantic analysis across HSML, CSSS, ShapeScript,
 * and StyleBot languages with physics consistency checking and cross-language
 * validation.
 */
class HSMLSemanticAnalyzer {
public:
    HSMLSemanticAnalyzer();
    ~HSMLSemanticAnalyzer() = default;
    
    /**
     * @brief Analyze parsed HSML element tree
     * @param root_element Root element from parser
     * @return True if analysis passed without errors
     */
    bool analyze(std::shared_ptr<HSMLElement> root_element);
    
    /**
     * @brief Analyze raw HSML content (will parse first)
     * @param content Raw HSML content string
     * @return True if analysis passed without errors
     */
    bool analyze_content(const std::string& content);
    
    /**
     * @brief Get all semantic errors found during analysis
     */
    const std::vector<SemanticError>& get_errors() const { return errors_; }
    
    /**
     * @brief Get only error-level issues (excluding warnings)
     */
    std::vector<SemanticError> get_errors_only() const;
    
    /**
     * @brief Get only warning-level issues
     */
    std::vector<SemanticError> get_warnings_only() const;
    
    /**
     * @brief Check if analysis found any errors
     */
    bool has_errors() const;
    
    /**
     * @brief Check if analysis found any warnings
     */
    bool has_warnings() const;
    
    /**
     * @brief Clear all stored errors and warnings
     */
    void clear_errors() { errors_.clear(); }
    
    /**
     * @brief Get the current symbol table
     */
    const SymbolTable& get_symbol_table() const { return symbol_table_; }
    
    /**
     * @brief Configuration methods
     */
    void set_strict_mode(bool strict) { strict_mode_ = strict; }
    bool is_strict_mode() const { return strict_mode_; }
    
    void set_physics_validation_enabled(bool enabled) { physics_validation_ = enabled; }
    bool is_physics_validation_enabled() const { return physics_validation_; }
    
    void set_cross_language_validation_enabled(bool enabled) { cross_language_validation_ = enabled; }
    bool is_cross_language_validation_enabled() const { return cross_language_validation_; }

private:
    SymbolTable symbol_table_;
    std::vector<SemanticError> errors_;
    std::vector<std::string> scope_stack_;
    
    // Configuration
    bool strict_mode_;
    bool physics_validation_;
    bool cross_language_validation_;
    
    // Analysis methods
    void analyze_element(std::shared_ptr<HSMLElement> element, AnalysisLanguage language);
    void analyze_hsml_element(std::shared_ptr<HSMLElement> element);
    void analyze_csss_element(std::shared_ptr<HSMLElement> element);
    void analyze_shape_element(std::shared_ptr<HSMLElement> element);
    void analyze_styb_element(std::shared_ptr<HSMLElement> element);
    
    // Validation methods
    void validate_spherical_coordinates(const AttributeMap& attributes, std::shared_ptr<HSMLElement> element);
    void validate_solid_angle_attributes(const AttributeMap& attributes, std::shared_ptr<HSMLElement> element);
    void validate_matter_state_consistency(const std::string& matter_state, std::shared_ptr<HSMLElement> element);
    void validate_physics_properties(const MaterialSymbol& material, std::shared_ptr<HSMLElement> element);
    void validate_behavior_constraints(const BehaviorSymbol& behavior, std::shared_ptr<HSMLElement> element);
    void validate_cross_language_references();
    
    // Symbol table management
    void enter_scope(const std::string& scope_name);
    void exit_scope();
    std::string current_scope() const;
    
    void declare_variable(const std::string& name, VariableType type, bool is_constant = false, const std::string& value = "");
    void declare_function(const std::string& name, VariableType return_type, const std::vector<ParameterSymbol>& parameters);
    void declare_element(const std::string& id, const std::string& tag_name);
    void declare_material(const std::string& name, const std::string& matter_state);
    void declare_behavior(const std::string& name);
    void declare_bot(const std::string& name, const std::string& agent_type);
    
    // Lookup methods
    VariableSymbol* lookup_variable(const std::string& name);
    FunctionSymbol* lookup_function(const std::string& name);
    ElementSymbol* lookup_element(const std::string& id);
    MaterialSymbol* lookup_material(const std::string& name);
    BehaviorSymbol* lookup_behavior(const std::string& name);
    BotSymbol* lookup_bot(const std::string& name);
    
    // Error reporting
    void report_error(const std::string& message, std::shared_ptr<HSMLElement> element, AnalysisLanguage language);
    void report_warning(const std::string& message, std::shared_ptr<HSMLElement> element, AnalysisLanguage language);
    
    // Utility methods
    AnalysisLanguage detect_language(std::shared_ptr<HSMLElement> element);
    VariableType parse_variable_type(const std::string& type_string);
    std::string variable_type_to_string(VariableType type);
    bool is_valid_spherical_coordinate(double r, double theta, double phi);
    bool is_valid_solid_angle(double omega);
    bool is_valid_matter_state(const std::string& state);
    
    // Physics validation helpers
    bool validate_material_physics_consistency(const MaterialSymbol& material);
    bool validate_behavior_physics_laws(const BehaviorSymbol& behavior);
    bool validate_constraint_feasibility(const ConstraintSymbol& constraint);
};

// Factory methods for creating analyzers with different configurations
namespace factory {

std::unique_ptr<HSMLSemanticAnalyzer> create_strict_analyzer();
std::unique_ptr<HSMLSemanticAnalyzer> create_lenient_analyzer();
std::unique_ptr<HSMLSemanticAnalyzer> create_physics_analyzer();
std::unique_ptr<HSMLSemanticAnalyzer> create_cross_language_analyzer();

} // namespace factory

// Utility functions for semantic analysis
namespace semantic_utils {

std::string format_error_report(const std::vector<SemanticError>& errors);
bool validate_physics_consistency(const SymbolTable& symbol_table);
std::vector<std::string> get_undefined_references(const SymbolTable& symbol_table);
std::vector<std::string> get_unused_symbols(const SymbolTable& symbol_table);

} // namespace semantic_utils

} // namespace parsers
} // namespace hsml