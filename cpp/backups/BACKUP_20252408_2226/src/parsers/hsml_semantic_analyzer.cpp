#include "hsml/parsers/hsml_semantic_analyzer.h"
#include "hsml/parsers/hsml_parser.h"
#include <algorithm>
#include <sstream>
#include <regex>
#include <cmath>

namespace hsml {
namespace parsers {

HSMLSemanticAnalyzer::HSMLSemanticAnalyzer()
    : strict_mode_(false)
    , physics_validation_(true)
    , cross_language_validation_(true) {
    scope_stack_.push_back("global");
}

bool HSMLSemanticAnalyzer::analyze(std::shared_ptr<HSMLElement> root_element) {
    if (!root_element) {
        report_error("Root element is null", nullptr, AnalysisLanguage::HSML);
        return false;
    }
    
    // Clear previous analysis results
    errors_.clear();
    symbol_table_.clear();
    
    // Start analysis from root
    analyze_element(root_element, AnalysisLanguage::HSML);
    
    // Perform cross-language validation if enabled
    if (cross_language_validation_) {
        validate_cross_language_references();
    }
    
    return !has_errors();
}

bool HSMLSemanticAnalyzer::analyze_content(const std::string& content) {
    HSMLParser parser;
    auto root_element = parser.parse_to_elements(content);
    
    if (parser.has_errors()) {
        for (const auto& parse_error : parser.errors()) {
            report_error("Parse error: " + parse_error.message, nullptr, AnalysisLanguage::HSML);
        }
        return false;
    }
    
    return analyze(root_element);
}

std::vector<SemanticError> HSMLSemanticAnalyzer::get_errors_only() const {
    std::vector<SemanticError> error_only;
    std::copy_if(errors_.begin(), errors_.end(), std::back_inserter(error_only),
                 [](const SemanticError& e) { return e.type == SemanticErrorType::ERROR; });
    return error_only;
}

std::vector<SemanticError> HSMLSemanticAnalyzer::get_warnings_only() const {
    std::vector<SemanticError> warnings_only;
    std::copy_if(errors_.begin(), errors_.end(), std::back_inserter(warnings_only),
                 [](const SemanticError& e) { return e.type == SemanticErrorType::WARNING; });
    return warnings_only;
}

bool HSMLSemanticAnalyzer::has_errors() const {
    return std::any_of(errors_.begin(), errors_.end(),
                       [](const SemanticError& e) { return e.type == SemanticErrorType::ERROR; });
}

bool HSMLSemanticAnalyzer::has_warnings() const {
    return std::any_of(errors_.begin(), errors_.end(),
                       [](const SemanticError& e) { return e.type == SemanticErrorType::WARNING; });
}

void HSMLSemanticAnalyzer::analyze_element(std::shared_ptr<HSMLElement> element, AnalysisLanguage language) {
    if (!element) return;
    
    // Detect language if not specified
    if (language == AnalysisLanguage::HSML) {
        language = detect_language(element);
    }
    
    // Dispatch to language-specific analysis
    switch (language) {
        case AnalysisLanguage::HSML:
            analyze_hsml_element(element);
            break;
        case AnalysisLanguage::CSSS:
            analyze_csss_element(element);
            break;
        case AnalysisLanguage::SHAPE:
            analyze_shape_element(element);
            break;
        case AnalysisLanguage::STYB:
            analyze_styb_element(element);
            break;
    }
    
    // Recursively analyze children
    for (auto& child : element->children) {
        analyze_element(child, language);
    }
}

void HSMLSemanticAnalyzer::analyze_hsml_element(std::shared_ptr<HSMLElement> element) {
    const std::string& tag_name = element->tag_name;
    
    // Validate basic HSML element structure
    if (tag_name.empty()) {
        report_error("Empty tag name", element, AnalysisLanguage::HSML);
        return;
    }
    
    // Check for required attributes based on element type
    if (tag_name == "bubble") {
        // Validate spherical coordinates
        validate_spherical_coordinates(element->attributes, element);
        
        // Declare element in symbol table
        std::string element_id = element->attributes.get("id", "anonymous_bubble");
        declare_element(element_id, tag_name);
        
        // Store position in symbol table
        auto* elem_symbol = lookup_element(element_id);
        if (elem_symbol) {
            double r = element->attributes.get_double("r", 1000.0);
            double theta = element->attributes.get_double("theta", 0.0);
            double phi = element->attributes.get_double("phi", 0.0);
            elem_symbol->position = SphericalCoords(r, theta, phi);
        }
        
    } else if (tag_name == "presence") {
        // Validate presence-specific attributes
        if (!element->attributes.has("entity_id")) {
            if (strict_mode_) {
                report_error("Presence element missing required entity_id attribute", element, AnalysisLanguage::HSML);
            } else {
                report_warning("Presence element missing entity_id attribute", element, AnalysisLanguage::HSML);
            }
        }
        
        // Validate influence radius
        double influence = element->attributes.get_double("influence", 1.0);
        if (influence <= 0.0) {
            report_error("Presence influence must be positive", element, AnalysisLanguage::HSML);
        }
        
    } else if (tag_name == "material") {
        // Analyze material properties
        std::string material_name = element->attributes.get("name", "default_material");
        std::string matter_state = element->attributes.get("state", "solid");
        
        if (!is_valid_matter_state(matter_state)) {
            report_error("Invalid matter state: " + matter_state, element, AnalysisLanguage::HSML);
        } else {
            declare_material(material_name, matter_state);
            
            if (physics_validation_) {
                auto* material_symbol = lookup_material(material_name);
                if (material_symbol) {
                    validate_physics_properties(*material_symbol, element);
                }
            }
        }
    }
}

void HSMLSemanticAnalyzer::analyze_csss_element(std::shared_ptr<HSMLElement> element) {
    const std::string& tag_name = element->tag_name;
    
    // Analyze CSSS (physics-based styling) elements
    if (tag_name == "rule") {
        // Validate CSS-like rule structure
        if (!element->attributes.has("selector")) {
            report_error("CSSS rule missing selector", element, AnalysisLanguage::CSSS);
        }
        
    } else if (tag_name == "material") {
        // Enhanced material analysis for CSSS
        std::string material_name = element->attributes.get("name");
        if (material_name.empty()) {
            report_error("CSSS material missing name", element, AnalysisLanguage::CSSS);
        }
        
        // Validate physics properties
        std::string density_str = element->attributes.get("density");
        std::string viscosity_str = element->attributes.get("viscosity");
        
        if (!density_str.empty()) {
            try {
                double density = std::stod(density_str);
                if (density <= 0.0) {
                    report_error("Material density must be positive", element, AnalysisLanguage::CSSS);
                }
            } catch (const std::exception&) {
                report_error("Invalid density value: " + density_str, element, AnalysisLanguage::CSSS);
            }
        }
        
    } else if (tag_name == "animation") {
        // Validate animation keyframes and physics
        if (!element->attributes.has("duration")) {
            report_warning("Animation missing duration", element, AnalysisLanguage::CSSS);
        }
        
        // Check for physics-based animation properties
        if (element->attributes.has("force") || element->attributes.has("momentum")) {
            if (!physics_validation_) {
                report_warning("Physics-based animation without physics validation", element, AnalysisLanguage::CSSS);
            }
        }
    }
}

void HSMLSemanticAnalyzer::analyze_shape_element(std::shared_ptr<HSMLElement> element) {
    const std::string& tag_name = element->tag_name;
    
    // Analyze ShapeScript (behavior definition) elements
    if (tag_name == "behavior") {
        std::string behavior_name = element->attributes.get("name");
        if (behavior_name.empty()) {
            report_error("ShapeScript behavior missing name", element, AnalysisLanguage::SHAPE);
        } else {
            declare_behavior(behavior_name);
        }
        
    } else if (tag_name == "physics") {
        // Validate physics constraints and forces
        std::string physics_type = element->attributes.get("type", "force");
        std::string matter_state = element->attributes.get("matter_state", "solid");
        
        if (!is_valid_matter_state(matter_state)) {
            report_error("Invalid matter state in physics: " + matter_state, element, AnalysisLanguage::SHAPE);
        }
        
        // Validate force/constraint parameters
        if (physics_type == "force") {
            if (!element->attributes.has("magnitude") && !element->attributes.has("vector")) {
                report_error("Physics force missing magnitude or vector", element, AnalysisLanguage::SHAPE);
            }
        } else if (physics_type == "constraint") {
            if (!element->attributes.has("type")) {
                report_error("Physics constraint missing type", element, AnalysisLanguage::SHAPE);
            }
        }
        
    } else if (tag_name == "event") {
        // Validate event triggers and responses
        if (!element->attributes.has("trigger")) {
            report_error("ShapeScript event missing trigger", element, AnalysisLanguage::SHAPE);
        }
        if (!element->attributes.has("response")) {
            report_error("ShapeScript event missing response", element, AnalysisLanguage::SHAPE);
        }
    }
}

void HSMLSemanticAnalyzer::analyze_styb_element(std::shared_ptr<HSMLElement> element) {
    const std::string& tag_name = element->tag_name;
    
    // Analyze StyleBot (autonomous agent) elements
    if (tag_name == "bot") {
        std::string bot_name = element->attributes.get("name");
        std::string agent_type = element->attributes.get("type", "renderer");
        
        if (bot_name.empty()) {
            report_error("StyleBot missing name", element, AnalysisLanguage::STYB);
        } else {
            declare_bot(bot_name, agent_type);
        }
        
        // Validate agent type
        std::vector<std::string> valid_types = {"renderer", "physics", "optimizer", "coordinator"};
        if (std::find(valid_types.begin(), valid_types.end(), agent_type) == valid_types.end()) {
            report_warning("Unknown StyleBot agent type: " + agent_type, element, AnalysisLanguage::STYB);
        }
        
    } else if (tag_name == "parallel") {
        // Validate parallel task definitions
        if (!element->attributes.has("tasks")) {
            report_error("StyleBot parallel block missing tasks", element, AnalysisLanguage::STYB);
        }
        
        // Check for potential race conditions
        std::string sync_mode = element->attributes.get("sync", "async");
        if (sync_mode == "sync" && element->children.size() > 1) {
            report_warning("Synchronous parallel tasks may cause performance issues", element, AnalysisLanguage::STYB);
        }
        
    } else if (tag_name == "render") {
        // Validate rendering agent configuration
        std::string render_mode = element->attributes.get("mode", "realtime");
        if (render_mode != "realtime" && render_mode != "batch" && render_mode != "adaptive") {
            report_error("Invalid StyleBot render mode: " + render_mode, element, AnalysisLanguage::STYB);
        }
    }
}

void HSMLSemanticAnalyzer::validate_spherical_coordinates(const AttributeMap& attributes, std::shared_ptr<HSMLElement> element) {
    double r = attributes.get_double("r", 1000.0);
    double theta = attributes.get_double("theta", 0.0);
    double phi = attributes.get_double("phi", 0.0);
    
    if (!is_valid_spherical_coordinate(r, theta, phi)) {
        report_error("Invalid spherical coordinates: r=" + std::to_string(r) + 
                    ", theta=" + std::to_string(theta) + 
                    ", phi=" + std::to_string(phi), element, AnalysisLanguage::HSML);
    }
    
    // Additional validation for viewer distance considerations
    if (r < 1.0) {
        report_warning("Radial distance very small (r=" + std::to_string(r) + ")", element, AnalysisLanguage::HSML);
    }
    if (r > 100000.0) {
        report_warning("Radial distance very large (r=" + std::to_string(r) + ")", element, AnalysisLanguage::HSML);
    }
}

bool HSMLSemanticAnalyzer::is_valid_spherical_coordinate(double r, double theta, double phi) {
    return r >= 0.0 && 
           theta >= 0.0 && theta <= M_PI && 
           phi >= -M_PI && phi <= M_PI &&
           std::isfinite(r) && std::isfinite(theta) && std::isfinite(phi);
}

bool HSMLSemanticAnalyzer::is_valid_matter_state(const std::string& state) {
    return state == "solid" || state == "liquid" || state == "gas" || state == "plasma";
}

AnalysisLanguage HSMLSemanticAnalyzer::detect_language(std::shared_ptr<HSMLElement> element) {
    const std::string& tag = element->tag_name;
    
    // CSSS elements
    if (tag == "rule" || tag == "selector" || tag == "declaration" || 
        tag == "animation" || tag == "keyframe" || tag.find("csss") != std::string::npos) {
        return AnalysisLanguage::CSSS;
    }
    
    // ShapeScript elements
    if (tag == "behavior" || tag == "physics" || tag == "force" || 
        tag == "constraint" || tag == "event" || tag.find("shape") != std::string::npos) {
        return AnalysisLanguage::SHAPE;
    }
    
    // StyleBot elements
    if (tag == "bot" || tag == "agent" || tag == "parallel" || 
        tag == "render" || tag.find("styb") != std::string::npos) {
        return AnalysisLanguage::STYB;
    }
    
    // Default to HSML
    return AnalysisLanguage::HSML;
}

void HSMLSemanticAnalyzer::declare_element(const std::string& id, const std::string& tag_name) {
    symbol_table_.elements[id] = ElementSymbol(id, tag_name);
}

void HSMLSemanticAnalyzer::declare_material(const std::string& name, const std::string& matter_state) {
    symbol_table_.materials[name] = MaterialSymbol(name, matter_state);
}

void HSMLSemanticAnalyzer::declare_behavior(const std::string& name) {
    symbol_table_.behaviors[name] = BehaviorSymbol(name);
}

void HSMLSemanticAnalyzer::declare_bot(const std::string& name, const std::string& agent_type) {
    symbol_table_.bots[name] = BotSymbol(name, agent_type);
}

ElementSymbol* HSMLSemanticAnalyzer::lookup_element(const std::string& id) {
    auto it = symbol_table_.elements.find(id);
    return (it != symbol_table_.elements.end()) ? &it->second : nullptr;
}

MaterialSymbol* HSMLSemanticAnalyzer::lookup_material(const std::string& name) {
    auto it = symbol_table_.materials.find(name);
    return (it != symbol_table_.materials.end()) ? &it->second : nullptr;
}

void HSMLSemanticAnalyzer::validate_physics_properties(const MaterialSymbol& material, std::shared_ptr<HSMLElement> element) {
    // Validate material physics consistency based on matter state
    if (material.matter_state == "solid") {
        // Solids should have high density, low compressibility
        auto density_it = material.physics_properties.find("density");
        if (density_it != material.physics_properties.end() && density_it->second < 100.0) {
            report_warning("Solid material with unusually low density", element, AnalysisLanguage::CSSS);
        }
    } else if (material.matter_state == "gas") {
        // Gases should have low density, high compressibility
        auto density_it = material.physics_properties.find("density");
        if (density_it != material.physics_properties.end() && density_it->second > 10.0) {
            report_warning("Gas material with unusually high density", element, AnalysisLanguage::CSSS);
        }
    }
}

void HSMLSemanticAnalyzer::validate_cross_language_references() {
    // Check for undefined material references in HSML elements
    for (const auto& [id, element] : symbol_table_.elements) {
        if (!element.material.empty()) {
            if (symbol_table_.materials.find(element.material) == symbol_table_.materials.end()) {
                report_error("Undefined material reference: " + element.material, nullptr, AnalysisLanguage::HSML);
            }
        }
        
        if (!element.behavior.empty()) {
            if (symbol_table_.behaviors.find(element.behavior) == symbol_table_.behaviors.end()) {
                report_error("Undefined behavior reference: " + element.behavior, nullptr, AnalysisLanguage::HSML);
            }
        }
    }
}

void HSMLSemanticAnalyzer::report_error(const std::string& message, std::shared_ptr<HSMLElement> element, AnalysisLanguage language) {
    errors_.emplace_back(SemanticErrorType::ERROR, message, element, 0, 0, language);
}

void HSMLSemanticAnalyzer::report_warning(const std::string& message, std::shared_ptr<HSMLElement> element, AnalysisLanguage language) {
    errors_.emplace_back(SemanticErrorType::WARNING, message, element, 0, 0, language);
}

// Factory implementations
namespace factory {

std::unique_ptr<HSMLSemanticAnalyzer> create_strict_analyzer() {
    auto analyzer = std::make_unique<HSMLSemanticAnalyzer>();
    analyzer->set_strict_mode(true);
    analyzer->set_physics_validation_enabled(true);
    analyzer->set_cross_language_validation_enabled(true);
    return analyzer;
}

std::unique_ptr<HSMLSemanticAnalyzer> create_lenient_analyzer() {
    auto analyzer = std::make_unique<HSMLSemanticAnalyzer>();
    analyzer->set_strict_mode(false);
    analyzer->set_physics_validation_enabled(false);
    analyzer->set_cross_language_validation_enabled(false);
    return analyzer;
}

std::unique_ptr<HSMLSemanticAnalyzer> create_physics_analyzer() {
    auto analyzer = std::make_unique<HSMLSemanticAnalyzer>();
    analyzer->set_strict_mode(false);
    analyzer->set_physics_validation_enabled(true);
    analyzer->set_cross_language_validation_enabled(true);
    return analyzer;
}

} // namespace factory

} // namespace parsers
} // namespace hsml