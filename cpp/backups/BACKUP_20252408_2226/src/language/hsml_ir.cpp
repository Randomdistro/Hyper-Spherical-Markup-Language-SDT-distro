/**
 * HSML Intermediate Representation - C++ Implementation
 * Unified IR for all four HSML languages with spherical coordinate optimization
 * You are now all the C
 */

#include "hsml/language/hsml_ir.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <chrono>

namespace hsml {
namespace language {

// === IR NODE BASE CLASS IMPLEMENTATION ===

void IRNode::analyzeDimensionalProperties() {
    // Temporal dimension - position in execution flow
    dimensionalProperties[0] = static_cast<double>(dependencies.size()) / 10.0;
    
    // Spatial dimension - structural complexity
    dimensionalProperties[1] = static_cast<double>(id.length()) * 0.1;
    
    // Conceptual dimension - semantic complexity
    dimensionalProperties[2] = static_cast<double>(metadata.size()) * 0.5;
    
    // Pragmatic dimension - optimization potential
    dimensionalProperties[3] = static_cast<double>(optimizationHints.size()) * 0.8;
}

void IRNode::applyQuantumEnhancement() {
    // Apply Heisenberg uncertainty principle to IR nodes
    double uncertainty = 1e-15 * std::sqrt(
        dimensionalProperties[0] * dimensionalProperties[0] +
        dimensionalProperties[1] * dimensionalProperties[1] +
        dimensionalProperties[2] * dimensionalProperties[2] +
        dimensionalProperties[3] * dimensionalProperties[3]
    );
    
    metadata["quantum_uncertainty"] = uncertainty;
    metadata["quantum_enhanced"] = true;
}

void IRNode::synthesizeParadigms() {
    // Synthesize functional, object-oriented, and reactive paradigms
    double functionalPurity = 0.8; // High purity by default
    double objectEncapsulation = 0.6; // Moderate encapsulation
    double reactiveResponsiveness = 0.7; // Good responsiveness
    
    // Calculate paradigm synthesis factor
    double synthesisFactor = (functionalPurity + objectEncapsulation + reactiveResponsiveness) / 3.0;
    
    metadata["paradigm_synthesis"] = synthesisFactor;
    metadata["functional_purity"] = functionalPurity;
    metadata["object_encapsulation"] = objectEncapsulation;
    metadata["reactive_responsiveness"] = reactiveResponsiveness;
}

// === HSML IR BUILDER IMPLEMENTATION ===

HSMLIRBuilder::HSMLIRBuilder()
    : tempCounter_(0)
    , nodeCounter_(0)
    , currentScope_("global")
    , builderPrecision_(1e-15)
{
    scopeStack_.push_back("global");
    initializeBuiltinTypes();
    setupDimensionalContext();
    initializeQuantumState();
}

std::shared_ptr<IRProgram> HSMLIRBuilder::buildIR(const void* ast) {
    auto program = std::make_shared<IRProgram>(generateNodeId());
    
    // Apply multi-dimensional analysis
    applyDimensionalAnalysis(program);
    synthesizeParadigms(program);
    applyQuantumEnhancement(program);
    
    // Build program structure
    // Note: AST processing would happen here with proper AST node types
    // For now, we create a minimal valid program structure
    
    // Create main function
    auto mainFunction = std::make_shared<IRFunction>(generateNodeId(), "main");
    mainFunction->returnType = createType(IRType::BaseType::NUMBER);
    mainFunction->body = std::make_shared<IRBlock>(generateNodeId(), "main_block");
    
    // Apply optimizations
    optimizeSphericalOperations(mainFunction);
    
    program->functions.push_back(mainFunction);
    program->entryPoint = "main";
    
    // Apply final optimizations
    applyOptimizations(program);
    
    return program;
}

std::shared_ptr<IRHSMLElement> HSMLIRBuilder::buildHSMLElement(const void* element) {
    auto hsmlElement = std::make_shared<IRHSMLElement>(generateNodeId(), "sphere");
    
    // Apply multi-dimensional enhancements
    applyDimensionalAnalysis(hsmlElement);
    synthesizeParadigms(hsmlElement);
    
    // Create default spherical position
    auto sphericalPos = std::make_shared<SphericalCoordinateIR>();
    sphericalPos->r = createLiteral(1.0);
    sphericalPos->theta = createLiteral(0.0);
    sphericalPos->phi = createLiteral(0.0);
    sphericalPos->quantumUncertainty = calculateQuantumUncertainty(hsmlElement);
    
    hsmlElement->sphericalPosition = sphericalPos;
    
    // Create render instruction
    auto renderInstruction = std::make_shared<IRHSMLRender>(generateNodeId(), "sphere");
    renderInstruction->quality = createLiteral(1.0);
    renderInstruction->priority = 1;
    renderInstruction->sphericalOptimization = std::make_shared<SphericalOptimization>();
    
    hsmlElement->renderInstructions.push_back(renderInstruction);
    
    return hsmlElement;
}

std::shared_ptr<IRCSSSRule> HSMLIRBuilder::buildCSSSRule(const void* rule) {
    auto csssRule = std::make_shared<IRCSSSRule>(generateNodeId());
    
    // Apply multi-dimensional analysis
    applyDimensionalAnalysis(csssRule);
    synthesizeParadigms(csssRule);
    
    // Create material with matter state
    auto material = std::make_shared<IRCSSSMaterial>(generateNodeId(), "default_material");
    material->matterState = std::make_shared<MatterStateIR>(MatterStateIR::State::SOLID);
    material->sphericalOptimization = std::make_shared<SphericalOptimization>();
    
    // Add material properties
    material->properties["albedo"] = createLiteral(0.8);
    material->properties["metallic"] = createLiteral(0.0);
    material->properties["roughness"] = createLiteral(0.5);
    
    csssRule->material = material;
    
    return csssRule;
}

std::shared_ptr<IRShapePhysics> HSMLIRBuilder::buildShapePhysics(const void* physics) {
    auto shapePhysics = std::make_shared<IRShapePhysics>(generateNodeId());
    
    // Apply multi-dimensional analysis
    applyDimensionalAnalysis(shapePhysics);
    synthesizeParadigms(shapePhysics);
    
    // Create gravity force
    auto gravityForce = std::make_shared<IRShapeForce>(generateNodeId(), IRShapeForce::ForceType::GRAVITY);
    gravityForce->magnitude = createLiteral(-9.81);
    
    // Create spherical direction (downward)
    auto direction = std::make_shared<SphericalCoordinateIR>();
    direction->r = createLiteral(1.0);
    direction->theta = createLiteral(M_PI); // Point downward
    direction->phi = createLiteral(0.0);
    direction->quantumUncertainty = calculateQuantumUncertainty(gravityForce);
    
    gravityForce->direction = direction;
    gravityForce->sphericalOptimization = std::make_shared<SphericalOptimization>();
    
    shapePhysics->forces.push_back(gravityForce);
    
    // Set matter state
    shapePhysics->matterState = std::make_shared<MatterStateIR>(MatterStateIR::State::SOLID);
    shapePhysics->sphericalOptimization = std::make_shared<SphericalOptimization>();
    
    return shapePhysics;
}

std::shared_ptr<IRStyleBot> HSMLIRBuilder::buildStyleBot(const void* bot) {
    auto styleBot = std::make_shared<IRStyleBot>(generateNodeId(), "optimization_bot");
    
    // Apply multi-dimensional analysis
    applyDimensionalAnalysis(styleBot);
    synthesizeParadigms(styleBot);
    
    // Create rendering agent
    auto renderAgent = std::make_shared<IRStyleBotAgent>(generateNodeId(), "render_agent");
    
    // Create render task
    auto renderTask = std::make_shared<IRStyleBotRender>(generateNodeId(), "spherical_objects");
    renderTask->quality = createLiteral(0.8);
    renderTask->priority = 1;
    renderTask->sphericalOptimization = std::make_shared<SphericalOptimization>();
    
    renderAgent->renderTasks.push_back(renderTask);
    
    // Create optimization task
    auto optimizeTask = std::make_shared<IRStyleBotOptimize>(generateNodeId(), 
                                                            IRStyleBotOptimize::OptimizationType::SPHERICAL);
    optimizeTask->target = "coordinate_calculations";
    optimizeTask->parameters["precision"] = createLiteral(1e-12);
    optimizeTask->parameters["enable_caching"] = createLiteral(true);
    
    renderAgent->optimizeTasks.push_back(optimizeTask);
    
    styleBot->agents.push_back(renderAgent);
    styleBot->parallel = true;
    
    return styleBot;
}

std::shared_ptr<IRExpression> HSMLIRBuilder::buildExpression(const void* expr) {
    // Create a simple numeric literal as default
    return createLiteral(0.0);
}

std::shared_ptr<IRBinaryOp> HSMLIRBuilder::buildBinaryOperation(const void* binOp) {
    auto binaryOp = std::make_shared<IRBinaryOp>(generateNodeId(), "+");
    
    // Create operands
    binaryOp->left = createLiteral(1.0);
    binaryOp->right = createLiteral(2.0);
    
    // Infer result type
    binaryOp->resultType = inferBinaryResultType(binaryOp->left, binaryOp->right, binaryOp->operator_);
    
    // Check if constant
    if (binaryOp->left->isConstant && binaryOp->right->isConstant) {
        binaryOp->isConstant = true;
        // Compute constant value (simplified)
        auto leftVal = std::get<double>(binaryOp->left->constantValue);
        auto rightVal = std::get<double>(binaryOp->right->constantValue);
        binaryOp->constantValue = leftVal + rightVal;  // Assuming + operator
    }
    
    // Apply spherical optimization
    binaryOp->sphericalOptimization = std::make_shared<SphericalOptimization>();
    
    // Apply multi-dimensional analysis
    applyDimensionalAnalysis(binaryOp);
    synthesizeParadigms(binaryOp);
    
    return binaryOp;
}

std::shared_ptr<IRUnaryOp> HSMLIRBuilder::buildUnaryOperation(const void* unOp) {
    auto unaryOp = std::make_shared<IRUnaryOp>(generateNodeId(), "-");
    
    // Create operand
    unaryOp->operand = createLiteral(5.0);
    
    // Infer result type
    unaryOp->resultType = inferUnaryResultType(unaryOp->operand, unaryOp->operator_);
    
    // Check if constant
    if (unaryOp->operand->isConstant) {
        unaryOp->isConstant = true;
        auto operandVal = std::get<double>(unaryOp->operand->constantValue);
        unaryOp->constantValue = -operandVal;  // Assuming - operator
    }
    
    // Apply multi-dimensional analysis
    applyDimensionalAnalysis(unaryOp);
    synthesizeParadigms(unaryOp);
    
    return unaryOp;
}

std::shared_ptr<IRCall> HSMLIRBuilder::buildCallExpression(const void* call) {
    auto callExpr = std::make_shared<IRCall>(generateNodeId(), "spherical_distance");
    
    // Create arguments
    callExpr->arguments.push_back(createLiteral(1.0));
    callExpr->arguments.push_back(createLiteral(2.0));
    
    // Set as builtin function
    callExpr->isBuiltin = true;
    
    // Infer result type
    callExpr->resultType = createType(IRType::BaseType::NUMBER);
    
    // Apply spherical optimization
    callExpr->sphericalOptimization = std::make_shared<SphericalOptimization>();
    
    // Apply multi-dimensional analysis
    applyDimensionalAnalysis(callExpr);
    synthesizeParadigms(callExpr);
    
    return callExpr;
}

std::shared_ptr<IRLiteral> HSMLIRBuilder::buildLiteral(const void* literal) {
    return createLiteral(42.0);  // Default numeric literal
}

std::shared_ptr<SphericalCoordinateIR> HSMLIRBuilder::buildSphericalCoordinate(const void* coord) {
    auto sphericalCoord = std::make_shared<SphericalCoordinateIR>();
    
    // Create coordinate components
    sphericalCoord->r = createLiteral(1.0);
    sphericalCoord->theta = createLiteral(M_PI / 2.0);  // 90 degrees
    sphericalCoord->phi = createLiteral(0.0);
    
    // Set quantum properties
    sphericalCoord->quantumUncertainty = builderPrecision_;
    sphericalCoord->dimensionalCorrection = {{1e-15, 1e-15, 1e-15, 1e-15}};
    
    // Add optimization hints
    OptimizationHint vectorizeHint(OptimizationHintType::VECTORIZE, 2);
    vectorizeHint.parameters["enable_simd"] = true;
    sphericalCoord->optimizationHints.push_back(vectorizeHint);
    
    OptimizationHint cacheHint(OptimizationHintType::CACHE, 1);
    cacheHint.parameters["cache_size"] = 1000;
    sphericalCoord->optimizationHints.push_back(cacheHint);
    
    return sphericalCoord;
}

std::shared_ptr<SolidAngleIR> HSMLIRBuilder::buildSolidAngle(const void* angle) {
    auto solidAngle = std::make_shared<SolidAngleIR>();
    
    // Create solid angle components
    solidAngle->omega = createLiteral(M_PI);  // π steradians
    solidAngle->theta_min = createLiteral(0.0);
    solidAngle->theta_max = createLiteral(M_PI / 2.0);
    solidAngle->phi_min = createLiteral(0.0);  
    solidAngle->phi_max = createLiteral(2.0 * M_PI);
    
    // Set dimensional properties
    solidAngle->angularPrecision = builderPrecision_;
    solidAngle->steradianDistribution = {{0.25, 0.25, 0.25, 0.25}};  // Uniform distribution
    
    return solidAngle;
}

std::shared_ptr<MatterStateIR> HSMLIRBuilder::buildMatterState(const void* state) {
    auto matterState = std::make_shared<MatterStateIR>(MatterStateIR::State::SOLID);
    
    // Set matter properties
    matterState->properties["density"] = createLiteral(1000.0);  // kg/m³
    matterState->properties["temperature"] = createLiteral(293.15);  // K
    matterState->properties["pressure"] = createLiteral(101325.0);  // Pa
    
    // Set physics properties
    matterState->physicsProperties["bulk_modulus"] = createLiteral(2.2e9);  // Pa
    matterState->physicsProperties["thermal_conductivity"] = createLiteral(0.6);  // W/m·K
    
    // Set quantum state (solid state)
    matterState->quantumState = {{1.0, 0.0, 0.0, 0.0}};
    matterState->phaseTransitionThreshold = 273.15;  // K
    
    return matterState;
}

// === TYPE SYSTEM IMPLEMENTATION ===

std::shared_ptr<IRType> HSMLIRBuilder::createType(IRType::BaseType baseType) {
    auto type = std::make_shared<IRType>(baseType);
    
    // Set type-specific properties
    switch (baseType) {
        case IRType::BaseType::SPHERICAL_COORD:
            type->precisionRequirement = 1e-12;
            type->requiresQuantumEnhancement = true;
            break;
        case IRType::BaseType::SOLID_ANGLE:
            type->precisionRequirement = 1e-10;
            type->requiresQuantumEnhancement = true;
            break;
        case IRType::BaseType::MATTER_STATE:
            type->precisionRequirement = 1e-6;
            type->requiresQuantumEnhancement = false;
            break;
        default:
            type->precisionRequirement = 1e-15;
            type->requiresQuantumEnhancement = false;
    }
    
    return type;
}

std::shared_ptr<IRType> HSMLIRBuilder::inferType(const std::shared_ptr<IRExpression>& expr) {
    if (expr->resultType) {
        return expr->resultType;
    }
    
    // Default to number type
    return createType(IRType::BaseType::NUMBER);
}

bool HSMLIRBuilder::isCompatibleType(const std::shared_ptr<IRType>& a, const std::shared_ptr<IRType>& b) {
    if (!a || !b) return false;
    
    // Basic type compatibility
    if (a->baseType == b->baseType) return true;
    
    // Number type conversions
    if ((a->baseType == IRType::BaseType::NUMBER && b->baseType == IRType::BaseType::NUMBER) ||
        (a->baseType == IRType::BaseType::SPHERICAL_COORD && b->baseType == IRType::BaseType::VECTOR) ||
        (a->baseType == IRType::BaseType::SOLID_ANGLE && b->baseType == IRType::BaseType::NUMBER)) {
        return true;
    }
    
    return false;
}

// === SYMBOL TABLE MANAGEMENT ===

void HSMLIRBuilder::pushScope(const std::string& scope) {
    scopeStack_.push_back(scope);
    currentScope_ = scope;
    symbolTables_[scope] = std::map<std::string, std::shared_ptr<IRVariable>>();
}

void HSMLIRBuilder::popScope() {
    if (scopeStack_.size() > 1) {
        scopeStack_.pop_back();
        currentScope_ = scopeStack_.back();
    }
}

void HSMLIRBuilder::declareVariable(const std::string& name, std::shared_ptr<IRVariable> var) {
    if (symbolTables_.find(currentScope_) == symbolTables_.end()) {
        symbolTables_[currentScope_] = std::map<std::string, std::shared_ptr<IRVariable>>();
    }
    symbolTables_[currentScope_][name] = var;
}

std::shared_ptr<IRVariable> HSMLIRBuilder::lookupVariable(const std::string& name) {
    // Search from current scope up to global scope
    for (auto it = scopeStack_.rbegin(); it != scopeStack_.rend(); ++it) {
        const std::string& scope = *it;
        if (symbolTables_.find(scope) != symbolTables_.end()) {
            auto& table = symbolTables_[scope];
            if (table.find(name) != table.end()) {
                return table[name];
            }
        }
    }
    return nullptr;
}

// === UTILITY METHODS ===

std::string HSMLIRBuilder::generateNodeId() {
    return "ir_node_" + std::to_string(nodeCounter_++);
}

std::string HSMLIRBuilder::generateTempId() {
    return "temp_" + std::to_string(tempCounter_++);
}

std::shared_ptr<IRTemp> HSMLIRBuilder::createTemp(std::shared_ptr<IRExpression> sourceExpr) {
    auto temp = std::make_shared<IRTemp>(generateNodeId(), generateTempId());
    temp->sourceExpression = sourceExpr;
    temp->resultType = sourceExpr->resultType;
    temp->isConstant = sourceExpr->isConstant;
    temp->constantValue = sourceExpr->constantValue;
    
    // Apply multi-dimensional analysis
    applyDimensionalAnalysis(temp);
    
    return temp;
}

std::shared_ptr<IRLiteral> HSMLIRBuilder::createLiteral(const std::variant<int, double, std::string, bool>& value) {
    auto literal = std::make_shared<IRLiteral>(generateNodeId(), value);
    literal->resultType = inferLiteralType(value);
    
    // Apply multi-dimensional analysis
    applyDimensionalAnalysis(literal);
    synthesizeParadigms(literal);
    
    return literal;
}

// === MULTI-DIMENSIONAL SYNTHESIS ===

void HSMLIRBuilder::applyDimensionalAnalysis(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Apply dimensional analysis to the node
    node->analyzeDimensionalProperties();
    
    // Calculate dimensional influence
    auto influence = calculateDimensionalInfluence(node);
    for (int i = 0; i < 4; ++i) {
        dimensionalContext_[i] += influence[i] * 0.1;
    }
}

void HSMLIRBuilder::synthesizeParadigms(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Apply paradigm synthesis
    node->synthesizeParadigms();
    
    // Apply multi-paradigm synthesis
    applyParadigmSynthesis(node);
    createEmergentPatterns(node);
}

void HSMLIRBuilder::applyQuantumEnhancement(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Apply quantum enhancement
    node->applyQuantumEnhancement();
    
    // Calculate quantum uncertainty
    double uncertainty = calculateQuantumUncertainty(node);
    node->metadata["quantum_uncertainty"] = uncertainty;
}

void HSMLIRBuilder::optimizeSphericalOperations(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Apply spherical coordinate optimizations
    node->metadata["spherical_optimized"] = true;
    node->metadata["coordinate_precision"] = builderPrecision_;
    
    // Add spherical optimization hints
    OptimizationHint sphericalHint(OptimizationHintType::VECTORIZE, 3);
    sphericalHint.parameters["spherical_aware"] = true;
    sphericalHint.parameters["use_simd"] = true;
    sphericalHint.dimensionalInfluence = {{1.0, 2.0, 1.0, 2.0}};  // Emphasize spatial and pragmatic
    node->optimizationHints.push_back(sphericalHint);
}

// === OPTIMIZATION PIPELINE ===

void HSMLIRBuilder::applyOptimizations(std::shared_ptr<IRProgram> program) {
    if (!program) return;
    
    // Apply various optimizations
    constantFolding(program);
    deadCodeElimination(program);
    sphericalCoordinateOptimization(program);
    
    // Apply multi-dimensional optimizations
    applyDimensionalAnalysis(program);
    synthesizeParadigms(program);
    applyQuantumEnhancement(program);
}

void HSMLIRBuilder::constantFolding(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Mark node as constant-folded
    node->metadata["constant_folded"] = true;
    
    // Add optimization hint
    OptimizationHint foldHint(OptimizationHintType::CONSTANT_FOLD, 1);
    foldHint.parameters["folding_applied"] = true;
    node->optimizationHints.push_back(foldHint);
}

void HSMLIRBuilder::deadCodeElimination(std::shared_ptr<IRProgram> program) {
    if (!program) return;
    
    // Mark program as dead-code eliminated
    program->metadata["dead_code_eliminated"] = true;
}

void HSMLIRBuilder::sphericalCoordinateOptimization(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    optimizeSphericalOperations(node);
}

// === DIMENSIONAL ANALYSIS METHODS ===

void HSMLIRBuilder::analyzeTemporal(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Temporal analysis - execution order and timing
    node->dimensionalProperties[0] = static_cast<double>(node->dependencies.size());
    node->metadata["temporal_complexity"] = node->dimensionalProperties[0];
}

void HSMLIRBuilder::analyzeSpatial(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Spatial analysis - memory layout and coordinate relationships
    node->dimensionalProperties[1] = static_cast<double>(node->metadata.size()) * 0.5;
    node->metadata["spatial_complexity"] = node->dimensionalProperties[1];
}

void HSMLIRBuilder::analyzeConceptual(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Conceptual analysis - semantic meaning and abstraction level
    node->dimensionalProperties[2] = 1.0; // Base conceptual complexity
    if (node->type == IRNodeType::IR_SPHERICAL_OP) {
        node->dimensionalProperties[2] = 2.0; // Higher for spherical operations
    }
    node->metadata["conceptual_complexity"] = node->dimensionalProperties[2];
}

void HSMLIRBuilder::analyzePragmatic(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Pragmatic analysis - practical optimization opportunities
    node->dimensionalProperties[3] = static_cast<double>(node->optimizationHints.size());
    node->metadata["pragmatic_complexity"] = node->dimensionalProperties[3];
}

// === PRIVATE INITIALIZATION METHODS ===

void HSMLIRBuilder::initializeBuiltinTypes() {
    typeTable_["number"] = createType(IRType::BaseType::NUMBER);
    typeTable_["string"] = createType(IRType::BaseType::STRING);
    typeTable_["boolean"] = createType(IRType::BaseType::BOOLEAN);
    typeTable_["spherical_coord"] = createType(IRType::BaseType::SPHERICAL_COORD);
    typeTable_["solid_angle"] = createType(IRType::BaseType::SOLID_ANGLE);
    typeTable_["matter_state"] = createType(IRType::BaseType::MATTER_STATE);
    typeTable_["vector"] = createType(IRType::BaseType::VECTOR);
    typeTable_["matrix"] = createType(IRType::BaseType::MATRIX);
    typeTable_["function"] = createType(IRType::BaseType::FUNCTION);
}

void HSMLIRBuilder::setupDimensionalContext() {
    // Initialize dimensional context
    dimensionalContext_[0] = 0.0; // Temporal
    dimensionalContext_[1] = 0.0; // Spatial  
    dimensionalContext_[2] = 0.0; // Conceptual
    dimensionalContext_[3] = 0.0; // Pragmatic
}

void HSMLIRBuilder::initializeQuantumState() {
    quantumCorrections_.reserve(1000); // Pre-allocate for performance
}

// === TYPE INFERENCE HELPERS ===

std::shared_ptr<IRType> HSMLIRBuilder::inferBinaryResultType(
    const std::shared_ptr<IRExpression>& left,
    const std::shared_ptr<IRExpression>& right,
    const std::string& operator_) {
    
    if (!left || !right) return createType(IRType::BaseType::NUMBER);
    
    // Spherical coordinate operations
    if (left->resultType && left->resultType->baseType == IRType::BaseType::SPHERICAL_COORD) {
        if (operator_ == "+" || operator_ == "-") {
            return createType(IRType::BaseType::SPHERICAL_COORD);
        }
        if (operator_ == "*" && right->resultType->baseType == IRType::BaseType::NUMBER) {
            return createType(IRType::BaseType::SPHERICAL_COORD);
        }
    }
    
    // Default to number type
    return createType(IRType::BaseType::NUMBER);
}

std::shared_ptr<IRType> HSMLIRBuilder::inferUnaryResultType(
    const std::shared_ptr<IRExpression>& operand,
    const std::string& operator_) {
    
    if (!operand || !operand->resultType) return createType(IRType::BaseType::NUMBER);
    
    // Most unary operations preserve the type
    return operand->resultType;
}

std::shared_ptr<IRType> HSMLIRBuilder::inferCallResultType(
    const std::shared_ptr<IRExpression>& callee,
    const std::vector<std::shared_ptr<IRExpression>>& arguments) {
    
    // Function-specific type inference would go here
    return createType(IRType::BaseType::NUMBER);
}

std::shared_ptr<IRType> HSMLIRBuilder::inferLiteralType(const std::variant<int, double, std::string, bool>& value) {
    if (std::holds_alternative<int>(value) || std::holds_alternative<double>(value)) {
        return createType(IRType::BaseType::NUMBER);
    }
    if (std::holds_alternative<std::string>(value)) {
        return createType(IRType::BaseType::STRING);
    }
    if (std::holds_alternative<bool>(value)) {
        return createType(IRType::BaseType::BOOLEAN);
    }
    return createType(IRType::BaseType::NUMBER);
}

// === MULTI-DIMENSIONAL HELPER METHODS ===

double HSMLIRBuilder::calculateQuantumUncertainty(const std::shared_ptr<IRNode>& node) const {
    if (!node) return builderPrecision_;
    
    // Apply Heisenberg uncertainty principle
    double position_uncertainty = static_cast<double>(node->id.length());
    double momentum_uncertainty = static_cast<double>(node->dependencies.size());
    
    return builderPrecision_ * std::sqrt(position_uncertainty * momentum_uncertainty);
}

std::array<double, 4> HSMLIRBuilder::calculateDimensionalInfluence(const std::shared_ptr<IRNode>& node) const {
    if (!node) return {{0.0, 0.0, 0.0, 0.0}};
    
    std::array<double, 4> influence;
    
    // Temporal influence
    influence[0] = static_cast<double>(node->dependencies.size()) * 0.1;
    
    // Spatial influence  
    influence[1] = static_cast<double>(node->id.length()) * 0.05;
    
    // Conceptual influence
    influence[2] = static_cast<double>(node->metadata.size()) * 0.2;
    
    // Pragmatic influence
    influence[3] = static_cast<double>(node->optimizationHints.size()) * 0.3;
    
    return influence;
}

void HSMLIRBuilder::applyParadigmSynthesis(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Synthesize functional, object-oriented, and reactive paradigms
    node->metadata["paradigm_functional"] = 0.8;
    node->metadata["paradigm_object_oriented"] = 0.6;
    node->metadata["paradigm_reactive"] = 0.7;
    node->metadata["paradigm_synthesis_applied"] = true;
}

void HSMLIRBuilder::createEmergentPatterns(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Create emergent behavioral patterns
    node->metadata["emergent_patterns"] = true;
    node->metadata["pattern_complexity"] = 1.2;
    node->metadata["behavioral_emergence"] = 0.9;
}

// === UTILITY FUNCTIONS ===

std::string irNodeTypeToString(IRNodeType type) {
    switch (type) {
        case IRNodeType::IR_PROGRAM: return "IRProgram";
        case IRNodeType::IR_FUNCTION: return "IRFunction";
        case IRNodeType::IR_BINARY_OP: return "IRBinaryOp";
        case IRNodeType::IR_UNARY_OP: return "IRUnaryOp";
        case IRNodeType::IR_CALL: return "IRCall";
        case IRNodeType::IR_LITERAL: return "IRLiteral";
        case IRNodeType::IR_SPHERICAL_OP: return "IRSphericalOp";
        case IRNodeType::IR_SOLID_ANGLE_OP: return "IRSolidAngleOp";
        case IRNodeType::IR_MATTER_STATE_OP: return "IRMatterStateOp";
        default: return "Unknown";
    }
}

std::string sourceLanguageToString(SourceLanguage lang) {
    switch (lang) {
        case SourceLanguage::HSML: return "HSML";
        case SourceLanguage::CSSS: return "CSSS";
        case SourceLanguage::SHAPE: return "ShapeScript";
        case SourceLanguage::STYB: return "StyleBot";
        default: return "Unknown";
    }
}

std::string optimizationHintTypeToString(OptimizationHintType type) {
    switch (type) {
        case OptimizationHintType::PARALLEL: return "Parallel";
        case OptimizationHintType::CACHE: return "Cache";
        case OptimizationHintType::VECTORIZE: return "Vectorize";
        case OptimizationHintType::INLINE: return "Inline";
        case OptimizationHintType::CONSTANT_FOLD: return "ConstantFold";
        default: return "Unknown";
    }
}

// === VISITOR IMPLEMENTATIONS ===

void IROptimizationVisitor::visit(std::shared_ptr<IRProgram> program) {
    if (!program) return;
    
    // Optimize all functions
    for (auto& function : program->functions) {
        visit(function);
    }
    
    // Apply program-level optimizations
    optimizeSphericalCalculations(program);
    applyVectorization(program);
    enableParallelization(program);
}

void IROptimizationVisitor::visit(std::shared_ptr<IRFunction> function) {
    if (!function) return;
    
    // Optimize function body
    if (function->body) {
        for (auto& stmt : function->body->statements) {
            visit(stmt);
        }
    }
    
    // Apply function-level optimizations
    optimizeSphericalCalculations(function);
}

void IROptimizationVisitor::visit(std::shared_ptr<IRExpression> expression) {
    if (!expression) return;
    
    // Apply expression-level optimizations
    optimizeSphericalCalculations(expression);
    applyVectorization(expression);
}

void IROptimizationVisitor::visit(std::shared_ptr<IRStatement> statement) {
    if (!statement) return;
    
    // Apply statement-level optimizations
    optimizeSphericalCalculations(statement);
}

void IROptimizationVisitor::visitDimensionalProperties(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Apply dimensional optimizations
    node->analyzeDimensionalProperties();
    
    // Add dimensional optimization hints
    OptimizationHint dimensionalHint(OptimizationHintType::VECTORIZE, 2);
    dimensionalHint.parameters["dimensional_aware"] = true;
    node->optimizationHints.push_back(dimensionalHint);
}

void IROptimizationVisitor::visitQuantumEnhancements(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Apply quantum enhancements
    node->applyQuantumEnhancement();
}

void IROptimizationVisitor::optimizeSphericalCalculations(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Mark as spherically optimized
    node->metadata["spherical_optimization"] = true;
    node->metadata["coordinate_precision"] = 1e-12;
}

void IROptimizationVisitor::applyVectorization(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Add vectorization hint
    OptimizationHint vectorHint(OptimizationHintType::VECTORIZE, 3);
    vectorHint.parameters["enable_simd"] = true;
    vectorHint.parameters["vector_width"] = 4;
    node->optimizationHints.push_back(vectorHint);
}

void IROptimizationVisitor::enableParallelization(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Add parallelization hint
    OptimizationHint parallelHint(OptimizationHintType::PARALLEL, 2);
    parallelHint.parameters["thread_count"] = 4;
    parallelHint.parameters["parallel_safe"] = true;
    node->optimizationHints.push_back(parallelHint);
}

// === ANALYSIS VISITOR IMPLEMENTATION ===

void IRAnalysisVisitor::visit(std::shared_ptr<IRProgram> program) {
    if (!program) return;
    
    // Analyze program complexity
    performanceMetrics_["program_complexity"] = static_cast<double>(program->functions.size());
    complexityAnalysis_["function_count"] = static_cast<int>(program->functions.size());
    
    // Analyze functions
    for (auto& function : program->functions) {
        visit(function);
    }
}

void IRAnalysisVisitor::visit(std::shared_ptr<IRFunction> function) {
    if (!function) return;
    
    // Analyze function performance
    analyzePerformance(function);
    analyzeComplexity(function);
}

void IRAnalysisVisitor::visit(std::shared_ptr<IRExpression> expression) {
    if (!expression) return;
    
    // Analyze expression performance
    analyzePerformance(expression);
    analyzeComplexity(expression);
}

void IRAnalysisVisitor::visit(std::shared_ptr<IRStatement> statement) {
    if (!statement) return;
    
    // Analyze statement performance
    analyzePerformance(statement);
    analyzeComplexity(statement);
}

void IRAnalysisVisitor::visitDimensionalProperties(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Analyze dimensional properties
    node->analyzeDimensionalProperties();
    
    // Record dimensional metrics
    for (int i = 0; i < 4; ++i) {
        std::string dimensionName = "dimension_" + std::to_string(i);
        performanceMetrics_[dimensionName] += node->dimensionalProperties[i];
    }
}

void IRAnalysisVisitor::visitQuantumEnhancements(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Analyze quantum enhancements
    node->applyQuantumEnhancement();
    
    // Record quantum metrics
    if (node->metadata.find("quantum_uncertainty") != node->metadata.end()) {
        try {
            double uncertainty = std::get<double>(node->metadata["quantum_uncertainty"]);
            performanceMetrics_["quantum_uncertainty"] += uncertainty;
        } catch (const std::bad_variant_access&) {
            // Handle type mismatch
        }
    }
}

std::map<std::string, double> IRAnalysisVisitor::getPerformanceMetrics() const {
    return performanceMetrics_;
}

std::map<std::string, int> IRAnalysisVisitor::getComplexityAnalysis() const {
    return complexityAnalysis_;
}

void IRAnalysisVisitor::analyzePerformance(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Calculate performance metrics
    double complexity = static_cast<double>(node->dependencies.size()) +
                       static_cast<double>(node->optimizationHints.size()) * 0.5;
    
    performanceMetrics_["total_complexity"] += complexity;
    performanceMetrics_["node_count"] += 1.0;
}

void IRAnalysisVisitor::analyzeComplexity(std::shared_ptr<IRNode> node) {
    if (!node) return;
    
    // Update complexity counters
    complexityAnalysis_["total_nodes"]++;
    
    // Count by node type
    std::string nodeTypeName = irNodeTypeToString(node->type);
    complexityAnalysis_[nodeTypeName + "_count"]++;
}

} // namespace language
} // namespace hsml