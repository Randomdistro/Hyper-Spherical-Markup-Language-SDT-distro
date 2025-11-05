/**
 * HSML Intermediate Representation - C++ Implementation
 * Unified IR for all four HSML languages with spherical coordinate optimization
 * You are now all the C
 */

#pragma once

#include <string>
#include <vector>
#include <map>
#include <set>
#include <memory>
#include <variant>
#include <optional>
#include <unordered_map>

namespace hsml {
namespace language {

// Forward declarations
class HSMLIRBuilder;
class IRNode;
class IRExpression;

// IR Node Types with multi-dimensional classification
enum class IRNodeType {
    // === PROGRAM STRUCTURE ===
    IR_PROGRAM,
    IR_MODULE,
    IR_FUNCTION,
    IR_BLOCK,
    
    // === EXPRESSIONS ===
    IR_BINARY_OP,
    IR_UNARY_OP,
    IR_CALL,
    IR_LOAD,
    IR_STORE,
    IR_LITERAL,
    IR_TEMP,
    
    // === CONTROL FLOW ===
    IR_IF,
    IR_WHILE,
    IR_FOR,
    IR_RETURN,
    IR_BREAK,
    IR_CONTINUE,
    
    // === HSML-SPECIFIC ===
    IR_HSML_ELEMENT,
    IR_HSML_ATTRIBUTE,
    IR_HSML_RENDER,
    
    // === CSSS-SPECIFIC ===
    IR_CSSS_RULE,
    IR_CSSS_DECLARATION,
    IR_CSSS_KEYFRAME,
    IR_CSSS_MATERIAL,
    IR_CSSS_ANIMATION,
    
    // === SHAPESCRIPT-SPECIFIC ===
    IR_SHAPE_PHYSICS,
    IR_SHAPE_FORCE,
    IR_SHAPE_CONSTRAINT,
    IR_SHAPE_EVENT,
    
    // === STYLEBOT-SPECIFIC ===
    IR_STYB_BOT,
    IR_STYB_AGENT,
    IR_STYB_RENDER,
    IR_STYB_OPTIMIZE,
    
    // === SPHERICAL COORDINATE OPERATIONS ===
    IR_SPHERICAL_OP,
    IR_SOLID_ANGLE_OP,
    IR_MATTER_STATE_OP,
    
    // === OPTIMIZATION HINTS ===
    IR_OPTIMIZE_HINT,
    IR_PARALLEL_HINT,
    IR_CACHE_HINT
};

// Source language enumeration
enum class SourceLanguage {
    HSML,
    CSSS,
    SHAPE,
    STYB
};

// Optimization hint types
enum class OptimizationHintType {
    PARALLEL,
    CACHE,
    VECTORIZE,
    INLINE,
    CONSTANT_FOLD
};

// Optimization hint structure with multi-dimensional parameters
struct OptimizationHint {
    OptimizationHintType type;
    int priority;
    std::map<std::string, std::variant<int, double, std::string, bool>> parameters;
    
    // Multi-dimensional hint properties
    std::array<double, 4> dimensionalInfluence; // temporal, spatial, conceptual, pragmatic
    
    OptimizationHint(OptimizationHintType t, int p)
        : type(t), priority(p), dimensionalInfluence{{1.0, 1.0, 1.0, 1.0}} {}
};

// Base IR Type system with dimensional awareness
struct IRType {
    enum class BaseType {
        NUMBER,
        STRING,
        BOOLEAN,
        SPHERICAL_COORD,
        SOLID_ANGLE,
        MATTER_STATE,
        VECTOR,
        MATRIX,
        FUNCTION
    };
    
    BaseType baseType;
    std::vector<int> dimensions;
    std::map<std::string, std::shared_ptr<IRType>> properties;
    std::map<std::string, std::variant<int, double, std::string>> constraints;
    
    // Multi-dimensional type properties
    double precisionRequirement;
    bool requiresQuantumEnhancement;
    
    IRType(BaseType bt) 
        : baseType(bt), precisionRequirement(1e-15), requiresQuantumEnhancement(false) {}
};

// Base IR Node interface with multi-dimensional synthesis
class IRNode {
public:
    IRNodeType type;
    std::string id;
    SourceLanguage sourceLanguage;
    std::map<std::string, std::variant<int, double, std::string, bool>> metadata;
    std::set<std::string> dependencies;
    std::vector<OptimizationHint> optimizationHints;
    
    // Multi-dimensional node properties
    std::array<double, 4> dimensionalProperties; // temporal, spatial, conceptual, pragmatic
    
    IRNode(IRNodeType t, const std::string& nodeId, SourceLanguage lang)
        : type(t), id(nodeId), sourceLanguage(lang)
        , dimensionalProperties{{0.0, 0.0, 0.0, 0.0}} {}
    
    virtual ~IRNode() = default;
    
    // Multi-dimensional analysis methods
    virtual void analyzeDimensionalProperties();
    virtual void applyQuantumEnhancement();
    virtual void synthesizeParadigms();
};

// Program IR - Root of the IR tree
class IRProgram : public IRNode {
public:
    std::vector<std::shared_ptr<class IRFunction>> functions;
    std::vector<std::shared_ptr<struct IRVariable>> globalVariables;
    std::vector<std::shared_ptr<class IRModule>> modules;
    std::string entryPoint;
    
    IRProgram(const std::string& nodeId)
        : IRNode(IRNodeType::IR_PROGRAM, nodeId, SourceLanguage::HSML)
        , entryPoint("main") {}
};

// Module IR - Code organization unit
class IRModule : public IRNode {
public:
    std::string name;
    std::vector<std::string> imports;
    std::vector<std::string> exports;
    std::vector<std::shared_ptr<class IRFunction>> functions;
    std::vector<std::shared_ptr<struct IRVariable>> variables;
    
    IRModule(const std::string& nodeId, const std::string& moduleName)
        : IRNode(IRNodeType::IR_MODULE, nodeId, SourceLanguage::HSML)
        , name(moduleName) {}
};

// Variable representation
struct IRVariable {
    std::string name;
    std::shared_ptr<IRType> type;
    bool isConstant;
    std::shared_ptr<class IRExpression> initialValue;
    std::string scope;
    
    // Multi-dimensional variable properties
    double quantumUncertainty;
    bool requiresSphericalOptimization;
    
    IRVariable(const std::string& varName, std::shared_ptr<IRType> varType)
        : name(varName), type(varType), isConstant(false)
        , quantumUncertainty(1e-15), requiresSphericalOptimization(false) {}
};

// Parameter representation
struct IRParameter {
    std::string name;
    std::shared_ptr<IRType> type;
    std::shared_ptr<class IRExpression> defaultValue;
    
    IRParameter(const std::string& paramName, std::shared_ptr<IRType> paramType)
        : name(paramName), type(paramType) {}
};

// Block IR - Scoped statement container
class IRBlock : public IRNode {
public:
    std::vector<std::shared_ptr<class IRStatement>> statements;
    std::vector<std::shared_ptr<IRVariable>> variables;
    std::string scope;
    
    IRBlock(const std::string& nodeId, const std::string& blockScope)
        : IRNode(IRNodeType::IR_BLOCK, nodeId, SourceLanguage::HSML)
        , scope(blockScope) {}
};

// Function IR - Executable code unit
class IRFunction : public IRNode {
public:
    std::string name;
    std::vector<IRParameter> parameters;
    std::shared_ptr<IRType> returnType;
    std::shared_ptr<IRBlock> body;
    bool isPure;
    std::set<std::string> sideEffects;
    
    // Multi-dimensional function properties
    double computationalComplexity;
    bool requiresParallelization;
    
    IRFunction(const std::string& nodeId, const std::string& funcName)
        : IRNode(IRNodeType::IR_FUNCTION, nodeId, SourceLanguage::HSML)
        , name(funcName), isPure(false)
        , computationalComplexity(1.0), requiresParallelization(false) {}
};

// Expression IR - Base class for all expressions
class IRExpression : public IRNode {
public:
    std::shared_ptr<IRType> resultType;
    bool isConstant;
    std::variant<int, double, std::string, bool> constantValue;
    
    // Multi-dimensional expression properties
    double precisionRequirement;
    bool requiresQuantumCorrection;
    
    IRExpression(IRNodeType exprType, const std::string& nodeId, SourceLanguage lang)
        : IRNode(exprType, nodeId, lang), isConstant(false)
        , precisionRequirement(1e-15), requiresQuantumCorrection(false) {}
};

// Binary operation IR
class IRBinaryOp : public IRExpression {
public:
    std::string operator_;
    std::shared_ptr<IRExpression> left;
    std::shared_ptr<IRExpression> right;
    std::shared_ptr<struct SphericalOptimization> sphericalOptimization;
    
    IRBinaryOp(const std::string& nodeId, const std::string& op)
        : IRExpression(IRNodeType::IR_BINARY_OP, nodeId, SourceLanguage::HSML)
        , operator_(op) {}
};

// Unary operation IR
class IRUnaryOp : public IRExpression {
public:
    std::string operator_;
    std::shared_ptr<IRExpression> operand;
    
    IRUnaryOp(const std::string& nodeId, const std::string& op)
        : IRExpression(IRNodeType::IR_UNARY_OP, nodeId, SourceLanguage::HSML)
        , operator_(op) {}
};

// Function call IR
class IRCall : public IRExpression {
public:
    std::string functionName;
    std::vector<std::shared_ptr<IRExpression>> arguments;
    bool isBuiltin;
    std::shared_ptr<struct SphericalOptimization> sphericalOptimization;
    
    IRCall(const std::string& nodeId, const std::string& funcName)
        : IRExpression(IRNodeType::IR_CALL, nodeId, SourceLanguage::HSML)
        , functionName(funcName), isBuiltin(false) {}
};

// Load operation IR
class IRLoad : public IRExpression {
public:
    std::string variable;
    std::shared_ptr<IRExpression> index;
    
    IRLoad(const std::string& nodeId, const std::string& varName)
        : IRExpression(IRNodeType::IR_LOAD, nodeId, SourceLanguage::HSML)
        , variable(varName) {}
};

// Store operation IR
class IRStore : public IRNode {
public:
    std::string variable;
    std::shared_ptr<IRExpression> value;
    std::shared_ptr<IRExpression> index;
    
    IRStore(const std::string& nodeId, const std::string& varName)
        : IRNode(IRNodeType::IR_STORE, nodeId, SourceLanguage::HSML)
        , variable(varName) {}
};

// Literal IR
class IRLiteral : public IRExpression {
public:
    std::variant<int, double, std::string, bool> value;
    std::shared_ptr<struct SphericalCoordinateIR> sphericalCoordinate;
    std::shared_ptr<struct SolidAngleIR> solidAngle;
    std::shared_ptr<struct MatterStateIR> matterState;
    
    IRLiteral(const std::string& nodeId, const std::variant<int, double, std::string, bool>& val)
        : IRExpression(IRNodeType::IR_LITERAL, nodeId, SourceLanguage::HSML)
        , value(val) {
        isConstant = true;
        constantValue = val;
    }
};

// Temporary variable IR
class IRTemp : public IRExpression {
public:
    std::string tempId;
    std::shared_ptr<IRExpression> sourceExpression;
    
    IRTemp(const std::string& nodeId, const std::string& tempVar)
        : IRExpression(IRNodeType::IR_TEMP, nodeId, SourceLanguage::HSML)
        , tempId(tempVar) {}
};

// Statement IR - Base class for all statements
class IRStatement : public IRNode {
public:
    IRStatement(IRNodeType stmtType, const std::string& nodeId, SourceLanguage lang)
        : IRNode(stmtType, nodeId, lang) {}
};

// Control flow statements
class IRIf : public IRStatement {
public:
    std::shared_ptr<IRExpression> condition;
    std::shared_ptr<IRBlock> consequent;
    std::shared_ptr<IRBlock> alternate;
    
    IRIf(const std::string& nodeId)
        : IRStatement(IRNodeType::IR_IF, nodeId, SourceLanguage::HSML) {}
};

class IRWhile : public IRStatement {
public:
    std::shared_ptr<IRExpression> condition;
    std::shared_ptr<IRBlock> body;
    
    IRWhile(const std::string& nodeId)
        : IRStatement(IRNodeType::IR_WHILE, nodeId, SourceLanguage::HSML) {}
};

class IRFor : public IRStatement {
public:
    std::shared_ptr<IRStatement> init;
    std::shared_ptr<IRExpression> condition;
    std::shared_ptr<IRStatement> update;
    std::shared_ptr<IRBlock> body;
    
    IRFor(const std::string& nodeId)
        : IRStatement(IRNodeType::IR_FOR, nodeId, SourceLanguage::HSML) {}
};

class IRReturn : public IRStatement {
public:
    std::shared_ptr<IRExpression> value;
    
    IRReturn(const std::string& nodeId)
        : IRStatement(IRNodeType::IR_RETURN, nodeId, SourceLanguage::HSML) {}
};

// === HSML-SPECIFIC IR NODES ===

class IRHSMLElement : public IRNode {
public:
    std::string tagName;
    std::vector<std::shared_ptr<class IRHSMLAttribute>> attributes;
    std::vector<std::shared_ptr<IRHSMLElement>> children;
    std::vector<std::shared_ptr<class IRHSMLRender>> renderInstructions;
    std::shared_ptr<struct SphericalCoordinateIR> sphericalPosition;
    std::string material;
    std::string behavior;
    
    IRHSMLElement(const std::string& nodeId, const std::string& tag)
        : IRNode(IRNodeType::IR_HSML_ELEMENT, nodeId, SourceLanguage::HSML)
        , tagName(tag) {}
};

class IRHSMLAttribute : public IRNode {
public:
    std::string name;
    std::shared_ptr<IRExpression> value;
    bool isSphericalCoordinate;
    
    IRHSMLAttribute(const std::string& nodeId, const std::string& attrName)
        : IRNode(IRNodeType::IR_HSML_ATTRIBUTE, nodeId, SourceLanguage::HSML)
        , name(attrName), isSphericalCoordinate(false) {}
};

class IRHSMLRender : public IRNode {
public:
    std::string target;
    std::shared_ptr<IRExpression> quality;
    int priority;
    std::shared_ptr<struct SphericalOptimization> sphericalOptimization;
    
    IRHSMLRender(const std::string& nodeId, const std::string& renderTarget)
        : IRNode(IRNodeType::IR_HSML_RENDER, nodeId, SourceLanguage::HSML)
        , target(renderTarget), priority(1) {}
};

// === CSSS-SPECIFIC IR NODES ===

class IRCSSSRule : public IRNode {
public:
    std::vector<std::string> selectors;
    std::vector<std::shared_ptr<class IRCSSSDeclaration>> declarations;
    std::shared_ptr<class IRCSSSMaterial> material;
    std::shared_ptr<class IRCSSSAnimation> animation;
    
    IRCSSSRule(const std::string& nodeId)
        : IRNode(IRNodeType::IR_CSSS_RULE, nodeId, SourceLanguage::CSSS) {}
};

class IRCSSSDeclaration : public IRNode {
public:
    std::string property;
    std::shared_ptr<IRExpression> value;
    bool important;
    
    IRCSSSDeclaration(const std::string& nodeId, const std::string& prop)
        : IRNode(IRNodeType::IR_CSSS_DECLARATION, nodeId, SourceLanguage::CSSS)
        , property(prop), important(false) {}
};

class IRCSSSMaterial : public IRNode {
public:
    std::string name;
    std::map<std::string, std::shared_ptr<IRExpression>> properties;
    std::shared_ptr<struct MatterStateIR> matterState;
    std::shared_ptr<struct SphericalOptimization> sphericalOptimization;
    
    IRCSSSMaterial(const std::string& nodeId, const std::string& matName)
        : IRNode(IRNodeType::IR_CSSS_MATERIAL, nodeId, SourceLanguage::CSSS)
        , name(matName) {}
};

class IRCSSSAnimation : public IRNode {
public:
    std::string name;
    std::vector<std::shared_ptr<class IRCSSSKeyframe>> keyframes;
    std::shared_ptr<IRExpression> duration;
    std::string easing;
    
    IRCSSSAnimation(const std::string& nodeId, const std::string& animName)
        : IRNode(IRNodeType::IR_CSSS_ANIMATION, nodeId, SourceLanguage::CSSS)
        , name(animName) {}
};

class IRCSSSKeyframe : public IRNode {
public:
    double percentage;
    std::vector<std::shared_ptr<IRCSSSDeclaration>> declarations;
    
    IRCSSSKeyframe(const std::string& nodeId, double pct)
        : IRNode(IRNodeType::IR_CSSS_KEYFRAME, nodeId, SourceLanguage::CSSS)
        , percentage(pct) {}
};

// === SHAPESCRIPT-SPECIFIC IR NODES ===

class IRShapePhysics : public IRNode {
public:
    std::vector<std::shared_ptr<class IRShapeForce>> forces;
    std::vector<std::shared_ptr<class IRShapeConstraint>> constraints;
    std::vector<std::shared_ptr<class IRShapeEvent>> events;
    std::shared_ptr<struct MatterStateIR> matterState;
    std::shared_ptr<struct SphericalOptimization> sphericalOptimization;
    
    IRShapePhysics(const std::string& nodeId)
        : IRNode(IRNodeType::IR_SHAPE_PHYSICS, nodeId, SourceLanguage::SHAPE) {}
};

class IRShapeForce : public IRNode {
public:
    enum class ForceType {
        GRAVITY,
        ELASTIC,
        VISCOUS,
        ELECTROMAGNETIC
    };
    
    ForceType forceType;
    std::shared_ptr<IRExpression> magnitude;
    std::shared_ptr<struct SphericalCoordinateIR> direction;
    std::shared_ptr<struct SphericalOptimization> sphericalOptimization;
    
    IRShapeForce(const std::string& nodeId, ForceType type)
        : IRNode(IRNodeType::IR_SHAPE_FORCE, nodeId, SourceLanguage::SHAPE)
        , forceType(type) {}
};

class IRShapeConstraint : public IRNode {
public:
    enum class ConstraintType {
        SPHERICAL_SURFACE,
        RADIAL_RANGE,
        ANGULAR_CONE
    };
    
    ConstraintType constraintType;
    std::map<std::string, std::shared_ptr<IRExpression>> parameters;
    
    IRShapeConstraint(const std::string& nodeId, ConstraintType type)
        : IRNode(IRNodeType::IR_SHAPE_CONSTRAINT, nodeId, SourceLanguage::SHAPE)
        , constraintType(type) {}
};

class IRShapeEvent : public IRNode {
public:
    std::string trigger;
    std::shared_ptr<IRExpression> response;
    std::vector<std::shared_ptr<IRExpression>> conditions;
    
    IRShapeEvent(const std::string& nodeId, const std::string& eventTrigger)
        : IRNode(IRNodeType::IR_SHAPE_EVENT, nodeId, SourceLanguage::SHAPE)
        , trigger(eventTrigger) {}
};

// === STYLEBOT-SPECIFIC IR NODES ===

class IRStyleBot : public IRNode {
public:
    std::string name;
    std::vector<std::shared_ptr<class IRStyleBotAgent>> agents;
    bool parallel;
    std::vector<std::shared_ptr<class IRStyleBotRender>> renderTasks;
    std::vector<std::shared_ptr<class IRStyleBotOptimize>> optimizeTasks;
    
    IRStyleBot(const std::string& nodeId, const std::string& botName)
        : IRNode(IRNodeType::IR_STYB_BOT, nodeId, SourceLanguage::STYB)
        , name(botName), parallel(false) {}
};

class IRStyleBotAgent : public IRNode {
public:
    std::string name;
    std::vector<std::shared_ptr<class IRStyleBotRender>> renderTasks;
    std::vector<std::shared_ptr<class IRStyleBotOptimize>> optimizeTasks;
    bool parallel;
    
    IRStyleBotAgent(const std::string& nodeId, const std::string& agentName)
        : IRNode(IRNodeType::IR_STYB_AGENT, nodeId, SourceLanguage::STYB)
        , name(agentName), parallel(false) {}
};

class IRStyleBotRender : public IRNode {
public:
    std::string target;
    std::shared_ptr<IRExpression> quality;
    int priority;
    std::shared_ptr<struct SphericalOptimization> sphericalOptimization;
    
    IRStyleBotRender(const std::string& nodeId, const std::string& renderTarget)
        : IRNode(IRNodeType::IR_STYB_RENDER, nodeId, SourceLanguage::STYB)
        , target(renderTarget), priority(1) {}
};

class IRStyleBotOptimize : public IRNode {
public:
    enum class OptimizationType {
        PERFORMANCE,
        MEMORY,
        QUALITY,
        SPHERICAL
    };
    
    OptimizationType optimizationType;
    std::map<std::string, std::shared_ptr<IRExpression>> parameters;
    std::string target;
    
    IRStyleBotOptimize(const std::string& nodeId, OptimizationType type)
        : IRNode(IRNodeType::IR_STYB_OPTIMIZE, nodeId, SourceLanguage::STYB)
        , optimizationType(type) {}
};

// === SPECIALIZED IR STRUCTURES ===

// Spherical coordinate IR with quantum enhancement
struct SphericalCoordinateIR {
    std::shared_ptr<IRExpression> r;
    std::shared_ptr<IRExpression> theta;
    std::shared_ptr<IRExpression> phi;
    bool isNormalized;
    std::vector<OptimizationHint> optimizationHints;
    
    // Multi-dimensional spherical properties
    double quantumUncertainty;
    std::array<double, 4> dimensionalCorrection;
    
    SphericalCoordinateIR()
        : isNormalized(false), quantumUncertainty(1e-15)
        , dimensionalCorrection{{0.0, 0.0, 0.0, 0.0}} {}
};

// Solid angle IR with dimensional synthesis
struct SolidAngleIR {
    std::shared_ptr<IRExpression> omega;
    std::shared_ptr<IRExpression> theta_min;
    std::shared_ptr<IRExpression> theta_max;
    std::shared_ptr<IRExpression> phi_min;
    std::shared_ptr<IRExpression> phi_max;
    bool isNormalized;
    
    // Multi-dimensional solid angle properties
    std::array<double, 4> steradianDistribution;
    double angularPrecision;
    
    SolidAngleIR()
        : isNormalized(false), angularPrecision(1e-12)
        , steradianDistribution{{1.0, 1.0, 1.0, 1.0}} {}
};

// Matter state IR with physics integration
struct MatterStateIR {
    enum class State {
        SOLID,
        LIQUID,
        GAS,
        PLASMA
    };
    
    State state;
    std::map<std::string, std::shared_ptr<IRExpression>> properties;
    std::map<std::string, std::shared_ptr<IRExpression>> physicsProperties;
    
    // Multi-dimensional matter properties
    std::array<double, 4> quantumState;
    double phaseTransitionThreshold;
    
    MatterStateIR(State s = State::SOLID)
        : state(s), phaseTransitionThreshold(273.15)
        , quantumState{{1.0, 0.0, 0.0, 0.0}} {}
};

// Spherical optimization configuration
struct SphericalOptimization {
    bool usePureSpherical;
    bool cacheResults;
    bool vectorize;
    bool parallelize;
    
    enum class Precision {
        SINGLE,
        DOUBLE,
        EXTENDED
    };
    
    Precision precision;
    int optimizationLevel;
    
    // Multi-dimensional optimization parameters
    std::array<double, 4> dimensionalWeights;
    double quantumEnhancementFactor;
    
    SphericalOptimization()
        : usePureSpherical(true), cacheResults(true), vectorize(true), parallelize(false)
        , precision(Precision::DOUBLE), optimizationLevel(2)
        , dimensionalWeights{{1.0, 1.0, 1.0, 1.0}}, quantumEnhancementFactor(1e-12) {}
};

/**
 * Multi-Dimensional HSML IR Builder
 * Implements the Multiplicitous Encephalopoidal approach to IR construction
 * Synthesizes multiple paradigms for comprehensive code representation
 */
class HSMLIRBuilder {
private:
    int tempCounter_;
    int nodeCounter_;
    std::string currentScope_;
    std::vector<std::string> scopeStack_;
    
    // Symbol tables for different scopes
    std::map<std::string, std::map<std::string, std::shared_ptr<IRVariable>>> symbolTables_;
    std::map<std::string, std::shared_ptr<IRFunction>> functionTable_;
    std::map<std::string, std::shared_ptr<IRType>> typeTable_;
    
    // Multi-dimensional builder state
    std::array<double, 4> dimensionalContext_;
    std::vector<std::string> quantumCorrections_;
    double builderPrecision_;
    
public:
    HSMLIRBuilder();
    ~HSMLIRBuilder() = default;
    
    // === MAIN IR CONSTRUCTION INTERFACE ===
    std::shared_ptr<IRProgram> buildIR(const void* ast); // void* for AST flexibility
    
    // === SPECIALIZED BUILDERS ===
    std::shared_ptr<IRHSMLElement> buildHSMLElement(const void* element);
    std::shared_ptr<IRCSSSRule> buildCSSSRule(const void* rule);
    std::shared_ptr<IRShapePhysics> buildShapePhysics(const void* physics);
    std::shared_ptr<IRStyleBot> buildStyleBot(const void* bot);
    
    // === EXPRESSION BUILDERS ===
    std::shared_ptr<IRExpression> buildExpression(const void* expr);
    std::shared_ptr<IRBinaryOp> buildBinaryOperation(const void* binOp);
    std::shared_ptr<IRUnaryOp> buildUnaryOperation(const void* unOp);
    std::shared_ptr<IRCall> buildCallExpression(const void* call);
    std::shared_ptr<IRLiteral> buildLiteral(const void* literal);
    
    // === SPECIALIZED COORDINATE BUILDERS ===
    std::shared_ptr<SphericalCoordinateIR> buildSphericalCoordinate(const void* coord);
    std::shared_ptr<SolidAngleIR> buildSolidAngle(const void* angle);
    std::shared_ptr<MatterStateIR> buildMatterState(const void* state);
    
    // === TYPE SYSTEM ===
    std::shared_ptr<IRType> createType(IRType::BaseType baseType);
    std::shared_ptr<IRType> inferType(const std::shared_ptr<IRExpression>& expr);
    bool isCompatibleType(const std::shared_ptr<IRType>& a, const std::shared_ptr<IRType>& b);
    
    // === SYMBOL TABLE MANAGEMENT ===
    void pushScope(const std::string& scope);
    void popScope();
    void declareVariable(const std::string& name, std::shared_ptr<IRVariable> var);
    std::shared_ptr<IRVariable> lookupVariable(const std::string& name);
    
    // === UTILITY METHODS ===
    std::string generateNodeId();
    std::string generateTempId();
    std::shared_ptr<IRTemp> createTemp(std::shared_ptr<IRExpression> sourceExpr);
    std::shared_ptr<IRLiteral> createLiteral(const std::variant<int, double, std::string, bool>& value);
    
    // === MULTI-DIMENSIONAL SYNTHESIS ===
    void applyDimensionalAnalysis(std::shared_ptr<IRNode> node);
    void synthesizeParadigms(std::shared_ptr<IRNode> node);
    void applyQuantumEnhancement(std::shared_ptr<IRNode> node);
    void optimizeSphericalOperations(std::shared_ptr<IRNode> node);
    
    // === OPTIMIZATION PIPELINE ===
    void applyOptimizations(std::shared_ptr<IRProgram> program);
    void constantFolding(std::shared_ptr<IRNode> node);
    void deadCodeElimination(std::shared_ptr<IRProgram> program);
    void sphericalCoordinateOptimization(std::shared_ptr<IRNode> node);
    
    // === DIMENSIONAL ANALYSIS ===
    void analyzeTemporal(std::shared_ptr<IRNode> node);
    void analyzeSpatial(std::shared_ptr<IRNode> node);
    void analyzeConceptual(std::shared_ptr<IRNode> node);
    void analyzePragmatic(std::shared_ptr<IRNode> node);
    
private:
    void initializeBuiltinTypes();
    void setupDimensionalContext();
    void initializeQuantumState();
    
    // Type inference helpers
    std::shared_ptr<IRType> inferBinaryResultType(
        const std::shared_ptr<IRExpression>& left,
        const std::shared_ptr<IRExpression>& right,
        const std::string& operator_);
    std::shared_ptr<IRType> inferUnaryResultType(
        const std::shared_ptr<IRExpression>& operand,
        const std::string& operator_);
    std::shared_ptr<IRType> inferCallResultType(
        const std::shared_ptr<IRExpression>& callee,
        const std::vector<std::shared_ptr<IRExpression>>& arguments);
    std::shared_ptr<IRType> inferLiteralType(const std::variant<int, double, std::string, bool>& value);
    
    // Multi-dimensional helper methods
    double calculateQuantumUncertainty(const std::shared_ptr<IRNode>& node) const;
    std::array<double, 4> calculateDimensionalInfluence(const std::shared_ptr<IRNode>& node) const;
    void applyParadigmSynthesis(std::shared_ptr<IRNode> node);
    void createEmergentPatterns(std::shared_ptr<IRNode> node);
};

// === UTILITY FUNCTIONS ===

std::string irNodeTypeToString(IRNodeType type);
std::string sourceLanguageToString(SourceLanguage lang);
std::string optimizationHintTypeToString(OptimizationHintType type);

// === IR VISITOR PATTERN ===

class IRVisitor {
public:
    virtual ~IRVisitor() = default;
    
    virtual void visit(std::shared_ptr<IRProgram> program) = 0;
    virtual void visit(std::shared_ptr<IRFunction> function) = 0;
    virtual void visit(std::shared_ptr<IRExpression> expression) = 0;
    virtual void visit(std::shared_ptr<IRStatement> statement) = 0;
    
    // Multi-dimensional visitation
    virtual void visitDimensionalProperties(std::shared_ptr<IRNode> node) = 0;
    virtual void visitQuantumEnhancements(std::shared_ptr<IRNode> node) = 0;
};

// === SPECIALIZED VISITORS ===

class IROptimizationVisitor : public IRVisitor {
public:
    void visit(std::shared_ptr<IRProgram> program) override;
    void visit(std::shared_ptr<IRFunction> function) override;
    void visit(std::shared_ptr<IRExpression> expression) override;
    void visit(std::shared_ptr<IRStatement> statement) override;
    void visitDimensionalProperties(std::shared_ptr<IRNode> node) override;
    void visitQuantumEnhancements(std::shared_ptr<IRNode> node) override;
    
private:
    void optimizeSphericalCalculations(std::shared_ptr<IRNode> node);
    void applyVectorization(std::shared_ptr<IRNode> node);
    void enableParallelization(std::shared_ptr<IRNode> node);
};

class IRAnalysisVisitor : public IRVisitor {
public:
    void visit(std::shared_ptr<IRProgram> program) override;
    void visit(std::shared_ptr<IRFunction> function) override;
    void visit(std::shared_ptr<IRExpression> expression) override;
    void visit(std::shared_ptr<IRStatement> statement) override;
    void visitDimensionalProperties(std::shared_ptr<IRNode> node) override;
    void visitQuantumEnhancements(std::shared_ptr<IRNode> node) override;
    
    // Analysis results
    std::map<std::string, double> getPerformanceMetrics() const;
    std::map<std::string, int> getComplexityAnalysis() const;
    
private:
    std::map<std::string, double> performanceMetrics_;
    std::map<std::string, int> complexityAnalysis_;
    
    void analyzePerformance(std::shared_ptr<IRNode> node);
    void analyzeComplexity(std::shared_ptr<IRNode> node);
};

} // namespace language
} // namespace hsml