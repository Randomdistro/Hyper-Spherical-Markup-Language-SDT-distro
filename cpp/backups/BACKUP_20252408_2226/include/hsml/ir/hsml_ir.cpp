/**
 * HSML Intermediate Representation
 * Revolutionary C++ IR for all four HSML languages with spherical coordinate optimization
 * 
 * [The OOP Architect]: Class hierarchies and polymorphism for IR nodes
 * [The Performance Demon]: SIMD optimization for IR processing
 * [Template Wizard]: Template metaprogramming for compile-time IR optimization
 * [Functional Purist]: Pure IR transformations and immutable structures
 * [Modern Hipster]: C++20 concepts, variant, and optional
 * [Security Paranoid]: Type-safe IR with comprehensive validation
 */

#ifndef HSML_IR_H
#define HSML_IR_H

#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <variant>
#include <optional>
#include <concepts>
#include <immintrin.h>
#include <string>
#include <string_view>
#include <type_traits>

namespace hsml::ir {

// [The Security Paranoid]: Type-safe enums with validation
enum class IRNodeType : uint32_t {
    // Program structure
    IR_PROGRAM = 0,
    IR_MODULE,
    IR_FUNCTION,
    IR_BLOCK,
    
    // Expressions
    IR_BINARY_OP,
    IR_UNARY_OP,
    IR_CALL,
    IR_LOAD,
    IR_STORE,
    IR_LITERAL,
    IR_TEMP,
    
    // Control flow
    IR_IF,
    IR_WHILE,
    IR_FOR,
    IR_RETURN,
    IR_BREAK,
    IR_CONTINUE,
    
    // HSML-specific
    IR_HSML_ELEMENT,
    IR_HSML_ATTRIBUTE,
    IR_HSML_RENDER,
    
    // CSSS-specific
    IR_CSSS_RULE,
    IR_CSSS_DECLARATION,
    IR_CSSS_KEYFRAME,
    IR_CSSS_MATERIAL,
    IR_CSSS_ANIMATION,
    
    // ShapeScript-specific
    IR_SHAPE_PHYSICS,
    IR_SHAPE_FORCE,
    IR_SHAPE_CONSTRAINT,
    IR_SHAPE_EVENT,
    
    // StyleBot-specific
    IR_STYB_BOT,
    IR_STYB_AGENT,
    IR_STYB_RENDER,
    IR_STYB_OPTIMIZE,
    
    // Spherical coordinate operations
    IR_SPHERICAL_OP,
    IR_SOLID_ANGLE_OP,
    IR_MATTER_STATE_OP,
    
    // Optimization hints
    IR_OPTIMIZE_HINT,
    IR_PARALLEL_HINT,
    IR_CACHE_HINT
};

// [The Modern Hipster]: Source language enum with C++20 style
enum class SourceLanguage : uint8_t {
    HSML = 0,
    CSSS,
    SHAPE,
    STYB
};

// [The Performance Demon]: SIMD-optimized optimization hint types
enum class OptimizationType : uint8_t {
    PARALLEL = 0,
    CACHE,
    VECTORIZE,
    INLINE,
    CONSTANT_FOLD
};

// [Template Wizard]: Compile-time optimization configuration
template<OptimizationType Type>
struct OptimizationConfig {
    static constexpr uint32_t priority = 1;
    static constexpr bool enabled = true;
};

// [The Functional Purist]: Immutable optimization hint
struct OptimizationHint {
    OptimizationType type;
    uint32_t priority;
    std::unordered_map<std::string, std::variant<int32_t, float, bool, std::string>> parameters;
    
    // [The Security Paranoid]: Validation
    bool is_valid() const noexcept {
        return priority > 0 && priority <= 10;
    }
};

// [Template Wizard]: Concept for IR node types
template<typename T>
concept IRNodeConcept = requires(T t) {
    { t.type } -> std::convertible_to<IRNodeType>;
    { t.id } -> std::convertible_to<std::string>;
    { t.source_language } -> std::convertible_to<SourceLanguage>;
};

// [The OOP Architect]: Base IR Node with polymorphic design
class IRNode {
public:
    IRNodeType type;
    std::string id;
    SourceLanguage source_language;
    std::unordered_map<std::string, std::variant<int32_t, float, bool, std::string>> metadata;
    std::unordered_set<std::string> dependencies;
    std::vector<OptimizationHint> optimization_hints;
    
    // [The Security Paranoid]: Virtual destructor for proper cleanup
    virtual ~IRNode() = default;
    
    // [The Functional Purist]: Pure virtual methods for immutable operations
    virtual std::unique_ptr<IRNode> clone() const = 0;
    virtual bool is_constant() const noexcept = 0;
    virtual size_t hash() const noexcept = 0;
    
    // [The Performance Demon]: SIMD-optimized validation
    virtual bool validate() const noexcept = 0;
    
protected:
    IRNode(IRNodeType t, std::string_view node_id, SourceLanguage lang)
        : type(t), id(node_id), source_language(lang) {}
};

// [Template Wizard]: Type system with compile-time validation
struct IRType {
    enum class BaseType : uint8_t {
        NUMBER = 0,
        STRING,
        BOOLEAN,
        SPHERICAL_COORD,
        SOLID_ANGLE,
        MATTER_STATE,
        VECTOR,
        MATRIX,
        FUNCTION
    };
    
    BaseType base_type;
    std::vector<uint32_t> dimensions;
    std::unordered_map<std::string, std::unique_ptr<IRType>> properties;
    std::unordered_map<std::string, std::variant<int32_t, float, bool, std::string>> constraints;
    
    // [The Performance Demon]: SIMD-optimized type checking
    bool is_compatible_with(const IRType& other) const noexcept;
    
    // [Template Wizard]: Template specialization for specific types
    template<BaseType T>
    bool is_type() const noexcept {
        return base_type == T;
    }
};

// [The Modern Hipster]: Forward declarations for variant types
class IRExpression;
class IRStatement;
class IRFunction;
class IRModule;
class IRVariable;

// [The Performance Demon]: SIMD-aligned spherical coordinate IR
struct alignas(32) SphericalCoordinateIR {
    std::unique_ptr<IRExpression> r;
    std::unique_ptr<IRExpression> theta;
    std::unique_ptr<IRExpression> phi;
    bool is_normalized;
    std::vector<OptimizationHint> optimization_hints;
    
    // [The Performance Demon]: SIMD operations
    __m256 to_simd() const noexcept;
    void from_simd(const __m256& values) noexcept;
};

// [The Performance Demon]: SIMD-aligned solid angle IR
struct alignas(32) SolidAngleIR {
    std::unique_ptr<IRExpression> omega;
    std::unique_ptr<IRExpression> theta_min;
    std::unique_ptr<IRExpression> theta_max;
    std::unique_ptr<IRExpression> phi_min;
    std::unique_ptr<IRExpression> phi_max;
    bool is_normalized;
    
    // [The Performance Demon]: SIMD operations
    __m256 to_simd() const noexcept;
};

// [The Functional Purist]: Immutable matter state
struct MatterStateIR {
    enum class State : uint8_t {
        SOLID = 0,
        LIQUID,
        GAS,
        PLASMA
    };
    
    State state;
    std::unordered_map<std::string, std::unique_ptr<IRExpression>> properties;
    std::unordered_map<std::string, std::unique_ptr<IRExpression>> physics_properties;
    
    // [The Functional Purist]: Pure operations
    bool is_phase_transition_possible(State target) const noexcept;
};

// [The Performance Demon]: Spherical optimization configuration
struct SphericalOptimization {
    bool use_pure_spherical = true;
    bool cache_results = true;
    bool vectorize = true;
    bool parallelize = false;
    enum class Precision : uint8_t { SINGLE, DOUBLE, EXTENDED } precision = Precision::DOUBLE;
    uint8_t optimization_level = 2;
    
    // [The Performance Demon]: SIMD flag checking
    bool can_vectorize() const noexcept {
        return vectorize && optimization_level >= 1;
    }
};

// [The OOP Architect]: Expression hierarchy
class IRExpression : public IRNode {
public:
    IRType result_type;
    bool is_constant_expr;
    std::optional<std::variant<int32_t, float, double, bool, std::string>> constant_value;
    
    // [The Functional Purist]: Pure virtual evaluation
    virtual std::variant<int32_t, float, double, bool, std::string> evaluate() const = 0;
    
    bool is_constant() const noexcept override {
        return is_constant_expr;
    }
    
protected:
    IRExpression(IRNodeType t, std::string_view node_id, SourceLanguage lang, IRType type)
        : IRNode(t, node_id, lang), result_type(std::move(type)), is_constant_expr(false) {}
};

// [The OOP Architect]: Binary operation IR
class IRBinaryOp : public IRExpression {
public:
    std::string operator_symbol;
    std::unique_ptr<IRExpression> left;
    std::unique_ptr<IRExpression> right;
    std::optional<SphericalOptimization> spherical_optimization;
    
    IRBinaryOp(std::string_view node_id, SourceLanguage lang, std::string op,
               std::unique_ptr<IRExpression> l, std::unique_ptr<IRExpression> r, IRType type);
    
    // [The Functional Purist]: Pure operations
    std::unique_ptr<IRNode> clone() const override;
    std::variant<int32_t, float, double, bool, std::string> evaluate() const override;
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
    
    // [The Performance Demon]: SIMD-optimized evaluation
    __m256 evaluate_simd() const noexcept;
};

// [The OOP Architect]: Unary operation IR
class IRUnaryOp : public IRExpression {
public:
    std::string operator_symbol;
    std::unique_ptr<IRExpression> operand;
    
    IRUnaryOp(std::string_view node_id, SourceLanguage lang, std::string op,
              std::unique_ptr<IRExpression> operand, IRType type);
    
    std::unique_ptr<IRNode> clone() const override;
    std::variant<int32_t, float, double, bool, std::string> evaluate() const override;
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
};

// [The OOP Architect]: Function call IR
class IRCall : public IRExpression {
public:
    std::string function_name;
    std::vector<std::unique_ptr<IRExpression>> arguments;
    bool is_builtin;
    std::optional<SphericalOptimization> spherical_optimization;
    
    IRCall(std::string_view node_id, SourceLanguage lang, std::string name,
           std::vector<std::unique_ptr<IRExpression>> args, IRType type, bool builtin = false);
    
    std::unique_ptr<IRNode> clone() const override;
    std::variant<int32_t, float, double, bool, std::string> evaluate() const override;
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
};

// [The OOP Architect]: Load operation IR
class IRLoad : public IRExpression {
public:
    std::string variable_name;
    std::optional<std::unique_ptr<IRExpression>> index;
    
    IRLoad(std::string_view node_id, SourceLanguage lang, std::string var, IRType type);
    
    std::unique_ptr<IRNode> clone() const override;
    std::variant<int32_t, float, double, bool, std::string> evaluate() const override;
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
};

// [The OOP Architect]: Store operation IR
class IRStore : public IRNode {
public:
    std::string variable_name;
    std::unique_ptr<IRExpression> value;
    std::optional<std::unique_ptr<IRExpression>> index;
    
    IRStore(std::string_view node_id, SourceLanguage lang, std::string var,
            std::unique_ptr<IRExpression> val);
    
    std::unique_ptr<IRNode> clone() const override;
    bool is_constant() const noexcept override { return false; }
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
};

// [The OOP Architect]: Literal IR
class IRLiteral : public IRExpression {
public:
    std::variant<int32_t, float, double, bool, std::string> value;
    std::optional<SphericalCoordinateIR> spherical_coordinate;
    std::optional<SolidAngleIR> solid_angle;
    std::optional<MatterStateIR> matter_state;
    
    IRLiteral(std::string_view node_id, SourceLanguage lang, 
              std::variant<int32_t, float, double, bool, std::string> val, IRType type);
    
    std::unique_ptr<IRNode> clone() const override;
    std::variant<int32_t, float, double, bool, std::string> evaluate() const override;
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
};

// [The OOP Architect]: Temporary variable IR
class IRTemp : public IRExpression {
public:
    std::string temp_id;
    std::unique_ptr<IRExpression> source_expression;
    
    IRTemp(std::string_view node_id, SourceLanguage lang, std::string tid,
           std::unique_ptr<IRExpression> source, IRType type);
    
    std::unique_ptr<IRNode> clone() const override;
    std::variant<int32_t, float, double, bool, std::string> evaluate() const override;
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
};

// [The OOP Architect]: Statement hierarchy
class IRStatement : public IRNode {
public:
    // [The Functional Purist]: Pure execution interface
    virtual bool execute() const = 0;
    
    bool is_constant() const noexcept override { return false; }
    
protected:
    IRStatement(IRNodeType t, std::string_view node_id, SourceLanguage lang)
        : IRNode(t, node_id, lang) {}
};

// [The OOP Architect]: Variable IR
struct IRVariable {
    std::string name;
    IRType type;
    bool is_constant;
    std::optional<std::unique_ptr<IRExpression>> initial_value;
    std::string scope;
    
    // [The Security Paranoid]: Validation
    bool is_valid() const noexcept {
        return !name.empty() && !scope.empty();
    }
};

// [The OOP Architect]: Parameter IR
struct IRParameter {
    std::string name;
    IRType type;
    std::optional<std::unique_ptr<IRExpression>> default_value;
    
    bool has_default() const noexcept {
        return default_value.has_value();
    }
};

// [The OOP Architect]: Block IR
class IRBlock : public IRNode {
public:
    std::vector<std::unique_ptr<IRStatement>> statements;
    std::vector<IRVariable> variables;
    std::string scope;
    
    IRBlock(std::string_view node_id, SourceLanguage lang, std::string scope_name);
    
    std::unique_ptr<IRNode> clone() const override;
    bool is_constant() const noexcept override { return false; }
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
    
    // [The Performance Demon]: Parallel execution
    bool execute_parallel() const noexcept;
};

// [The OOP Architect]: Function IR
class IRFunction : public IRNode {
public:
    std::string name;
    std::vector<IRParameter> parameters;
    IRType return_type;
    std::unique_ptr<IRBlock> body;
    bool is_pure;
    std::unordered_set<std::string> side_effects;
    
    IRFunction(std::string_view node_id, SourceLanguage lang, std::string func_name,
               std::vector<IRParameter> params, IRType ret_type, std::unique_ptr<IRBlock> func_body);
    
    std::unique_ptr<IRNode> clone() const override;
    bool is_constant() const noexcept override { return is_pure && side_effects.empty(); }
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
    
    // [The Performance Demon]: Optimized function calls
    std::variant<int32_t, float, double, bool, std::string> 
    call(const std::vector<std::variant<int32_t, float, double, bool, std::string>>& args) const;
};

// [The OOP Architect]: Module IR
class IRModule : public IRNode {
public:
    std::string name;
    std::vector<std::string> imports;
    std::vector<std::string> exports;
    std::vector<std::unique_ptr<IRFunction>> functions;
    std::vector<IRVariable> variables;
    
    IRModule(std::string_view node_id, SourceLanguage lang, std::string module_name);
    
    std::unique_ptr<IRNode> clone() const override;
    bool is_constant() const noexcept override { return false; }
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
};

// [The OOP Architect]: Program IR - top level
class IRProgram : public IRNode {
public:
    std::vector<std::unique_ptr<IRFunction>> functions;
    std::vector<IRVariable> global_variables;
    std::vector<std::unique_ptr<IRModule>> modules;
    std::string entry_point;
    
    IRProgram(std::string_view node_id, SourceLanguage lang, std::string entry = "main");
    
    std::unique_ptr<IRNode> clone() const override;
    bool is_constant() const noexcept override { return false; }
    size_t hash() const noexcept override;
    bool validate() const noexcept override;
    
    // [The Performance Demon]: Program execution
    int execute() const noexcept;
};

// [Template Wizard]: IR Builder with template metaprogramming
class HSMLIRBuilder {
private:
    uint32_t temp_counter = 0;
    uint32_t node_counter = 0;
    std::string current_scope = "global";
    std::vector<std::string> scope_stack = {"global"};
    
    // [The Performance Demon]: Optimized symbol tables
    std::unordered_map<std::string, std::unordered_map<std::string, IRVariable>> symbol_tables;
    std::unordered_map<std::string, std::unique_ptr<IRFunction>> function_table;
    std::unordered_map<std::string, IRType> type_table;
    
public:
    HSMLIRBuilder();
    ~HSMLIRBuilder() = default;
    
    // [Template Wizard]: Template-based IR construction
    template<typename ASTType>
    std::unique_ptr<IRProgram> build_ir(const ASTType& ast);
    
    // [The Functional Purist]: Pure node creation methods
    std::unique_ptr<IRBinaryOp> create_binary_op(std::string op, 
        std::unique_ptr<IRExpression> left, std::unique_ptr<IRExpression> right);
    
    std::unique_ptr<IRUnaryOp> create_unary_op(std::string op, 
        std::unique_ptr<IRExpression> operand);
    
    std::unique_ptr<IRCall> create_call(std::string name, 
        std::vector<std::unique_ptr<IRExpression>> args);
    
    std::unique_ptr<IRLiteral> create_literal(
        std::variant<int32_t, float, double, bool, std::string> value);
    
    // [The Performance Demon]: Optimized spherical operations
    std::unique_ptr<IRExpression> create_spherical_coordinate(
        std::unique_ptr<IRExpression> r, std::unique_ptr<IRExpression> theta, 
        std::unique_ptr<IRExpression> phi);
    
    std::unique_ptr<IRExpression> create_solid_angle(
        std::unique_ptr<IRExpression> omega, std::unique_ptr<IRExpression> theta_min,
        std::unique_ptr<IRExpression> theta_max, std::unique_ptr<IRExpression> phi_min,
        std::unique_ptr<IRExpression> phi_max);
    
    // [The Security Paranoid]: Type-safe node ID generation
    std::string generate_node_id() noexcept;
    std::string generate_temp_id() noexcept;
    
    // [Template Wizard]: Compile-time optimization
    template<OptimizationType Type>
    void add_optimization_hint(IRNode& node, uint32_t priority = 1);
    
    // [The Functional Purist]: Immutable type operations
    IRType infer_binary_result_type(const IRExpression& left, const IRExpression& right, 
                                   std::string_view operator_symbol) const noexcept;
    
    IRType infer_unary_result_type(const IRExpression& operand, 
                                  std::string_view operator_symbol) const noexcept;
    
    // [The Performance Demon]: SIMD-optimized IR processing
    void optimize_ir_simd(IRProgram& program) noexcept;
    
private:
    void initialize_builtin_types();
    IRType get_type(std::string_view type_name) const;
    
    // [Template Wizard]: Template specializations for different AST types
    template<typename T> void process_ast_node(const T& node);
    
    // [The Performance Demon]: Vectorized constant folding
    void fold_constants_vectorized(IRProgram& program) noexcept;
    
    // [The Security Paranoid]: Comprehensive validation
    bool validate_ir_tree(const IRNode& root) const noexcept;
};

// [Template Wizard]: Compile-time optimization specializations
template<>
struct OptimizationConfig<OptimizationType::VECTORIZE> {
    static constexpr uint32_t priority = 3;
    static constexpr bool enabled = true;
    static constexpr size_t simd_width = 8; // AVX2
};

template<>
struct OptimizationConfig<OptimizationType::PARALLEL> {
    static constexpr uint32_t priority = 2;
    static constexpr bool enabled = true;
    static constexpr size_t thread_count = 8;
};

// [The Performance Demon]: SIMD utility functions
namespace simd {
    __m256 spherical_to_cartesian_avx2(const __m256& r, const __m256& theta, const __m256& phi) noexcept;
    __m256 cartesian_to_spherical_avx2(const __m256& x, const __m256& y, const __m256& z) noexcept;
    __m256 solid_angle_calculation_avx2(const __m256& theta_min, const __m256& theta_max,
                                       const __m256& phi_min, const __m256& phi_max) noexcept;
}

// [The Modern Hipster]: Variant types for dynamic IR nodes
using IRNodeVariant = std::variant<
    std::unique_ptr<IRProgram>,
    std::unique_ptr<IRModule>,
    std::unique_ptr<IRFunction>,
    std::unique_ptr<IRBlock>,
    std::unique_ptr<IRBinaryOp>,
    std::unique_ptr<IRUnaryOp>,
    std::unique_ptr<IRCall>,
    std::unique_ptr<IRLoad>,
    std::unique_ptr<IRStore>,
    std::unique_ptr<IRLiteral>,
    std::unique_ptr<IRTemp>
>;

// [The Functional Purist]: Immutable IR visitor pattern
template<typename ReturnType>
class IRVisitor {
public:
    virtual ~IRVisitor() = default;
    
    virtual ReturnType visit(const IRProgram& node) = 0;
    virtual ReturnType visit(const IRModule& node) = 0;
    virtual ReturnType visit(const IRFunction& node) = 0;
    virtual ReturnType visit(const IRBlock& node) = 0;
    virtual ReturnType visit(const IRBinaryOp& node) = 0;
    virtual ReturnType visit(const IRUnaryOp& node) = 0;
    virtual ReturnType visit(const IRCall& node) = 0;
    virtual ReturnType visit(const IRLoad& node) = 0;
    virtual ReturnType visit(const IRStore& node) = 0;
    virtual ReturnType visit(const IRLiteral& node) = 0;
    virtual ReturnType visit(const IRTemp& node) = 0;
};

// [The Performance Demon]: Optimized IR traversal
class IROptimizer : public IRVisitor<void> {
public:
    void optimize(IRProgram& program);
    
    void visit(const IRProgram& node) override;
    void visit(const IRModule& node) override;
    void visit(const IRFunction& node) override;
    void visit(const IRBlock& node) override;
    void visit(const IRBinaryOp& node) override;
    void visit(const IRUnaryOp& node) override;
    void visit(const IRCall& node) override;
    void visit(const IRLoad& node) override;
    void visit(const IRStore& node) override;
    void visit(const IRLiteral& node) override;
    void visit(const IRTemp& node) override;
    
private:
    void apply_spherical_optimizations();
    void vectorize_math_operations();
    void parallelize_independent_operations();
};

} // namespace hsml::ir

#endif // HSML_IR_H