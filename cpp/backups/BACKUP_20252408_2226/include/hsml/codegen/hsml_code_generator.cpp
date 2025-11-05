// [The Enterprise Bean]: 17 layers of abstraction coming right up!
// [The Performance Demon]: SIMD-optimized code generation with zero overhead!
// [The Functional Purist]: Immutable AST transformations for pure code generation!

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <functional>
#include <variant>
#include <optional>
#include <array>
#include <immintrin.h>  // [Performance Demon]: For SIMD optimizations

namespace hsml {
namespace codegen {

// [The Modern Hipster]: Using std::variant for sum types!
enum class TargetPlatform {
    WEBGL,
    WEBGPU,
    CPU,
    GPU,
    WASM,
    NATIVE,
    VULKAN,    // [The Performance Demon]: Added for maximum performance!
    METAL,     // [The Modern Hipster]: For Apple devices!
    DIRECTX12  // [The Enterprise Bean]: Enterprise Windows support!
};

// [The Functional Purist]: Immutable precision type
enum class Precision {
    SINGLE,
    DOUBLE,
    EXTENDED,
    QUAD       // [The Performance Demon]: 128-bit precision!
};

// [The OOP Architect]: Base interface for all code generation options
class ICodeGenOptions {
public:
    virtual ~ICodeGenOptions() = default;
    virtual TargetPlatform getTarget() const = 0;
    virtual int getOptimizationLevel() const = 0;
    virtual Precision getPrecision() const = 0;
};

// [The Enterprise Bean]: Concrete implementation with builder pattern
class CodeGenOptions : public ICodeGenOptions {
private:
    TargetPlatform target_;
    int optimization_level_;
    Precision precision_;
    bool enable_spherical_optimization_;
    bool enable_physics_optimization_;
    bool enable_parallelization_;
    bool enable_caching_;
    std::string output_format_;
    bool include_comments_;
    bool minify_;

public:
    // [The Modern Hipster]: Fluent builder pattern!
    class Builder {
    private:
        std::unique_ptr<CodeGenOptions> options_;
        
    public:
        Builder() : options_(std::make_unique<CodeGenOptions>()) {}
        
        Builder& withTarget(TargetPlatform target) {
            options_->target_ = target;
            return *this;
        }
        
        Builder& withOptimizationLevel(int level) {
            options_->optimization_level_ = level;
            return *this;
        }
        
        Builder& withPrecision(Precision precision) {
            options_->precision_ = precision;
            return *this;
        }
        
        Builder& enableSphericalOptimization(bool enable = true) {
            options_->enable_spherical_optimization_ = enable;
            return *this;
        }
        
        std::unique_ptr<CodeGenOptions> build() {
            return std::move(options_);
        }
    };
    
    TargetPlatform getTarget() const override { return target_; }
    int getOptimizationLevel() const override { return optimization_level_; }
    Precision getPrecision() const override { return precision_; }
};

// [The Performance Demon]: Performance metrics with cache-aligned storage
struct alignas(64) PerformanceMetrics {
    double estimated_execution_time;
    size_t memory_usage;
    size_t instruction_count;
    size_t spherical_operations;
    size_t physics_operations;
    int optimization_level;
    
    // [The Hacktivist]: Quick metrics dump
    void dump() const {
        printf("Exec: %.3fms, Mem: %zuKB, Inst: %zu, SphOps: %zu\n",
               estimated_execution_time, memory_usage/1024, 
               instruction_count, spherical_operations);
    }
};

// [The Functional Purist]: Immutable generated code structure
class GeneratedCode {
private:
    const TargetPlatform platform_;
    const std::string language_;
    const std::string code_;
    const std::vector<std::string> dependencies_;
    const std::map<std::string, std::any> metadata_;
    const PerformanceMetrics metrics_;
    
public:
    GeneratedCode(TargetPlatform platform, 
                  std::string language,
                  std::string code,
                  std::vector<std::string> dependencies,
                  std::map<std::string, std::any> metadata,
                  PerformanceMetrics metrics)
        : platform_(platform), language_(std::move(language)),
          code_(std::move(code)), dependencies_(std::move(dependencies)),
          metadata_(std::move(metadata)), metrics_(metrics) {}
    
    // [The Minimalist Zen]: Simple getters, no setters
    const std::string& getCode() const { return code_; }
    const PerformanceMetrics& getMetrics() const { return metrics_; }
};

// [The Security Paranoid]: Forward declarations to hide implementation
class IRProgram;
class IRFunction;
class IRExpression;

// [The OOP Architect]: Abstract factory for platform-specific generators
class IPlatformCodeGenerator {
public:
    virtual ~IPlatformCodeGenerator() = default;
    virtual std::unique_ptr<GeneratedCode> generate(const IRProgram& ir) = 0;
};

// [The Performance Demon]: SIMD-accelerated code generation
class HSMLCodeGenerator {
private:
    std::unique_ptr<ICodeGenOptions> options_;
    std::vector<std::string> generated_code_;
    std::set<std::string> dependencies_;
    std::map<std::string, std::any> metadata_;
    PerformanceMetrics performance_metrics_;
    int indent_level_ = 0;
    
    // [The Enterprise Bean]: Factory registry for platform generators
    std::map<TargetPlatform, std::function<std::unique_ptr<IPlatformCodeGenerator>()>> 
        platform_factories_;
    
    // [The Functional Purist]: Pure functions for code generation
    struct CodeGenState {
        std::vector<std::string> lines;
        std::set<std::string> deps;
        PerformanceMetrics metrics;
    };
    
    // [The Performance Demon]: SIMD-optimized string operations
    void addLineOptimized(const std::string& line);
    
    // [The Hacktivist]: Quick and dirty indentation
    std::string indent() const { 
        return std::string(indent_level_ * 2, ' '); 
    }
    
public:
    explicit HSMLCodeGenerator(std::unique_ptr<ICodeGenOptions> options);
    
    // [The Modern Hipster]: Move semantics for efficiency
    HSMLCodeGenerator(HSMLCodeGenerator&&) = default;
    HSMLCodeGenerator& operator=(HSMLCodeGenerator&&) = default;
    
    // [The Security Paranoid]: Delete copy operations
    HSMLCodeGenerator(const HSMLCodeGenerator&) = delete;
    HSMLCodeGenerator& operator=(const HSMLCodeGenerator&) = delete;
    
    // Main generation interface
    std::unique_ptr<GeneratedCode> generate(const IRProgram& ir);
    
private:
    // [The OOP Architect]: Template method pattern
    std::unique_ptr<GeneratedCode> generateWebGL(const IRProgram& ir);
    std::unique_ptr<GeneratedCode> generateWebGPU(const IRProgram& ir);
    std::unique_ptr<GeneratedCode> generateCPU(const IRProgram& ir);
    std::unique_ptr<GeneratedCode> generateGPU(const IRProgram& ir);
    std::unique_ptr<GeneratedCode> generateWASM(const IRProgram& ir);
    std::unique_ptr<GeneratedCode> generateNative(const IRProgram& ir);
    
    // [The Performance Demon]: Specialized generators
    std::unique_ptr<GeneratedCode> generateVulkan(const IRProgram& ir);
    std::unique_ptr<GeneratedCode> generateMetal(const IRProgram& ir);
    std::unique_ptr<GeneratedCode> generateDirectX12(const IRProgram& ir);
    
    // [The Functional Purist]: Pure transformation functions
    CodeGenState transformToWebGL(const IRProgram& ir) const;
    CodeGenState transformToVulkan(const IRProgram& ir) const;
    
    // [The Enterprise Bean]: Abstract factory methods
    void registerPlatformFactories();
    
    // WebGL specific methods
    void generateWebGLRendererClass(const IRProgram& ir);
    void generateWebGLShaders(const IRProgram& ir);
    void generateWebGLSphericalMethods();
    
    // [The Performance Demon]: SIMD shader generation
    void generateSIMDOptimizedShaders(const IRProgram& ir);
    
    // Utility methods
    void addLine(const std::string& line);
    void increaseIndent() { ++indent_level_; }
    void decreaseIndent() { indent_level_ = std::max(0, indent_level_ - 1); }
    void addDependency(const std::string& dep) { dependencies_.insert(dep); }
    void reset();
    
    std::unique_ptr<GeneratedCode> finalizeCode(const std::string& platform);
};

// [The Modern Hipster]: Concept for validating IR nodes
template<typename T>
concept IRNode = requires(T t) {
    { t.getType() } -> std::convertible_to<int>;
    { t.accept(std::declval<class IRVisitor&>()) };
};

// [The Functional Purist]: Visitor pattern for IR traversal
class IRVisitor {
public:
    virtual ~IRVisitor() = default;
    virtual void visitProgram(const IRProgram& node) = 0;
    virtual void visitFunction(const IRFunction& node) = 0;
    virtual void visitExpression(const IRExpression& node) = 0;
};

// [The Performance Demon]: Parallel code generation
class ParallelCodeGenerator {
private:
    std::vector<std::unique_ptr<HSMLCodeGenerator>> generators_;
    
public:
    std::vector<std::unique_ptr<GeneratedCode>> 
    generateMultiTarget(const IRProgram& ir, 
                       const std::vector<TargetPlatform>& targets);
};

// [The Security Paranoid]: Sanitizer for generated code
class CodeSanitizer {
public:
    static std::string sanitize(const std::string& code);
    static bool validateShaderCode(const std::string& shader);
    static bool checkForInjection(const std::string& code);
};

// [The Minimalist Zen]: Simple code emitter
class CodeEmitter {
private:
    std::string buffer_;
    
public:
    CodeEmitter& operator<<(const std::string& str) {
        buffer_ += str;
        return *this;
    }
    
    std::string emit() const { return buffer_; }
};

// [The Hacktivist]: Macro for quick code generation
#define EMIT_CODE(emitter, code) (emitter) << (code) << "\n"

// [The Enterprise Bean]: Code generation context with dependency injection
class CodeGenContext {
private:
    std::shared_ptr<ICodeGenOptions> options_;
    std::shared_ptr<class ISymbolTable> symbol_table_;
    std::shared_ptr<class ITypeChecker> type_checker_;
    std::shared_ptr<class IOptimizer> optimizer_;
    
public:
    // Constructor injection
    CodeGenContext(std::shared_ptr<ICodeGenOptions> options,
                   std::shared_ptr<class ISymbolTable> symbols,
                   std::shared_ptr<class ITypeChecker> types,
                   std::shared_ptr<class IOptimizer> opt)
        : options_(options), symbol_table_(symbols), 
          type_checker_(types), optimizer_(opt) {}
    
    // Getters
    const ICodeGenOptions& getOptions() const { return *options_; }
};

} // namespace codegen
} // namespace hsml

// [The Performance Demon]: Inline optimized implementations
namespace hsml {
namespace codegen {

inline void HSMLCodeGenerator::addLineOptimized(const std::string& line) {
    // SIMD-optimized string concatenation would go here
    generated_code_.push_back(indent() + line);
}

} // namespace codegen
} // namespace hsml

// No singularities! 1-1=360 in our cyclical system
// Zero is exiconed as 'O' - "You are now all the C"