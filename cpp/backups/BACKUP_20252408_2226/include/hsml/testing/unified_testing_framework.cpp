/**
 * HSML Unified Testing Framework - Comprehensive Multi-Component Validation
 * Integrates unit tests, benchmarks, integration tests, and validation suites
 * Provides coherent testing infrastructure across all HSML components
 */

#pragma once

#include <memory>
#include <string>
#include <vector>
#include <functional>
#include <chrono>
#include <unordered_map>
#include <atomic>
#include <mutex>
#include <shared_mutex>
#include <fstream>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <concepts>

// Include component-specific testing headers
#include "hsml_testing_framework.h"
#include "comprehensive_spatial_validator.h"
#include "mpd_spatial_testing_framework.h"
#include "memory_leak_detector.h"

// Include components being tested
#include "../core/spatial_indexer_coordinator.h"
#include "../gui/unified_gui_coordinator.h"
#include "../rendering/unified_rendering_coordinator.h"
#include "../core/spherical_physics_engine.h"

namespace hsml::testing {

// Test result status
enum class TestStatus : uint8_t {
    PASSED,
    FAILED,
    SKIPPED,
    TIMEOUT,
    ERROR
};

// Test categories for organization
enum class TestCategory : uint8_t {
    UNIT,
    INTEGRATION,
    PERFORMANCE,
    STRESS,
    REGRESSION,
    VALIDATION,
    BENCHMARK
};

// Test priority levels
enum class TestPriority : uint8_t {
    CRITICAL,   // Must pass for release
    HIGH,       // Should pass for release
    MEDIUM,     // Important but not blocking
    LOW         // Nice to have
};

// Performance benchmarking metrics
struct PerformanceBenchmark {
    std::string metric_name;
    double value = 0.0;
    std::string unit;
    double baseline_value = 0.0;
    double tolerance_percent = 10.0;
    bool meets_baseline = false;
};

// Test result structure
struct TestResult {
    std::string test_name;
    TestCategory category;
    TestPriority priority;
    TestStatus status = TestStatus::FAILED;
    std::chrono::milliseconds execution_time{0};
    std::string message;
    std::string details;
    std::vector<PerformanceBenchmark> benchmarks;
    size_t memory_usage_kb = 0;
    std::chrono::system_clock::time_point timestamp;
    
    // Test configuration
    std::unordered_map<std::string, std::string> configuration;
    
    TestResult(const std::string& name, TestCategory cat, TestPriority prio)
        : test_name(name), category(cat), priority(prio), 
          timestamp(std::chrono::system_clock::now()) {}
};

// Test suite statistics
struct TestSuiteStats {
    size_t total_tests = 0;
    size_t passed_tests = 0;
    size_t failed_tests = 0;
    size_t skipped_tests = 0;
    size_t timeout_tests = 0;
    size_t error_tests = 0;
    std::chrono::milliseconds total_execution_time{0};
    double pass_rate = 0.0;
    std::string fastest_test;
    std::string slowest_test; 
    std::chrono::milliseconds fastest_time{0};
    std::chrono::milliseconds slowest_time{0};
};

// Concept for testable components
template<typename T>
concept TestableComponent = requires(T component) {
    { component.initialize() } -> std::convertible_to<bool>;
    { component.shutdown() } -> std::same_as<void>;
    { component.is_healthy() } -> std::convertible_to<bool>;
};

// Base test case interface
class TestCase {
public:
    virtual ~TestCase() = default;
    
    virtual TestResult run() = 0;
    virtual std::string get_name() const = 0;
    virtual TestCategory get_category() const = 0;
    virtual TestPriority get_priority() const = 0;
    virtual std::chrono::milliseconds get_timeout() const { return std::chrono::milliseconds{5000}; }
    virtual void setup() {}
    virtual void teardown() {}
};

// Template-based test case for specific components
template<TestableComponent ComponentType>
class ComponentTestCase : public TestCase {
protected:
    std::unique_ptr<ComponentType> component_;
    std::string test_name_;
    TestCategory category_;
    TestPriority priority_;

public:
    ComponentTestCase(const std::string& name, TestCategory category, TestPriority priority)
        : test_name_(name), category_(category), priority_(priority) {}
    
    std::string get_name() const override { return test_name_; }
    TestCategory get_category() const override { return category_; }
    TestPriority get_priority() const override { return priority_; }
    
    void setup() override {
        component_ = std::make_unique<ComponentType>();
    }
    
    void teardown() override {
        if (component_) {
            component_->shutdown();
            component_.reset();
        }
    }
    
protected:
    ComponentType* get_component() { return component_.get(); }
};

// Spatial Indexer Test Suite
class SpatialIndexerTestSuite {
public:
    static std::vector<std::unique_ptr<TestCase>> create_test_cases() {
        std::vector<std::unique_ptr<TestCase>> tests;
        
        // Unit tests for each spatial indexer personality
        tests.push_back(std::make_unique<SpatialIndexerBasicFunctionalityTest>());
        tests.push_back(std::make_unique<SpatialIndexerPerformanceTest>());
        tests.push_back(std::make_unique<SpatialIndexerCoordinatorTest>());
        tests.push_back(std::make_unique<SpatialIndexerAdaptiveSelectionTest>());
        
        return tests;
    }
    
private:
    // Specific test implementations
    class SpatialIndexerBasicFunctionalityTest : public ComponentTestCase<hsml::core::SpatialIndexerCoordinator> {
    public:
        SpatialIndexerBasicFunctionalityTest() 
            : ComponentTestCase("SpatialIndexer Basic Functionality", TestCategory::UNIT, TestPriority::CRITICAL) {}
        
        TestResult run() override {
            TestResult result(get_name(), get_category(), get_priority());
            auto start_time = std::chrono::steady_clock::now();
            
            try {
                auto* indexer = get_component();
                
                // Test basic insertion and querying
                hsml::core::SphericalCoords test_coord{1.0, 0.5, 1.0};
                size_t id = indexer->insert(test_coord);
                
                if (id == 0) {
                    result.status = TestStatus::FAILED;
                    result.message = "Failed to insert coordinate";
                    return result;
                }
                
                auto query_results = indexer->query_sphere(test_coord, 0.1);
                if (query_results.empty()) {
                    result.status = TestStatus::FAILED;
                    result.message = "Query failed to find inserted coordinate";
                    return result;
                }
                
                // Test removal
                bool removed = indexer->remove(id);
                if (!removed) {
                    result.status = TestStatus::FAILED;
                    result.message = "Failed to remove coordinate";
                    return result;
                }
                
                result.status = TestStatus::PASSED;
                result.message = "All basic functionality tests passed";
                
            } catch (const std::exception& e) {
                result.status = TestStatus::ERROR;
                result.message = "Exception thrown: " + std::string(e.what());
            }
            
            auto end_time = std::chrono::steady_clock::now();
            result.execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time);
            
            return result;
        }
    };
    
    class SpatialIndexerPerformanceTest : public ComponentTestCase<hsml::core::SpatialIndexerCoordinator> {
    public:
        SpatialIndexerPerformanceTest()
            : ComponentTestCase("SpatialIndexer Performance", TestCategory::PERFORMANCE, TestPriority::HIGH) {}
        
        TestResult run() override {
            TestResult result(get_name(), get_category(), get_priority());
            auto start_time = std::chrono::steady_clock::now();
            
            try {
                auto* indexer = get_component();
                
                // Performance test with large dataset
                const size_t test_size = 10000;
                std::vector<size_t> ids;
                ids.reserve(test_size);
                
                // Insertion benchmark
                auto insert_start = std::chrono::steady_clock::now();
                for (size_t i = 0; i < test_size; ++i) {
                    hsml::core::SphericalCoords coord{
                        1.0 + (i % 100) * 0.01,
                        (i % 180) * M_PI / 180.0,
                        (i % 360) * M_PI / 180.0
                    };
                    ids.push_back(indexer->insert(coord));
                }
                auto insert_end = std::chrono::steady_clock::now();
                
                // Query benchmark
                auto query_start = std::chrono::steady_clock::now();
                for (size_t i = 0; i < 1000; ++i) {
                    hsml::core::SphericalCoords query_coord{1.5, M_PI/4, M_PI/2};
                    auto results = indexer->query_sphere(query_coord, 0.5);
                }
                auto query_end = std::chrono::steady_clock::now();
                
                // Calculate benchmarks
                auto insert_time = std::chrono::duration_cast<std::chrono::microseconds>(
                    insert_end - insert_start).count() / 1000.0;
                auto query_time = std::chrono::duration_cast<std::chrono::microseconds>(
                    query_end - query_start).count() / 1000.0;
                
                result.benchmarks.push_back({
                    "Insertion Time", insert_time, "ms", 100.0, 20.0,
                    insert_time <= 120.0 // 20% tolerance
                });
                
                result.benchmarks.push_back({
                    "Query Time", query_time, "ms", 50.0, 30.0,
                    query_time <= 65.0 // 30% tolerance
                });
                
                // Check if all benchmarks meet baseline
                bool all_benchmarks_pass = std::all_of(result.benchmarks.begin(), result.benchmarks.end(),
                    [](const PerformanceBenchmark& bench) { return bench.meets_baseline; });
                
                result.status = all_benchmarks_pass ? TestStatus::PASSED : TestStatus::FAILED;
                result.message = all_benchmarks_pass ? 
                    "Performance benchmarks met" : "Some benchmarks failed to meet baseline";
                
            } catch (const std::exception& e) {
                result.status = TestStatus::ERROR;
                result.message = "Exception thrown: " + std::string(e.what());
            }
            
            auto end_time = std::chrono::steady_clock::now();
            result.execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time);
            
            return result;
        }
        
        std::chrono::milliseconds get_timeout() const override { 
            return std::chrono::milliseconds{30000}; // Longer timeout for performance tests
        }
    };
    
    class SpatialIndexerCoordinatorTest : public ComponentTestCase<hsml::core::SpatialIndexerCoordinator> {
    public:
        SpatialIndexerCoordinatorTest()
            : ComponentTestCase("SpatialIndexer Coordinator", TestCategory::INTEGRATION, TestPriority::HIGH) {}
        
        TestResult run() override {
            TestResult result(get_name(), get_category(), get_priority());
            auto start_time = std::chrono::steady_clock::now();
            
            try {
                auto* coordinator = get_component();
                
                // Test personality switching
                coordinator->set_personality(hsml::core::IndexerPersonalityType::MINIMAL);
                if (coordinator->get_current_personality() != hsml::core::IndexerPersonalityType::MINIMAL) {
                    result.status = TestStatus::FAILED;
                    result.message = "Failed to switch to MINIMAL personality";
                    return result;
                }
                
                coordinator->set_personality(hsml::core::IndexerPersonalityType::SIMD);
                if (coordinator->get_current_personality() != hsml::core::IndexerPersonalityType::SIMD) {
                    result.status = TestStatus::FAILED;
                    result.message = "Failed to switch to SIMD personality";
                    return result;
                }
                
                // Test adaptive mode
                coordinator->set_adaptive_mode(true);
                if (!coordinator->is_adaptive_mode_enabled()) {
                    result.status = TestStatus::FAILED;
                    result.message = "Failed to enable adaptive mode";
                    return result;
                }
                
                result.status = TestStatus::PASSED;
                result.message = "Coordinator functionality verified";
                
            } catch (const std::exception& e) {
                result.status = TestStatus::ERROR;
                result.message = "Exception thrown: " + std::string(e.what());
            }
            
            auto end_time = std::chrono::steady_clock::now();
            result.execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time);
            
            return result;
        }
    };
    
    class SpatialIndexerAdaptiveSelectionTest : public ComponentTestCase<hsml::core::SpatialIndexerCoordinator> {
    public:
        SpatialIndexerAdaptiveSelectionTest()
            : ComponentTestCase("SpatialIndexer Adaptive Selection", TestCategory::VALIDATION, TestPriority::MEDIUM) {}
        
        TestResult run() override {
            TestResult result(get_name(), get_category(), get_priority());
            auto start_time = std::chrono::steady_clock::now();
            
            try {
                auto* coordinator = get_component();
                
                // Test workload-based selection
                hsml::core::WorkloadCharacteristics low_workload{
                    .dataset_size = 100,
                    .query_frequency = 1.0,
                    .insertion_frequency = 0.1,
                    .has_simd_capability = false
                };
                
                coordinator->update_workload_characteristics(low_workload);
                coordinator->set_personality(hsml::core::IndexerPersonalityType::ADAPTIVE);
                
                // Should select MINIMAL for low workload
                if (coordinator->get_current_personality() != hsml::core::IndexerPersonalityType::MINIMAL) {
                    result.status = TestStatus::FAILED;
                    result.message = "Adaptive selection failed for low workload";
                    return result;
                }
                
                hsml::core::WorkloadCharacteristics high_workload{
                    .dataset_size = 50000,
                    .query_frequency = 100.0,
                    .insertion_frequency = 10.0,
                    .has_simd_capability = true
                };
                
                coordinator->update_workload_characteristics(high_workload);
                coordinator->set_personality(hsml::core::IndexerPersonalityType::ADAPTIVE);
                
                // Should select SIMD or PERFORMANCE_DEMON for high workload
                auto selected = coordinator->get_current_personality();
                if (selected != hsml::core::IndexerPersonalityType::SIMD &&
                    selected != hsml::core::IndexerPersonalityType::PERFORMANCE_DEMON) {
                    result.status = TestStatus::FAILED;
                    result.message = "Adaptive selection failed for high workload";
                    return result;
                }
                
                result.status = TestStatus::PASSED;
                result.message = "Adaptive selection working correctly";
                
            } catch (const std::exception& e) {
                result.status = TestStatus::ERROR;
                result.message = "Exception thrown: " + std::string(e.what());
            }
            
            auto end_time = std::chrono::steady_clock::now();
            result.execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time);
            
            return result;
        }
    };
};

// Main unified testing framework
class UnifiedTestingFramework {
private:
    std::vector<std::unique_ptr<TestCase>> test_cases_;
    std::vector<TestResult> test_results_;
    std::atomic<size_t> tests_completed_{0};
    std::atomic<bool> stop_requested_{false};
    mutable std::shared_mutex framework_mutex_;
    
    // Configuration
    bool parallel_execution_enabled_ = true;
    size_t max_parallel_tests_ = std::thread::hardware_concurrency();
    bool continue_on_failure_ = true;
    bool generate_reports_ = true;
    std::string output_directory_ = "./test_results/";
    
    // Memory leak detection
    std::unique_ptr<MemoryLeakDetector> memory_detector_;

public:
    UnifiedTestingFramework() : memory_detector_(std::make_unique<MemoryLeakDetector>()) {
        register_all_test_suites();
    }
    
    // Register test suites from all components
    void register_all_test_suites() {
        // Spatial Indexer tests
        auto spatial_tests = SpatialIndexerTestSuite::create_test_cases();
        for (auto& test : spatial_tests) {
            test_cases_.push_back(std::move(test));
        }
        
        // TODO: Add other component test suites here
        // - GUI Framework tests
        // - Rendering Pipeline tests  
        // - Physics Engine tests
        // - SDT Integration tests
        // - P0rt3r Browser tests
    }
    
    // Run all tests
    TestSuiteStats run_all_tests() {
        std::lock_guard<std::shared_mutex> lock(framework_mutex_);
        
        memory_detector_->start_monitoring();
        
        test_results_.clear();
        test_results_.reserve(test_cases_.size());
        tests_completed_.store(0);
        stop_requested_.store(false);
        
        auto start_time = std::chrono::steady_clock::now();
        
        if (parallel_execution_enabled_) {
            run_tests_parallel();
        } else {
            run_tests_sequential();
        }
        
        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - start_time);
        
        memory_detector_->stop_monitoring();
        
        // Generate test statistics
        TestSuiteStats stats = calculate_statistics(total_time);
        
        if (generate_reports_) {
            generate_test_reports(stats);
        }
        
        return stats;
    }
    
    // Run specific test category
    TestSuiteStats run_category_tests(TestCategory category) {
        std::lock_guard<std::shared_mutex> lock(framework_mutex_);
        
        test_results_.clear();
        tests_completed_.store(0);
        
        auto start_time = std::chrono::steady_clock::now();
        
        for (auto& test_case : test_cases_) {
            if (test_case->get_category() == category) {
                if (stop_requested_.load()) break;
                
                TestResult result = execute_single_test(*test_case);
                test_results_.push_back(std::move(result));
                tests_completed_.fetch_add(1);
            }
        }
        
        auto end_time = std::chrono::steady_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            end_time - start_time);
        
        return calculate_statistics(total_time);
    }
    
    // Configuration methods
    void set_parallel_execution(bool enabled) { parallel_execution_enabled_ = enabled; }
    void set_max_parallel_tests(size_t max_tests) { max_parallel_tests_ = max_tests; }
    void set_continue_on_failure(bool continue_on_fail) { continue_on_failure_ = continue_on_fail; }
    void set_generate_reports(bool generate) { generate_reports_ = generate; }
    void set_output_directory(const std::string& dir) { output_directory_ = dir; }
    
    // Get test results
    const std::vector<TestResult>& get_test_results() const { return test_results_; }
    
    // Stop test execution
    void stop_tests() { stop_requested_.store(true); }

private:
    void run_tests_sequential() {
        for (auto& test_case : test_cases_) {
            if (stop_requested_.load()) break;
            
            TestResult result = execute_single_test(*test_case);
            test_results_.push_back(std::move(result));
            tests_completed_.fetch_add(1);
            
            if (!continue_on_failure_ && result.status == TestStatus::FAILED) {
                break;
            }
        }
    }
    
    void run_tests_parallel() {
        // Implement parallel test execution
        // For now, fall back to sequential
        run_tests_sequential();
    }
    
    TestResult execute_single_test(TestCase& test_case) {
        TestResult result(test_case.get_name(), test_case.get_category(), test_case.get_priority());
        
        try {
            // Setup phase
            test_case.setup();
            
            // Execute test with timeout
            auto timeout = test_case.get_timeout();
            auto start_time = std::chrono::steady_clock::now();
            
            result = test_case.run();
            
            auto end_time = std::chrono::steady_clock::now();
            auto execution_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time);
            
            if (execution_time > timeout) {
                result.status = TestStatus::TIMEOUT;
                result.message = "Test exceeded timeout of " + std::to_string(timeout.count()) + "ms";
            }
            
            result.execution_time = execution_time;
            result.memory_usage_kb = memory_detector_->get_current_usage_kb();
            
            // Teardown phase
            test_case.teardown();
            
        } catch (const std::exception& e) {
            result.status = TestStatus::ERROR;
            result.message = "Exception during test execution: " + std::string(e.what());
            
            try {
                test_case.teardown();
            } catch (...) {
                result.details += "\nAdditional exception during teardown";
            }
        }
        
        return result;
    }
    
    TestSuiteStats calculate_statistics(std::chrono::milliseconds total_time) {
        TestSuiteStats stats;
        stats.total_tests = test_results_.size();
        stats.total_execution_time = total_time;
        
        std::chrono::milliseconds fastest_time{std::chrono::milliseconds::max()};
        std::chrono::milliseconds slowest_time{0};
        
        for (const auto& result : test_results_) {
            switch (result.status) {
                case TestStatus::PASSED: stats.passed_tests++; break;
                case TestStatus::FAILED: stats.failed_tests++; break;
                case TestStatus::SKIPPED: stats.skipped_tests++; break;
                case TestStatus::TIMEOUT: stats.timeout_tests++; break;
                case TestStatus::ERROR: stats.error_tests++; break;
            }
            
            if (result.execution_time < fastest_time) {
                fastest_time = result.execution_time;
                stats.fastest_test = result.test_name;
                stats.fastest_time = fastest_time;
            }
            
            if (result.execution_time > slowest_time) {
                slowest_time = result.execution_time;
                stats.slowest_test = result.test_name;
                stats.slowest_time = slowest_time;
            }
        }
        
        stats.pass_rate = stats.total_tests > 0 ? 
            (static_cast<double>(stats.passed_tests) / stats.total_tests) * 100.0 : 0.0;
        
        return stats;
    }
    
    void generate_test_reports(const TestSuiteStats& stats) {
        // Generate comprehensive test reports
        // This would include HTML, XML, and JSON formats
        std::string report_path = output_directory_ + "test_report.txt";
        std::ofstream report_file(report_path);
        
        if (report_file.is_open()) {
            report_file << "HSML Unified Testing Framework Report\n";
            report_file << "=====================================\n\n";
            report_file << "Total Tests: " << stats.total_tests << "\n";
            report_file << "Passed: " << stats.passed_tests << "\n";
            report_file << "Failed: " << stats.failed_tests << "\n";
            report_file << "Skipped: " << stats.skipped_tests << "\n";
            report_file << "Timeout: " << stats.timeout_tests << "\n";
            report_file << "Error: " << stats.error_tests << "\n";
            report_file << "Pass Rate: " << stats.pass_rate << "%\n";
            report_file << "Total Execution Time: " << stats.total_execution_time.count() << "ms\n";
            report_file << "Fastest Test: " << stats.fastest_test << " (" << stats.fastest_time.count() << "ms)\n";
            report_file << "Slowest Test: " << stats.slowest_test << " (" << stats.slowest_time.count() << "ms)\n\n";
            
            // Detailed test results
            report_file << "Detailed Results:\n";
            report_file << "-----------------\n";
            for (const auto& result : test_results_) {
                report_file << result.test_name << ": ";
                switch (result.status) {
                    case TestStatus::PASSED: report_file << "PASSED"; break;
                    case TestStatus::FAILED: report_file << "FAILED"; break;
                    case TestStatus::SKIPPED: report_file << "SKIPPED"; break;
                    case TestStatus::TIMEOUT: report_file << "TIMEOUT"; break;
                    case TestStatus::ERROR: report_file << "ERROR"; break;
                }
                report_file << " (" << result.execution_time.count() << "ms)\n";
                if (!result.message.empty()) {
                    report_file << "  Message: " << result.message << "\n";
                }
                if (!result.benchmarks.empty()) {
                    report_file << "  Benchmarks:\n";
                    for (const auto& benchmark : result.benchmarks) {
                        report_file << "    " << benchmark.metric_name << ": " 
                                   << benchmark.value << " " << benchmark.unit
                                   << (benchmark.meets_baseline ? " [PASS]" : " [FAIL]") << "\n";
                    }
                }
                report_file << "\n";
            }
            
            report_file.close();
        }
    }
};

// Convenience functions for common test scenarios
inline TestSuiteStats run_critical_tests() {
    UnifiedTestingFramework framework;
    return framework.run_category_tests(TestCategory::UNIT);
}

inline TestSuiteStats run_performance_tests() {
    UnifiedTestingFramework framework;
    return framework.run_category_tests(TestCategory::PERFORMANCE);
}

inline TestSuiteStats run_integration_tests() {
    UnifiedTestingFramework framework;
    return framework.run_category_tests(TestCategory::INTEGRATION);
}

} // namespace hsml::testing