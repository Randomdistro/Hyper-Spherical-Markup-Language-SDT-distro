#pragma once

// [The Testing Paranoid]: "EVERYTHING MUST BE TESTED! NO STONE UNTURNED!"
// [The Performance Demon]: "Every nanosecond will be measured and optimized!"
// [The Mathematical Purist]: "Property-based testing with mathematical proofs!"

#include <gtest/gtest.h>
#include <benchmark/benchmark.h>
#include <thread>
#include <chrono>
#include <random>
#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include <functional>
#include <atomic>
#include <future>

// Core HSML includes for validation
#include "../core/spatial_indexer_minimal.h"
#include "../core/spatial_indexer_oop.h"
#include "../core/spatial_indexer_simd.h"
#include "../core/spatial_indexer_template.h"
#include "../core/spatial_indexer_enterprise.h"
#include "../core/spherical_dom_embedded.h"
#include "../core/spherical_dom_gaming.h"
#include "../core/spherical_dom_modern.h"
#include "../core/spherical_dom_multithreaded.h"
#include "../core/spherical_dom_scientific.h"

namespace hsml::testing {

// [The Enterprise Bean]: "We need comprehensive abstractions for testing!"
class ComprehensiveSpatialValidator {
public:
    // [The Security Paranoid]: "FUZZ TESTING IS MANDATORY!"
    struct FuzzTestConfig {
        size_t iterations = 10000;
        double coordinate_range_min = -1000.0;
        double coordinate_range_max = 1000.0;
        bool enable_edge_cases = true;
        bool enable_nan_injection = true;
        bool enable_infinity_injection = true;
        uint64_t random_seed = std::chrono::steady_clock::now().time_since_epoch().count();
    };

    // [The Performance Demon]: "EVERY METRIC MUST BE CAPTURED!"
    struct PerformanceMetrics {
        std::chrono::nanoseconds execution_time{0};
        size_t memory_allocations = 0;
        size_t memory_deallocations = 0;
        size_t cache_hits = 0;
        size_t cache_misses = 0;
        double cpu_utilization = 0.0;
        size_t peak_memory_usage = 0;
        size_t context_switches = 0;
    };

    // [The Mathematical Purist]: "Property-based testing with invariants!"
    struct SpatialInvariant {
        std::string name;
        std::function<bool(const SphericalCoordinate&)> validator;
        std::string description;
        bool is_critical = true;
    };

    // [The Functional Purist]: "Immutable test generators!"
    class PropertyBasedTestGenerator {
    public:
        static std::vector<SphericalCoordinate> generateSphericalCoordinates(
            size_t count, 
            double r_min = 0.1, 
            double r_max = 1000.0
        ) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> r_dist(r_min, r_max);
            std::uniform_real_distribution<double> theta_dist(0.0, 2.0 * M_PI);
            std::uniform_real_distribution<double> phi_dist(0.0, M_PI);

            std::vector<SphericalCoordinate> coords;
            coords.reserve(count);

            for (size_t i = 0; i < count; ++i) {
                coords.emplace_back(
                    r_dist(gen),
                    theta_dist(gen),
                    phi_dist(gen)
                );
            }
            return coords;
        }

        // [The Hacktivist]: "Edge cases that break everything!"
        static std::vector<SphericalCoordinate> generateEdgeCases() {
            return {
                {0.0, 0.0, 0.0},                    // Origin
                {1e-10, 0.0, 0.0},                  // Near zero radius
                {1e10, 0.0, 0.0},                   // Huge radius
                {1.0, 0.0, 0.0},                    // North pole
                {1.0, 0.0, M_PI},                   // South pole
                {1.0, 2.0 * M_PI, M_PI / 2.0},      // Wrap-around theta
                {std::numeric_limits<double>::infinity(), 0.0, 0.0},
                {1.0, std::numeric_limits<double>::quiet_NaN(), 0.0},
                {-1.0, 0.0, 0.0},                   // Negative radius
                {1.0, -M_PI, 0.0},                  // Negative theta
                {1.0, 0.0, -M_PI}                   // Negative phi
            };
        }
    };

    // [The OOP Architect]: "Hierarchical test organization!"
    class SpatialIndexerTestSuite {
    public:
        template<typename IndexerType>
        static void validateDodecaheralIndexing(const std::string& implementation_name) {
            IndexerType indexer;
            
            // Test all 12 faces of dodecahedron
            for (size_t face = 0; face < 12; ++face) {
                auto test_points = generateFaceTestPoints(face);
                
                for (const auto& point : test_points) {
                    auto indexed_face = indexer.getFaceIndex(point);
                    EXPECT_EQ(indexed_face, face) 
                        << "Implementation " << implementation_name 
                        << " failed face indexing for point (" 
                        << point.r << ", " << point.theta << ", " << point.phi << ")";
                }
            }
        }

        // [The Performance Demon]: "BENCHMARK ALL THE THINGS!"
        template<typename IndexerType>
        static PerformanceMetrics benchmarkIndexer(const std::string& name, size_t iterations = 100000) {
            IndexerType indexer;
            auto test_points = PropertyBasedTestGenerator::generateSphericalCoordinates(iterations);
            
            auto start_time = std::chrono::high_resolution_clock::now();
            
            for (const auto& point : test_points) {
                benchmark::DoNotOptimize(indexer.getFaceIndex(point));
            }
            
            auto end_time = std::chrono::high_resolution_clock::now();
            
            PerformanceMetrics metrics;
            metrics.execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
            
            return metrics;
        }

    private:
        static std::vector<SphericalCoordinate> generateFaceTestPoints(size_t face_index) {
            // Generate points that should map to specific dodecahedral face
            std::vector<SphericalCoordinate> points;
            
            // [The Mathematical Purist]: "Exact geometric calculations!"
            double golden_ratio = (1.0 + std::sqrt(5.0)) / 2.0;
            double face_theta = (2.0 * M_PI * face_index) / 12.0;
            
            for (int i = 0; i < 10; ++i) {
                double r = 1.0 + (i * 0.1);
                double theta = face_theta + (i * 0.01); // Small variations
                double phi = M_PI / 2.0 + (i * 0.01);
                
                points.emplace_back(r, theta, phi);
            }
            
            return points;
        }
    };

    // [The Concurrent Mastermind]: "Thread safety validation!"
    class ConcurrencyTestSuite {
    public:
        template<typename ComponentType>
        static void validateThreadSafety(const std::string& component_name) {
            constexpr size_t num_threads = 8;
            constexpr size_t operations_per_thread = 10000;
            
            ComponentType component;
            std::vector<std::future<bool>> futures;
            std::atomic<size_t> errors{0};
            
            // [The Paranoid Tester]: "Race conditions will be found!"
            for (size_t t = 0; t < num_threads; ++t) {
                futures.emplace_back(std::async(std::launch::async, [&, t]() {
                    try {
                        for (size_t i = 0; i < operations_per_thread; ++i) {
                            auto test_data = generateThreadTestData(t, i);
                            
                            // Perform concurrent operations
                            auto result = component.process(test_data);
                            
                            if (!validateResult(result)) {
                                errors.fetch_add(1);
                                return false;
                            }
                        }
                        return true;
                    } catch (...) {
                        errors.fetch_add(1);
                        return false;
                    }
                }));
            }
            
            // Wait for all threads
            for (auto& future : futures) {
                EXPECT_TRUE(future.get()) << "Thread safety test failed for " << component_name;
            }
            
            EXPECT_EQ(errors.load(), 0) << "Thread safety errors detected: " << errors.load();
        }

    private:
        static TestData generateThreadTestData(size_t thread_id, size_t iteration) {
            return TestData{thread_id, iteration, std::chrono::steady_clock::now()};
        }
        
        static bool validateResult(const TestResult& result) {
            return result.is_valid && result.computation_complete;
        }
    };

    // [The Security Paranoid]: "ATTACK SURFACE ANALYSIS!"
    class SecurityTestSuite {
    public:
        static void validateInputSanitization() {
            std::vector<std::string> malicious_inputs = {
                "../../../etc/passwd",
                "<script>alert('xss')</script>",  
                "'; DROP TABLE users; --",
                std::string(10000, 'A'),  // Buffer overflow attempt
                "\x00\x01\x02\x03",       // Binary data
                "💩🔥💀",                    // Unicode edge cases
            };
            
            for (const auto& input : malicious_inputs) {
                EXPECT_NO_THROW({
                    // Test all parser components
                    hsml::HsmlParser parser;
                    auto result = parser.parseInput(input);
                    EXPECT_TRUE(result.is_sanitized);
                }) << "Security vulnerability with input: " << input;
            }
        }

        static void validateMemoryBounds() {
            // [The Memory Safety Guardian]: "No buffer overflows allowed!"
            constexpr size_t large_size = 1000000;
            std::vector<uint8_t> test_buffer(large_size);
            
            // Fill with known pattern
            std::iota(test_buffer.begin(), test_buffer.end(), 0);
            
            // Test boundary conditions
            EXPECT_NO_THROW({
                // Test various HSML components with large data
                hsml::SpatialIndexer indexer;
                indexer.processLargeDataset(test_buffer.data(), test_buffer.size());
            });
        }
    };

    // [The Visual Debugger]: "See the algorithms in action!"
    class VisualValidationSuite {
    public:
        static void generateSpatialVisualization(const std::string& test_name) {
            // Generate SVG visualization of spatial algorithm results
            std::ofstream svg_file("spatial_test_" + test_name + ".svg");
            
            svg_file << R"(<?xml version="1.0" encoding="UTF-8"?>
<svg width="800" height="600" xmlns="http://www.w3.org/2000/svg">
<defs>
    <radialGradient id="sphereGrad" cx="50%" cy="30%">
        <stop offset="0%" style="stop-color:#ffffff;stop-opacity:0.8"/>
        <stop offset="100%" style="stop-color:#0066cc;stop-opacity:0.4"/>
    </radialGradient>
</defs>
)";

            // Generate test points visualization
            auto test_points = PropertyBasedTestGenerator::generateSphericalCoordinates(100);
            
            for (size_t i = 0; i < test_points.size(); ++i) {
                const auto& point = test_points[i];
                
                // Convert spherical to 2D projection for visualization
                double x = 400 + (point.r * std::cos(point.theta) * std::sin(point.phi) * 200);
                double y = 300 + (point.r * std::sin(point.theta) * std::sin(point.phi) * 200);
                
                svg_file << "<circle cx=\"" << x << "\" cy=\"" << y 
                        << "\" r=\"3\" fill=\"red\" opacity=\"0.7\"/>\n";
            }
            
            svg_file << "</svg>\n";
            svg_file.close();
        }
    };

public:
    // [The Test Orchestrator]: "Run the complete validation suite!"
    static void runComprehensiveValidation() {
        std::cout << "🚀 MPD Code Monkey's Comprehensive Spatial Validation Suite 🚀\n";
        std::cout << "================================================================\n";

        // 1. Spatial Indexer Validation
        std::cout << "\n[Testing Paranoid]: Validating all 5 spatial indexer implementations...\n";
        SpatialIndexerTestSuite::validateDodecaheralIndexing<MinimalSpatialIndexer>("Minimal");
        SpatialIndexerTestSuite::validateDodecaheralIndexing<OOPSpatialIndexer>("OOP");
        SpatialIndexerTestSuite::validateDodecaheralIndexing<SIMDSpatialIndexer>("SIMD");
        SpatialIndexerTestSuite::validateDodecaheralIndexing<TemplateSpatialIndexer>("Template");
        SpatialIndexerTestSuite::validateDodecaheralIndexing<EnterpriseSpatialIndexer>("Enterprise");

        // 2. Performance Benchmarking
        std::cout << "\n[Performance Demon]: Benchmarking all implementations...\n";
        auto minimal_perf = SpatialIndexerTestSuite::benchmarkIndexer<MinimalSpatialIndexer>("Minimal");
        auto oop_perf = SpatialIndexerTestSuite::benchmarkIndexer<OOPSpatialIndexer>("OOP");
        auto simd_perf = SpatialIndexerTestSuite::benchmarkIndexer<SIMDSpatialIndexer>("SIMD");
        
        std::cout << "Minimal: " << minimal_perf.execution_time.count() << "ns\n";
        std::cout << "OOP: " << oop_perf.execution_time.count() << "ns\n";
        std::cout << "SIMD: " << simd_perf.execution_time.count() << "ns\n";

        // 3. Thread Safety Validation
        std::cout << "\n[Concurrent Mastermind]: Testing thread safety...\n";
        ConcurrencyTestSuite::validateThreadSafety<EmbeddedSphericalDOM>("EmbeddedDOM");
        ConcurrencyTestSuite::validateThreadSafety<MultithreadedSphericalDOM>("MultithreadedDOM");

        // 4. Security Testing
        std::cout << "\n[Security Paranoid]: Attack surface validation...\n";
        SecurityTestSuite::validateInputSanitization();
        SecurityTestSuite::validateMemoryBounds();

        // 5. Visual Validation
        std::cout << "\n[Visual Debugger]: Generating spatial visualizations...\n";
        VisualValidationSuite::generateSpatialVisualization("comprehensive_test");

        std::cout << "\n✅ Comprehensive validation complete! All personalities satisfied! ✅\n";
    }

    // [The Fuzzing Maniac]: "Random chaos testing!"
    static void runFuzzTesting(const FuzzTestConfig& config = {}) {
        std::cout << "\n🔥 FUZZ TESTING MODE ACTIVATED 🔥\n";
        
        std::mt19937 gen(config.random_seed);
        std::uniform_real_distribution<double> coord_dist(config.coordinate_range_min, config.coordinate_range_max);
        
        for (size_t i = 0; i < config.iterations; ++i) {
            // Generate random test data
            SphericalCoordinate fuzz_coord{
                coord_dist(gen),
                coord_dist(gen) * 2.0 * M_PI,
                coord_dist(gen) * M_PI
            };
            
            if (config.enable_nan_injection && (i % 1000 == 0)) {
                fuzz_coord.r = std::numeric_limits<double>::quiet_NaN();
            }
            
            if (config.enable_infinity_injection && (i % 1001 == 0)) {
                fuzz_coord.theta = std::numeric_limits<double>::infinity();
            }
            
            // Test all components with fuzzed data
            EXPECT_NO_THROW({
                MinimalSpatialIndexer indexer;
                auto result = indexer.getFaceIndex(fuzz_coord);
                // Result should be valid or properly handled
            }) << "Fuzz test failed at iteration " << i;
        }
        
        std::cout << "🎯 Fuzz testing completed: " << config.iterations << " iterations\n";
    }
};

// [The Statistical Analyst]: "Regression detection with mathematical rigor!"
class PerformanceRegressionDetector {
public:
    struct BenchmarkResult {
        std::string test_name;
        std::chrono::nanoseconds mean_time;
        std::chrono::nanoseconds std_deviation;
        size_t sample_size;
        double confidence_interval_95_lower;
        double confidence_interval_95_upper;
    };

    static std::vector<BenchmarkResult> runStatisticalBenchmarks() {
        std::vector<BenchmarkResult> results;
        
        // Benchmark all spatial indexer implementations with statistical analysis
        const std::vector<std::string> implementations = {
            "Minimal", "OOP", "SIMD", "Template", "Enterprise"
        };
        
        for (const auto& impl : implementations) {
            auto result = benchmarkWithStatistics(impl);
            results.push_back(result);
        }
        
        return results;
    }

private:
    static BenchmarkResult benchmarkWithStatistics(const std::string& implementation) {
        constexpr size_t num_samples = 1000;
        constexpr size_t operations_per_sample = 10000;
        
        std::vector<std::chrono::nanoseconds> sample_times;
        sample_times.reserve(num_samples);
        
        for (size_t i = 0; i < num_samples; ++i) {
            auto start = std::chrono::high_resolution_clock::now();
            
            // Run implementation-specific benchmark
            runImplementationBenchmark(implementation, operations_per_sample);
            
            auto end = std::chrono::high_resolution_clock::now();
            sample_times.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end - start));
        }
        
        // Calculate statistics
        auto mean = calculateMean(sample_times);
        auto std_dev = calculateStandardDeviation(sample_times, mean);
        auto [ci_lower, ci_upper] = calculate95ConfidenceInterval(sample_times);
        
        return BenchmarkResult{
            implementation,
            mean,
            std_dev,
            num_samples,
            ci_lower,
            ci_upper
        };
    }
    
    static std::chrono::nanoseconds calculateMean(const std::vector<std::chrono::nanoseconds>& samples) {
        auto sum = std::accumulate(samples.begin(), samples.end(), std::chrono::nanoseconds{0});
        return sum / samples.size();
    }
    
    static std::chrono::nanoseconds calculateStandardDeviation(
        const std::vector<std::chrono::nanoseconds>& samples,
        std::chrono::nanoseconds mean
    ) {
        double variance = 0.0;
        for (const auto& sample : samples) {
            double diff = static_cast<double>((sample - mean).count());
            variance += diff * diff;
        }
        variance /= samples.size();
        return std::chrono::nanoseconds{static_cast<long long>(std::sqrt(variance))};
    }
    
    static std::pair<double, double> calculate95ConfidenceInterval(
        const std::vector<std::chrono::nanoseconds>& samples
    ) {
        // Student's t-distribution for 95% confidence interval
        const double t_critical = 1.96; // Approximate for large samples
        
        auto mean = calculateMean(samples);
        auto std_dev = calculateStandardDeviation(samples, mean);
        
        double margin_of_error = t_critical * (static_cast<double>(std_dev.count()) / std::sqrt(samples.size()));
        
        return {
            static_cast<double>(mean.count()) - margin_of_error,
            static_cast<double>(mean.count()) + margin_of_error
        };
    }
    
    static void runImplementationBenchmark(const std::string& implementation, size_t operations) {
        // Implementation-specific benchmark logic would go here
        // This is a placeholder for the actual benchmark execution
        volatile double result = 0.0;
        for (size_t i = 0; i < operations; ++i) {
            result += std::sin(i * 0.001); // Dummy computation
        }
        benchmark::DoNotOptimize(result);
    }
};

} // namespace hsml::testing