#pragma once

#include <random>
#include <functional>
#include <vector>
#include <string>
#include <memory>
#include <chrono>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <type_traits>
#include <concepts>

namespace hsml {
namespace testing {

// Random number generator with reproducible seeds
class test_random_generator {
    std::mt19937_64 rng_;
    std::uniform_real_distribution<double> real_dist_;
    std::uniform_int_distribution<int> int_dist_;
    
public:
    explicit test_random_generator(uint64_t seed = 42) 
        : rng_(seed), real_dist_(-1000.0, 1000.0), int_dist_(-1000, 1000) {}
    
    void set_seed(uint64_t seed) {
        rng_.seed(seed);
    }
    
    double random_real() {
        return real_dist_(rng_);
    }
    
    int random_int() {
        return int_dist_(rng_);
    }
    
    template<typename T>
    T random_in_range(T min, T max) {
        if constexpr (std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> dist(min, max);
            return dist(rng_);
        } else {
            std::uniform_int_distribution<T> dist(min, max);
            return dist(rng_);
        }
    }
    
    // Generate random spherical coordinates
    struct spherical_coords_generator {
        test_random_generator& rng;
        
        struct spherical_coords {
            double r, theta, phi;
        };
        
        spherical_coords operator()() {
            return {
                rng.random_in_range(0.1, 1000.0),  // radius
                rng.random_in_range(0.0, M_PI),     // theta
                rng.random_in_range(-M_PI, M_PI)    // phi
            };
        }
    };
    
    spherical_coords_generator spherical_coords() {
        return {*this};
    }
    
    // Generate random vectors
    struct vector3_generator {
        test_random_generator& rng;
        
        struct vector3 {
            double x, y, z;
        };
        
        vector3 operator()() {
            return {
                rng.random_real(),
                rng.random_real(),
                rng.random_real()
            };
        }
    };
    
    vector3_generator vector3() {
        return {*this};
    }
    
    // Generate random colors
    struct color_generator {
        test_random_generator& rng;
        
        struct color {
            double r, g, b;
        };
        
        color operator()() {
            return {
                rng.random_in_range(0.0, 1.0),
                rng.random_in_range(0.0, 1.0),
                rng.random_in_range(0.0, 1.0)
            };
        }
    };
    
    color_generator color() {
        return {*this};
    }
};

// Test result structure
struct test_result {
    bool passed;
    std::string message;
    std::chrono::microseconds execution_time;
    size_t iteration;
    std::string test_data_description;
    
    test_result(bool p, std::string msg, std::chrono::microseconds time, size_t iter, std::string data_desc)
        : passed(p), message(std::move(msg)), execution_time(time), iteration(iter), test_data_description(std::move(data_desc)) {}
};

// Property-based tester with automatic test case generation
template<typename Generator>
class property_based_tester {
    Generator generator_;
    std::string property_name_;
    size_t max_iterations_;
    std::chrono::milliseconds timeout_;
    bool verbose_;
    
public:
    explicit property_based_tester(Generator generator, std::string property_name = "Unnamed Property")
        : generator_(std::move(generator))
        , property_name_(std::move(property_name))
        , max_iterations_(10000)
        , timeout_(std::chrono::milliseconds{5000})
        , verbose_(false) {}
    
    // Configure test parameters
    property_based_tester& max_iterations(size_t iterations) {
        max_iterations_ = iterations;
        return *this;
    }
    
    property_based_tester& timeout(std::chrono::milliseconds timeout) {
        timeout_ = timeout;
        return *this;
    }
    
    property_based_tester& verbose(bool enable) {
        verbose_ = enable;
        return *this;
    }
    
    // Test a property with automatic test case generation
    template<typename Property>
    bool test_property(Property&& prop, size_t iterations = 10000) {
        auto start_time = std::chrono::high_resolution_clock::now();
        size_t actual_iterations = std::min(iterations, max_iterations_);
        
        if (verbose_) {
            std::cout << "Testing property: " << property_name_ << std::endl;
            std::cout << "Iterations: " << actual_iterations << std::endl;
            std::cout << "Timeout: " << timeout_.count() << "ms" << std::endl;
        }
        
        for (size_t i = 0; i < actual_iterations; ++i) {
            auto iteration_start = std::chrono::high_resolution_clock::now();
            
            // Check timeout
            if (iteration_start - start_time > timeout_) {
                if (verbose_) {
                    std::cout << "Test timed out after " << i << " iterations" << std::endl;
                }
                return false;
            }
            
            // Generate test data
            auto test_data = generator_();
            
            // Test the property
            bool result = false;
            std::string error_message;
            
            try {
                result = prop(test_data);
            } catch (const std::exception& e) {
                error_message = std::string("Exception: ") + e.what();
                result = false;
            } catch (...) {
                error_message = "Unknown exception";
                result = false;
            }
            
            if (!result) {
                if (verbose_) {
                    std::cout << "Property failed at iteration " << i << std::endl;
                    std::cout << "Test data: " << describe_test_data(test_data) << std::endl;
                    if (!error_message.empty()) {
                        std::cout << "Error: " << error_message << std::endl;
                    }
                }
                return false;
            }
            
            if (verbose_ && (i + 1) % 1000 == 0) {
                auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::high_resolution_clock::now() - start_time);
                std::cout << "Completed " << (i + 1) << " iterations in " << elapsed.count() << "ms" << std::endl;
            }
        }
        
        if (verbose_) {
            auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::high_resolution_clock::now() - start_time);
            std::cout << "Property passed all " << actual_iterations << " iterations in " 
                     << total_time.count() << "ms" << std::endl;
        }
        
        return true;
    }
    
    // Test with detailed results
    template<typename Property>
    std::vector<test_result> test_property_detailed(Property&& prop, size_t iterations = 10000) {
        std::vector<test_result> results;
        size_t actual_iterations = std::min(iterations, max_iterations_);
        
        for (size_t i = 0; i < actual_iterations; ++i) {
            auto start_time = std::chrono::high_resolution_clock::now();
            
            auto test_data = generator_();
            bool passed = false;
            std::string message;
            
            try {
                passed = prop(test_data);
                message = passed ? "Passed" : "Failed";
            } catch (const std::exception& e) {
                passed = false;
                message = std::string("Exception: ") + e.what();
            } catch (...) {
                passed = false;
                message = "Unknown exception";
            }
            
            auto end_time = std::chrono::high_resolution_clock::now();
            auto execution_time = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
            
            results.emplace_back(passed, message, execution_time, i, describe_test_data(test_data));
            
            if (!passed) {
                break; // Stop on first failure
            }
        }
        
        return results;
    }
    
    // Shrink failing test cases
    template<typename Property>
    auto shrink_failing_case(Property&& prop, size_t max_shrink_attempts = 100) -> struct {
        bool found_smaller_case;
        std::string minimal_test_data;
        size_t shrink_attempts;
    } {
        // Generate initial failing case
        auto failing_data = generator_();
        if (prop(failing_data)) {
            return {false, "", 0}; // No failing case found
        }
        
        std::string current_data_desc = describe_test_data(failing_data);
        size_t attempts = 0;
        
        while (attempts < max_shrink_attempts) {
            // Try to generate a "smaller" test case
            auto smaller_data = generate_smaller_case(failing_data);
            
            if (prop(smaller_data)) {
                // Smaller case passes, keep current failing case
                break;
            } else {
                // Smaller case also fails, use it as new failing case
                failing_data = smaller_data;
                current_data_desc = describe_test_data(failing_data);
            }
            
            ++attempts;
        }
        
        return {true, current_data_desc, attempts};
    }
    
private:
    // Describe test data for debugging
    template<typename T>
    std::string describe_test_data(const T& data) {
        std::ostringstream oss;
        if constexpr (std::is_same_v<T, typename Generator::spherical_coords_generator::spherical_coords>) {
            oss << "SphericalCoords(r=" << data.r << ", theta=" << data.theta << ", phi=" << data.phi << ")";
        } else if constexpr (std::is_same_v<T, typename Generator::vector3_generator::vector3>) {
            oss << "Vector3(x=" << data.x << ", y=" << data.y << ", z=" << data.z << ")";
        } else if constexpr (std::is_same_v<T, typename Generator::color_generator::color>) {
            oss << "Color(r=" << data.r << ", g=" << data.g << ", b=" << data.b << ")";
        } else {
            oss << "Unknown test data type";
        }
        return oss.str();
    }
    
    // Generate a "smaller" test case for shrinking
    template<typename T>
    T generate_smaller_case(const T& original) {
        if constexpr (std::is_same_v<T, typename Generator::spherical_coords_generator::spherical_coords>) {
            return {
                original.r * 0.5,           // Smaller radius
                original.theta * 0.5,       // Smaller theta
                original.phi * 0.5          // Smaller phi
            };
        } else if constexpr (std::is_same_v<T, typename Generator::vector3_generator::vector3>) {
            return {
                original.x * 0.5,           // Smaller x
                original.y * 0.5,           // Smaller y
                original.z * 0.5            // Smaller z
            };
        } else if constexpr (std::is_same_v<T, typename Generator::color_generator::color>) {
            return {
                original.r * 0.5,           // Smaller red
                original.g * 0.5,           // Smaller green
                original.b * 0.5            // Smaller blue
            };
        } else {
            return original; // No shrinking for unknown types
        }
    }
};

// Predefined test generators
class test_generators {
public:
    // Spherical coordinates generator
    static auto spherical_coords(uint64_t seed = 42) {
        test_random_generator rng(seed);
        return property_based_tester(rng.spherical_coords(), "Spherical Coordinates Property");
    }
    
    // Vector3 generator
    static auto vector3(uint64_t seed = 42) {
        test_random_generator rng(seed);
        return property_based_tester(rng.vector3(), "Vector3 Property");
    }
    
    // Color generator
    static auto color(uint64_t seed = 42) {
        test_random_generator rng(seed);
        return property_based_tester(rng.color(), "Color Property");
    }
    
    // Custom generator
    template<typename CustomGenerator>
    static auto custom(CustomGenerator&& generator, std::string name = "Custom Property") {
        return property_based_tester(std::forward<CustomGenerator>(generator), std::move(name));
    }
};

// Property testing macros
#define HSML_PROPERTY_TEST(name, generator, property) \
    TEST(name, PropertyBased) { \
        auto tester = generator.verbose(true); \
        EXPECT_TRUE(tester.test_property(property)); \
    }

#define HSML_PROPERTY_TEST_DETAILED(name, generator, property) \
    TEST(name, PropertyBasedDetailed) { \
        auto tester = generator.verbose(true); \
        auto results = tester.test_property_detailed(property); \
        for (const auto& result : results) { \
            EXPECT_TRUE(result.passed) << "Failed at iteration " << result.iteration \
                                      << ": " << result.message \
                                      << " (Data: " << result.test_data_description << ")"; \
        } \
    }

// Example properties for testing
namespace example_properties {
    // Property: Spherical coordinates should be valid
    auto spherical_coords_valid = [](const test_random_generator::spherical_coords_generator::spherical_coords& coords) {
        return coords.r >= 0.0 && 
               coords.theta >= 0.0 && coords.theta <= M_PI &&
               coords.phi >= -M_PI && coords.phi <= M_PI;
    };
    
    // Property: Vector3 magnitude should be non-negative
    auto vector3_magnitude_nonnegative = [](const test_random_generator::vector3_generator::vector3& vec) {
        double magnitude = std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
        return magnitude >= 0.0;
    };
    
    // Property: Color components should be in [0,1] range
    auto color_components_in_range = [](const test_random_generator::color_generator::color& color) {
        return color.r >= 0.0 && color.r <= 1.0 &&
               color.g >= 0.0 && color.g <= 1.0 &&
               color.b >= 0.0 && color.b <= 1.0;
    };
}

// Compile-time tests
namespace compile_time_tests {
    consteval bool test_property_based_tester_compilation() {
        test_random_generator rng;
        auto tester = property_based_tester(rng.spherical_coords(), "Test");
        return true;
    }
    
    consteval bool test_generators_compilation() {
        auto spherical_tester = test_generators::spherical_coords();
        auto vector_tester = test_generators::vector3();
        auto color_tester = test_generators::color();
        return true;
    }
    
    static_assert(test_property_based_tester_compilation());
    static_assert(test_generators_compilation());
}

} // namespace testing
} // namespace hsml
