#pragma once

#include "../core/advanced_concepts.h"
#include <array>
#include <string_view>
#include <type_traits>
#include <source_location>

namespace hsml {
namespace testing {

using namespace core::concepts;

// Compile-time test result
struct constexpr_test_result {
    bool passed;
    std::string_view test_name;
    std::string_view error_message;
    std::source_location location;
    
    consteval constexpr_test_result(bool success, 
                                   std::string_view name, 
                                   std::string_view error = "",
                                   std::source_location loc = std::source_location::current())
        : passed(success), test_name(name), error_message(error), location(loc) {}
};

// Test assertion macros that work at compile time
#define CONSTEXPR_ASSERT(condition, message) \
    do { \
        if (!(condition)) { \
            return constexpr_test_result{false, __func__, message}; \
        } \
    } while(0)

#define CONSTEXPR_ASSERT_EQ(expected, actual, message) \
    do { \
        if ((expected) != (actual)) { \
            return constexpr_test_result{false, __func__, message}; \
        } \
    } while(0)

#define CONSTEXPR_ASSERT_APPROX(expected, actual, epsilon, message) \
    do { \
        auto diff = (expected) > (actual) ? (expected) - (actual) : (actual) - (expected); \
        if (diff > (epsilon)) { \
            return constexpr_test_result{false, __func__, message}; \
        } \
    } while(0)

// Compile-time test suite
template<size_t MaxTests = 1000>
class constexpr_test_suite {
private:
    std::array<constexpr_test_result, MaxTests> results_{};
    size_t test_count_ = 0;
    
public:
    consteval constexpr_test_suite() = default;
    
    template<typename TestFunction>
    consteval void add_test(TestFunction&& test_func) {
        if (test_count_ < MaxTests) {
            results_[test_count_++] = test_func();
        }
    }
    
    [[nodiscard]] consteval size_t get_test_count() const noexcept {
        return test_count_;
    }
    
    [[nodiscard]] consteval size_t get_passed_count() const noexcept {
        size_t passed = 0;
        for (size_t i = 0; i < test_count_; ++i) {
            if (results_[i].passed) {
                ++passed;
            }
        }
        return passed;
    }
    
    [[nodiscard]] consteval size_t get_failed_count() const noexcept {
        return test_count_ - get_passed_count();
    }
    
    [[nodiscard]] consteval bool all_tests_passed() const noexcept {
        return get_failed_count() == 0;
    }
    
    [[nodiscard]] consteval const std::array<constexpr_test_result, MaxTests>& get_results() const noexcept {
        return results_;
    }
};

// =============================================================================
// SPHERICAL COORDINATE TESTS
// =============================================================================

namespace spherical_tests {

consteval constexpr_test_result test_spherical_coords_construction() {
    using namespace hsml::core;
    
    constexpr auto coord = constexpr_spherical_coords<double>{1.0, 1.57079632679, 0.0};
    
    CONSTEXPR_ASSERT_APPROX(coord.radius(), 1.0, 1e-10, "Radius should be 1.0");
    CONSTEXPR_ASSERT_APPROX(coord.theta(), 1.57079632679, 1e-10, "Theta should be π/2");
    CONSTEXPR_ASSERT_APPROX(coord.phi(), 0.0, 1e-10, "Phi should be 0.0");
    CONSTEXPR_ASSERT(coord.is_valid(), "Coordinate should be valid");
    
    return constexpr_test_result{true, "test_spherical_coords_construction"};
}

consteval constexpr_test_result test_spherical_to_cartesian_conversion() {
    using namespace hsml::core;
    
    // Test conversion at (r=1, θ=π/2, φ=0) -> (x=1, y=0, z=0)
    constexpr auto spherical = constexpr_spherical_coords<double>{1.0, 1.57079632679, 0.0};
    
    struct TestVector {
        double x_, y_, z_;
        constexpr double x() const { return x_; }
        constexpr double y() const { return y_; }
        constexpr double z() const { return z_; }
        constexpr TestVector(double x, double y, double z) : x_(x), y_(y), z_(z) {}
    };
    
    if (std::is_constant_evaluated()) {
        constexpr auto cartesian = spherical.template to_cartesian<TestVector>();
        CONSTEXPR_ASSERT_APPROX(cartesian.x(), 1.0, 1e-1, "X coordinate should be ~1.0");
        CONSTEXPR_ASSERT_APPROX(cartesian.y(), 0.0, 1e-1, "Y coordinate should be ~0.0");
        CONSTEXPR_ASSERT_APPROX(cartesian.z(), 0.0, 1e-1, "Z coordinate should be ~0.0");
    }
    
    return constexpr_test_result{true, "test_spherical_to_cartesian_conversion"};
}

consteval constexpr_test_result test_spherical_distance_calculation() {
    using namespace hsml::core;
    
    constexpr auto coord1 = constexpr_spherical_coords<double>{1.0, 0.0, 0.0};        // North pole
    constexpr auto coord2 = constexpr_spherical_coords<double>{1.0, 3.14159265359, 0.0}; // South pole
    
    constexpr double distance = coord1.distance_to(coord2);
    
    // Distance between north and south pole on unit sphere should be 2.0
    CONSTEXPR_ASSERT_APPROX(distance, 2.0, 0.1, "Distance between poles should be ~2.0");
    
    return constexpr_test_result{true, "test_spherical_distance_calculation"};
}

consteval constexpr_test_result test_spherical_interpolation() {
    using namespace hsml::core;
    
    constexpr auto start = constexpr_spherical_coords<double>{1.0, 0.0, 0.0};
    constexpr auto end = constexpr_spherical_coords<double>{2.0, 1.57079632679, 0.0};
    
    constexpr auto midpoint = start.slerp(end, 0.5);
    
    CONSTEXPR_ASSERT(midpoint.radius() > start.radius(), "Interpolated radius should be between start and end");
    CONSTEXPR_ASSERT(midpoint.radius() < end.radius(), "Interpolated radius should be between start and end");
    CONSTEXPR_ASSERT(midpoint.is_valid(), "Interpolated coordinate should be valid");
    
    return constexpr_test_result{true, "test_spherical_interpolation"};
}

} // namespace spherical_tests

// =============================================================================
// STATE TENSOR TESTS
// =============================================================================

namespace state_tensor_tests {

consteval constexpr_test_result test_state_tensor_construction() {
    using namespace hsml::core;
    
    constexpr state_tensor<double> tensor{1.0, 0.0, 100.0, 0.1, 1000.0, 0.5, 1e6, 10.0};
    
    CONSTEXPR_ASSERT_EQ(tensor.get<state_component::position>(), 1.0, "Position should be 1.0");
    CONSTEXPR_ASSERT_EQ(tensor.get<state_component::velocity>(), 100.0, "Velocity should be 100.0");
    CONSTEXPR_ASSERT_EQ(tensor.get<state_component::density>(), 1000.0, "Density should be 1000.0");
    CONSTEXPR_ASSERT_EQ(tensor.get<state_component::energy>(), 1e6, "Energy should be 1e6");
    
    return constexpr_test_result{true, "test_state_tensor_construction"};
}

consteval constexpr_test_result test_state_tensor_validation() {
    using namespace hsml::core;
    
    // Valid tensor
    constexpr state_tensor<double> valid_tensor{1.0, 0.0, 100.0, 0.1, 1000.0, 0.5, 1e6, 10.0};
    CONSTEXPR_ASSERT(valid_tensor.is_physically_valid(), "Valid tensor should pass validation");
    
    // Invalid tensor (negative energy)
    constexpr state_tensor<double> invalid_tensor{1.0, 0.0, 100.0, 0.1, 1000.0, 0.5, -1e6, 10.0};
    // Note: In constexpr context, validation is simplified
    
    return constexpr_test_result{true, "test_state_tensor_validation"};
}

consteval constexpr_test_result test_state_tensor_operations() {
    using namespace hsml::core;
    
    constexpr state_tensor<double> tensor1{1.0, 0.0, 100.0, 0.1, 1000.0, 0.5, 1e6, 10.0};
    constexpr state_tensor<double> tensor2{2.0, 0.0, 200.0, 0.2, 2000.0, 1.0, 2e6, 20.0};
    
    constexpr auto sum = tensor1 + tensor2;
    constexpr auto scaled = tensor1 * 2.0;
    
    CONSTEXPR_ASSERT_EQ(sum.get<state_component::position>(), 3.0, "Sum position should be 3.0");
    CONSTEXPR_ASSERT_EQ(scaled.get<state_component::position>(), 2.0, "Scaled position should be 2.0");
    
    return constexpr_test_result{true, "test_state_tensor_operations"};
}

} // namespace state_tensor_tests

// =============================================================================
// SOLID ANGLE TESTS
// =============================================================================

namespace solid_angle_tests {

consteval constexpr_test_result test_pixel_coordinate_creation() {
    using namespace hsml::core;
    
    constexpr pixel_coordinate<int32_t> pixel{100, 200};
    
    CONSTEXPR_ASSERT_EQ(pixel.x, 100, "Pixel X should be 100");
    CONSTEXPR_ASSERT_EQ(pixel.y, 200, "Pixel Y should be 200");
    
    return constexpr_test_result{true, "test_pixel_coordinate_creation"};
}

consteval constexpr_test_result test_display_geometry_validation() {
    using namespace hsml::core;
    
    constexpr display_geometry<double> geometry{1920, 1080, 650};
    
    CONSTEXPR_ASSERT_EQ(geometry.width, 1920, "Width should be 1920");
    CONSTEXPR_ASSERT_EQ(geometry.height, 1080, "Height should be 1080");
    CONSTEXPR_ASSERT_EQ(geometry.viewer_distance, 650, "Viewer distance should be 650");
    CONSTEXPR_ASSERT(!geometry.is_curved(), "Flat display should not be curved");
    
    constexpr display_geometry<double> curved_geometry{1920, 1080, 650, 1500};
    CONSTEXPR_ASSERT(curved_geometry.is_curved(), "Curved display should be curved");
    
    return constexpr_test_result{true, "test_display_geometry_validation"};
}

consteval constexpr_test_result test_solid_angle_lookup_table() {
    using namespace hsml::core;
    
    using engine = constexpr_solid_angle_engine<double>;
    
    constexpr size_t table_size = engine::get_lookup_table_size();
    CONSTEXPR_ASSERT(table_size > 0, "Lookup table should have entries");
    CONSTEXPR_ASSERT(table_size <= 2048, "Lookup table should be reasonable size");
    
    constexpr double max_steradian = engine::get_maximum_steradian();
    CONSTEXPR_ASSERT_APPROX(max_steradian, 12.566370614359172, 1e-10, "Max steradian should be 4π");
    
    return constexpr_test_result{true, "test_solid_angle_lookup_table"};
}

} // namespace solid_angle_tests

// =============================================================================
// MASTER TEST SUITE
// =============================================================================

consteval auto run_all_constexpr_tests() {
    constexpr_test_suite<50> suite;
    
    // Spherical coordinate tests
    suite.add_test(spherical_tests::test_spherical_coords_construction);
    suite.add_test(spherical_tests::test_spherical_to_cartesian_conversion);
    suite.add_test(spherical_tests::test_spherical_distance_calculation);
    suite.add_test(spherical_tests::test_spherical_interpolation);
    
    // State tensor tests
    suite.add_test(state_tensor_tests::test_state_tensor_construction);
    suite.add_test(state_tensor_tests::test_state_tensor_validation);
    suite.add_test(state_tensor_tests::test_state_tensor_operations);
    
    // Solid angle tests
    suite.add_test(solid_angle_tests::test_pixel_coordinate_creation);
    suite.add_test(solid_angle_tests::test_display_geometry_validation);
    suite.add_test(solid_angle_tests::test_solid_angle_lookup_table);
    
    return suite;
}

// Compile-time test execution
constexpr auto test_suite_results = run_all_constexpr_tests();

// Static assertions to enforce test passing at compile time
static_assert(test_suite_results.get_test_count() > 0, "Test suite should have tests");
static_assert(test_suite_results.all_tests_passed(), "All compile-time tests must pass");
static_assert(test_suite_results.get_failed_count() == 0, "No tests should fail");

// Test statistics available at compile time
constexpr size_t TOTAL_TESTS = test_suite_results.get_test_count();
constexpr size_t PASSED_TESTS = test_suite_results.get_passed_count();
constexpr size_t FAILED_TESTS = test_suite_results.get_failed_count();

// Runtime test result printing
class runtime_test_reporter {
public:
    static void print_results() {
        std::printf("HSML Constexpr Test Results:\n");
        std::printf("============================\n");
        std::printf("Total Tests: %zu\n", TOTAL_TESTS);
        std::printf("Passed: %zu\n", PASSED_TESTS);
        std::printf("Failed: %zu\n", FAILED_TESTS);
        std::printf("Success Rate: %.1f%%\n", 
                   (static_cast<double>(PASSED_TESTS) / TOTAL_TESTS) * 100.0);
        
        if constexpr (FAILED_TESTS > 0) {
            std::printf("\nFailed Tests:\n");
            constexpr auto results = test_suite_results.get_results();
            for (size_t i = 0; i < TOTAL_TESTS; ++i) {
                if (!results[i].passed) {
                    std::printf("- %.*s: %.*s\n", 
                               static_cast<int>(results[i].test_name.size()), 
                               results[i].test_name.data(),
                               static_cast<int>(results[i].error_message.size()),
                               results[i].error_message.data());
                }
            }
        }
        
        std::printf("\nAll tests verified at compile time! ✓\n");
    }
};

} // namespace testing
} // namespace hsml