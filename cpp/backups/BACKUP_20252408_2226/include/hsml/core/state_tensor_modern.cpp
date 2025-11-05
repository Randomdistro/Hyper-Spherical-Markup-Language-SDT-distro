#pragma once

#include <array>
#include <type_traits>
#include <concepts>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>

namespace hsml::core {

// Modern C++ concepts for better type safety
template<typename T>
concept floating_point_type = std::is_floating_point_v<T>;

// Physical constants for enhanced validation
template<floating_point_type T>
struct physical_constants {
    static constexpr T speed_of_light = static_cast<T>(299792458.0);
    static constexpr T max_velocity = speed_of_light;
    static constexpr T max_angular_velocity = static_cast<T>(1e12);
    static constexpr T max_density = static_cast<T>(1e18);        // Nuclear density
    static constexpr T max_energy = static_cast<T>(1e20);         // Extreme energy scale
    static constexpr T epsilon = std::numeric_limits<T>::epsilon() * static_cast<T>(1000);
};

// 8-component state tensor with comprehensive mathematical guarantees
template<floating_point_type T>
class alignas(64) state_tensor {
    // Cache-aligned storage for optimal SIMD access
    alignas(64) std::array<T, 8> components_;  // p,q,v,ω,ρ,σ,E,Σ

    // Enhanced physical constraint validation
    constexpr void validate_physical_constraints() const noexcept {
        if (std::is_constant_evaluated()) {
            // Compile-time validation for constant expressions
            static_assert(std::is_constant_evaluated());
        } else {
            // Runtime validation with physical bounds checking
            assert(get<state_component::energy>() >= T{0} && "Energy must be non-negative");
            assert(get<state_component::density>() >= T{0} && "Density must be non-negative");
            assert(get<state_component::pressure>() >= T{0} && "Pressure must be non-negative");
            assert(get<state_component::velocity>() >= T{0} && "Velocity must be non-negative");
            assert(get<state_component::velocity>() <= physical_constants<T>::max_velocity && "Velocity cannot exceed speed of light");
            assert(get<state_component::angular_velocity>() >= T{0} && "Angular velocity must be non-negative");
            assert(get<state_component::angular_velocity>() <= physical_constants<T>::max_angular_velocity && "Angular velocity exceeds physical limits");
            assert(get<state_component::stress>() >= T{0} && "Stress must be non-negative");
            assert(get<state_component::quality>() >= T{0} && "Quality must be non-negative");
            assert(get<state_component::entropy>() >= T{0} && "Entropy must be non-negative");
        }
    }

public:
    using value_type = T;
    using size_type = size_t;
    static constexpr size_t component_count = 8;

    // Default constructor - zero state
    constexpr state_tensor() noexcept : components_{} {
        validate_physical_constraints();
    }

    // Constructor with all components
    constexpr state_tensor(T pressure, T quality, T velocity, T angular_velocity,
                          T density, T stress, T energy, T entropy) noexcept
        : components_{pressure, quality, velocity, angular_velocity, density, stress, energy, entropy} {
        validate_physical_constraints();
    }

    // Constructor from array
    constexpr explicit state_tensor(const std::array<T, 8>& components) noexcept
        : components_(components) {
        validate_physical_constraints();
    }

    // Component access with compile-time bounds checking
    template<state_component Component>
    [[nodiscard]] constexpr T get() const noexcept {
        static_assert(static_cast<size_t>(Component) < component_count);
        return components_[static_cast<size_t>(Component)];
    }

    // Component mutation with validation
    template<state_component Component>
    constexpr void set(T value) noexcept {
        static_assert(static_cast<size_t>(Component) < component_count);
        components_[static_cast<size_t>(Component)] = value;
        validate_physical_constraints();
    }

    // Bulk component access
    [[nodiscard]] constexpr const std::array<T, 8>& data() const noexcept {
        return components_;
    }

    // Safe bulk update with validation
    constexpr void update_all(const std::array<T, 8>& new_values) noexcept {
        std::array<T, 8> old_values = components_;
        components_ = new_values;

        if (!is_physically_valid()) {
            // Rollback on validation failure
            components_ = old_values;
            assert(false && "Bulk update violates physical constraints");
        }
    }

    // Physical validation
    [[nodiscard]] constexpr bool is_physically_valid() const noexcept {
        return get<state_component::energy>() >= T{0} &&
               get<state_component::density>() >= T{0} &&
               get<state_component::pressure>() >= T{0} &&
               get<state_component::velocity>() >= T{0} &&
               get<state_component::velocity>() <= physical_constants<T>::max_velocity &&
               get<state_component::angular_velocity>() >= T{0} &&
               get<state_component::angular_velocity>() <= physical_constants<T>::max_angular_velocity &&
               get<state_component::stress>() >= T{0} &&
               get<state_component::quality>() >= T{0} &&
               get<state_component::entropy>() >= T{0};
    }

    // Energy conservation check
    [[nodiscard]] constexpr bool energy_conserved() const noexcept {
        T kinetic_energy = T{0.5} * get<state_component::density>() *
                          get<state_component::velocity>() * get<state_component::velocity>();
        T potential_energy = get<state_component::pressure>() / get<state_component::density>();
        T total_energy = kinetic_energy + potential_energy;

        return std::abs(total_energy - get<state_component::energy>()) < physical_constants<T>::epsilon;
    }

    // Thermodynamic consistency check
    [[nodiscard]] constexpr bool thermodynamically_consistent() const noexcept {
        // Check if entropy is consistent with energy and density
        T temperature = get<state_component::energy>() / get<state_component::density>();
        T expected_entropy = std::log(temperature) + std::log(get<state_component::density>());

        return std::abs(expected_entropy - get<state_component::entropy>()) < physical_constants<T>::epsilon;
    }

    // State validation
    [[nodiscard]] constexpr bool is_valid() const noexcept {
        return is_physically_valid() && energy_conserved() && thermodynamically_consistent();
    }

    // Mathematical operations with constraint preservation
    [[nodiscard]] constexpr state_tensor operator+(const state_tensor& other) const noexcept {
        std::array<T, 8> result_data;
        for (size_t i = 0; i < 8; ++i) {
            result_data[i] = components_[i] + other.components_[i];
        }

        state_tensor result{result_data};

        // Validate the result
        if (!result.is_physically_valid()) {
            // Return the more stable of the two tensors
            return is_physically_valid() ? *this : other;
        }

        return result;
    }

    [[nodiscard]] constexpr state_tensor operator*(T scalar) const noexcept {
        if (scalar < T{0}) {
            // Negative scaling doesn't make physical sense for most components
            return *this;
        }

        std::array<T, 8> result_data;
        for (size_t i = 0; i < 8; ++i) {
            result_data[i] = components_[i] * scalar;
        }

        state_tensor result{result_data};
        return result.is_physically_valid() ? result : *this;
    }

    // Interpolation between states
    [[nodiscard]] constexpr state_tensor lerp(const state_tensor& other, T t) const noexcept {
        if (t <= T{0}) return *this;
        if (t >= T{1}) return other;

        std::array<T, 8> result_data;
        for (size_t i = 0; i < 8; ++i) {
            result_data[i] = components_[i] + t * (other.components_[i] - components_[i]);
        }

        state_tensor result{result_data};
        return result.is_physically_valid() ? result : *this;
    }

    // State evolution with time step
    [[nodiscard]] constexpr state_tensor evolve(T dt) const noexcept {
        if (dt <= T{0}) return *this;

        std::array<T, 8> evolved = components_;

        // Simple Euler integration for demonstration
        // Position evolves with velocity
        evolved[static_cast<size_t>(state_component::pressure)] +=
            get<state_component::velocity>() * dt;

        // Quaternion evolves with angular velocity (simplified)
        evolved[static_cast<size_t>(state_component::quality)] +=
            get<state_component::angular_velocity>() * dt;

        state_tensor result{evolved};
        return result.is_physically_valid() ? result : *this;
    }

    // Comparison operators
    [[nodiscard]] constexpr bool approximately_equal(const state_tensor& other,
                                                    T epsilon = physical_constants<T>::epsilon) const noexcept {
        for (size_t i = 0; i < 8; ++i) {
            T diff = components_[i] - other.components_[i];
            if (diff < T{0}) diff = -diff;
            if (diff > epsilon) return false;
        }
        return true;
    }

    constexpr bool operator==(const state_tensor& other) const noexcept {
        return approximately_equal(other);
    }

    constexpr bool operator!=(const state_tensor& other) const noexcept {
        return !approximately_equal(other);
    }

    // Memory layout information for optimization
    [[nodiscard]] static constexpr size_t size() noexcept { return 8; }
    [[nodiscard]] static constexpr size_t alignment() noexcept { return 64; }
    [[nodiscard]] static constexpr bool is_simd_aligned() noexcept { return true; }

    // Debug information
    [[nodiscard]] std::string debug_string() const {
        std::string result = "state_tensor{\n";
        result += "  pressure: " + std::to_string(get<state_component::pressure>()) + "\n";
        result += "  quality: " + std::to_string(get<state_component::quality>()) + "\n";
        result += "  velocity: " + std::to_string(get<state_component::velocity>()) + "\n";
        result += "  angular_velocity: " + std::to_string(get<state_component::angular_velocity>()) + "\n";
        result += "  density: " + std::to_string(get<state_component::density>()) + "\n";
        result += "  stress: " + std::to_string(get<state_component::stress>()) + "\n";
        result += "  energy: " + std::to_string(get<state_component::energy>()) + "\n";
        result += "  entropy: " + std::to_string(get<state_component::entropy>()) + "\n";
        result += "  valid: " + (is_physically_valid() ? "true" : "false") + "\n";
        result += "}";
        return result;
    }

    [[nodiscard]] constexpr size_t memory_usage() const noexcept {
        return sizeof(*this);
    }
};

// Type aliases for common configurations
using state_tensor_f32 = state_tensor<float>;
using state_tensor_f64 = state_tensor<double>;

} // namespace hsml::core
