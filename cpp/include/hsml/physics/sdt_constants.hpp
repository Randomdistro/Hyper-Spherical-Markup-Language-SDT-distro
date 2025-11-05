#pragma once

#include <cmath>

namespace hsml::physics {

/**
 * SDT Fundamental Constants
 * Pure Spatial Displacement Theory - No Cartesian contamination
 */
struct SDTConstants {
    // Fundamental SDT constants
    static constexpr double K_SDT = 1.2700e-4;      // m³/kg - SDT displacement constant
    static constexpr double EPSILON = 2.3e-20;      // m³/kg - Non-linear coefficient
    static constexpr double A0 = 866.0;             // m/s² - Universal scaling factor
    static constexpr double C_SDT = 299792458.0;    // m/s - Speed of SDT wave propagation
    static constexpr double P0 = 1.0e10;            // Pa - Baseline spation pressure
    static constexpr double LAMBDA0 = 1.0e6;        // m - Base length scale
    static constexpr double M0 = 1.0e30;            // kg - Reference mass

    // Resonance constants
    static constexpr double RESONANCE_432HZ = 432.0;  // Hz - Universal resonance frequency
    static constexpr double RESONANCE_528HZ = 528.0;  // Hz - DNA repair frequency

    // Physical constants (for reference)
    static constexpr double PI = M_PI;
    static constexpr double TWO_PI = 2.0 * M_PI;
    static constexpr double FOUR_PI = 4.0 * M_PI;

    // Zero-division safety
    static constexpr double EPSILON_DISTANCE = 1e-10;  // Minimum safe distance
    static constexpr double EPSILON_MASS = 1e-20;      // Minimum safe mass
};

/**
 * SDT Field Calculations - Pure Spherical
 */
template<typename T = double>
class SDTFieldCalculations {
public:
    /**
     * Calculate spatial displacement field at a point
     * D(r) = k * M / r³ * f(r/λ) + ε * (k * M)² / r⁶ * g(r/λ)
     *
     * PURE SPHERICAL - uses only radial distance
     */
    static T displacement_field(T mass, T distance) noexcept {
        if (distance < SDTConstants::EPSILON_DISTANCE) {
            distance = SDTConstants::EPSILON_DISTANCE;
        }

        const T k = static_cast<T>(SDTConstants::K_SDT);
        const T epsilon = static_cast<T>(SDTConstants::EPSILON);
        const T lambda = static_cast<T>(SDTConstants::LAMBDA0) *
                        std::pow(mass / SDTConstants::M0, T(1)/T(3));

        // Scale functions
        const T f = T(1) - std::exp(-distance / lambda);
        const T g = (T(1) - std::exp(-distance / lambda)) * (T(1) - std::exp(-distance / lambda));

        // Primary and non-linear terms
        const T r3 = distance * distance * distance;
        const T r6 = r3 * r3;
        const T primary = (k * mass / r3) * f;
        const T nonLinear = epsilon * (k * mass) * (k * mass) / r6 * g;

        return primary + nonLinear;
    }

    /**
     * Calculate pressure gradient from displacement field
     * ∇P(r) = -α * ∇D(r) * h(ρ(r))
     */
    static T pressure_gradient(T displacement, T density, T alpha = T(1)) noexcept {
        const T density_response = std::pow(density / T(1), T(1)/T(3));
        return -alpha * displacement * density_response;
    }

    /**
     * Calculate eclipsing function for multi-body interactions
     * E(r, M₁, M₂) = (Ω₁₂ + Ω₂₁) / (4π) * [1 - exp(-r/λ_eff)] * F(M₁, M₂, r)
     */
    static T eclipsing(T r, T m1, T m2, T omega12, T omega21, T beta = T(0.1)) noexcept {
        const T lambda_eff = static_cast<T>(SDTConstants::LAMBDA0) *
                            std::pow((m1 + m2) / SDTConstants::M0, T(1)/T(3));
        const T F = T(1) + beta * (m1 * m2) / (r * SDTConstants::M0 * SDTConstants::M0);

        return (omega12 + omega21) / (T(4) * M_PI) * (T(1) - std::exp(-r / lambda_eff)) * F;
    }

    /**
     * Calculate resonance interaction strength
     * Pure harmonic calculation - no Cartesian coordinates
     */
    static T resonance_strength(T freq1, T freq2) noexcept {
        if (freq2 < SDTConstants::EPSILON_DISTANCE) {
            freq2 = SDTConstants::EPSILON_DISTANCE;
        }

        const T ratio = freq1 / freq2;
        const T harmonic = std::sin(ratio * static_cast<T>(SDTConstants::TWO_PI));

        return harmonic * harmonic;
    }
};

} // namespace hsml::physics
