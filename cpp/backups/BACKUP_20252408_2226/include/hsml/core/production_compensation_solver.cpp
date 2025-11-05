#pragma once

#include "vector3.h"
#include "matrix4.h"
#include <array>
#include <vector>
#include <string>
#include <limits>

namespace hsml::core {

/**
 * @brief Maximum Error Negation Trust-Region Compensation Solver v1.1.2
 * 
 * Implements full nonlinear constraint solving with hierarchical redundancy,
 * post-commit verification, and deterministic S-conservation to < 1e-12 tolerance.
 * 
 * Mathematical Foundation:
 * Full constraint: w_k²(2v_k·δ·step + ‖δ‖²·step²) + Σ_j w_j²(2v_j·Δ_j + ‖Δ_j‖²) = 0
 * 
 * Error Negation Features:
 * - Trust-region with step halving and relinearization
 * - Scaled Jacobian and variable normalization  
 * - Adaptive ridge + condition-number based damping
 * - Residual-based rollback with max retries
 * - Deterministic tolerance: eps = max(1e-12 * S_total, 1e-18)
 * - Hierarchical redundancy: Moore-Penrose → Constrained LS → Proportional
 */
class ProductionCompensationSolver {
public:
    static constexpr size_t NUM_LEVELS = 7;
    static constexpr int DEFAULT_MAX_OUTER_ITERATIONS = 10;
    static constexpr double MIN_STEP_SCALE = 1e-6;
    static constexpr double STEP_HALVING_FACTOR = 0.5;
    static constexpr double CONDITION_THRESHOLD = 1e12;
    
    enum class SolverMethod {
        TRUST_REGION_PRIMARY = 0,        // Moore-Penrose with trust region
        CONSTRAINED_LEAST_SQUARES = 1,   // Lagrange multiplier closed form
        PROPORTIONAL_FALLBACK = 2        // Emergency proportional distribution
    };
    
    enum class SolveStatus {
        SUCCESS = 0,
        FELL_BACK_TO_SECONDARY = 1,
        FELL_BACK_TO_EMERGENCY = 2,
        SINGULAR_JACOBIAN = 3,
        STEP_SIZE_UNDERFLOW = 4,
        POST_COMMIT_DRIFT = 5,
        CATASTROPHIC_FAILURE = 6
    };

    /**
     * @brief Maximum error negation solve result with comprehensive diagnostics
     */
    struct SolveResult {
        std::vector<Vector3> compensation_deltas;  // Δ_j for j > k
        double final_residual_error;               // |ΔS| after solve
        int iterations_used;                       // Outer iterations consumed
        SolveStatus status;                        // Final result status
        SolverMethod successful_method;            // Which method succeeded
        std::string diagnostics;                   // Human-readable details
        
        // Error negation metrics
        double condition_number;                   // Jacobian conditioning
        double step_scale_final;                  // Final step scaling factor
        int halving_count;                        // Number of step halvings
        bool used_fallback;                       // Whether fallback was needed
        
        // Performance tracking
        double solve_time_ms;                     // Total solve time
        std::vector<double> residual_history;     // Per-iteration residuals
    };

    /**
     * @brief Primary trust-region solver with full nonlinear constraint handling
     * 
     * Implements the core algorithm:
     * 1. Build scaled Jacobian J_normalized = J / max(|J|)
     * 2. Solve J·x = b via Moore-Penrose with adaptive ridge
     * 3. Apply step with trust-region scaling (step halving)
     * 4. Evaluate |ΔS| residual and backtrack if necessary
     * 5. Relinearize and repeat until convergence
     */
    static SolveResult solve_translate(
        size_t translation_level,                     // k: level being moved
        const Vector3& delta,                         // δ: desired translation
        const std::array<Vector3, NUM_LEVELS>& v,     // Current state vectors
        const std::array<double, NUM_LEVELS>& w,      // Level weights
        double S_before,                              // S-invariant before move
        double epsilon_absolute,                      // |ΔS| tolerance
        int max_outer_iterations = DEFAULT_MAX_OUTER_ITERATIONS
    );

    /**
     * @brief Constrained least-squares fallback (Lagrange multiplier closed form)
     * 
     * When Moore-Penrose fails due to singularity, uses analytical solution:
     * min ||x||² subject to J·x = b via Lagrangian method
     */
    static SolveResult solve_constrained_least_squares(
        size_t translation_level,
        const Vector3& delta,
        const std::array<Vector3, NUM_LEVELS>& v,
        const std::array<double, NUM_LEVELS>& w,
        double S_before,
        double epsilon_absolute
    );

    /**
     * @brief Emergency proportional fallback 
     * 
     * Distributes S-change proportionally by w_j² when all else fails.
     * Uses relaxed tolerance (100x normal) for emergency conditions.
     */
    static SolveResult solve_proportional_emergency(
        size_t translation_level,
        const Vector3& delta,
        const std::array<Vector3, NUM_LEVELS>& v,
        const std::array<double, NUM_LEVELS>& w,
        double S_before,
        double epsilon_absolute
    );

    /**
     * @brief Calculate deterministic tolerance: eps = max(1e-12 * S_total, 1e-18)
     */
    static double calculate_deterministic_tolerance(
        const std::array<Vector3, NUM_LEVELS>& v,
        const std::array<double, NUM_LEVELS>& w
    );

    /**
     * @brief Evaluate S-invariant with compensation deltas applied
     */
    static double evaluate_S_with_deltas(
        size_t translation_level,
        const Vector3& delta,
        const std::vector<Vector3>& compensation_deltas,
        const std::array<Vector3, NUM_LEVELS>& v,
        const std::array<double, NUM_LEVELS>& w
    );

    /**
     * @brief Calculate current S-invariant value
     */
    static double calculate_S_total(
        const std::array<Vector3, NUM_LEVELS>& v,
        const std::array<double, NUM_LEVELS>& w
    );

    /**
     * @brief Execute translation with post-commit verification and rollback
     * 
     * This is the main entry point that combines solve + apply + verify + rollback.
     * Implements maximum error negation by automatically rolling back if post-commit
     * S-conservation is violated.
     */
    static SolveResult translate_with_verification(
        size_t translation_level,
        const Vector3& delta,
        std::array<Vector3, NUM_LEVELS>& v,  // Modified in-place
        const std::array<double, NUM_LEVELS>& w,
        int max_outer_iterations = DEFAULT_MAX_OUTER_ITERATIONS
    );

private:
    /**
     * @brief Build scaled Jacobian matrix for linearized constraint
     * Returns J_normalized = J / max(|J|) and the scaling factor
     */
    static std::pair<std::vector<double>, double> build_scaled_jacobian(
        size_t translation_level,
        double step_scale,
        const std::array<Vector3, NUM_LEVELS>& v,
        const std::array<double, NUM_LEVELS>& w
    );

    /**
     * @brief Calculate adaptive ridge parameter based on condition number
     */
    static double calculate_adaptive_ridge(
        const std::vector<double>& jacobian_normalized
    );

    /**
     * @brief Solve Moore-Penrose pseudoinverse: x = J^T (JJ^T + ridge)^-1 b
     */
    static std::vector<double> solve_moore_penrose(
        const std::vector<double>& jacobian_normalized,
        double constraint_rhs,
        double ridge_parameter
    );

    /**
     * @brief Convert flat solution vector back to Vector3 deltas
     */
    static std::vector<Vector3> unflatten_compensation_deltas(
        const std::vector<double>& solution,
        size_t translation_level
    );

    /**
     * @brief Check Jacobian condition number via simple norm ratio
     */
    static double estimate_condition_number(
        const std::vector<double>& jacobian
    );
};

} // namespace hsml::core