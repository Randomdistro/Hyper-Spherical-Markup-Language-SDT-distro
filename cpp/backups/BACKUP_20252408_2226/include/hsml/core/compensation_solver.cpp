#pragma once

#include <vector>
#include <array>
#include <memory>
#include "vector3.h"
#include "matrix4.h"

namespace hsml::core {

/**
 * @brief Moore-Penrose Pseudoinverse Compensation Solver for HCS-21 Translations
 * 
 * Implements the critical S-invariant conservation algorithm for legal translations
 * in the HCS-21 21-dimensional hierarchical coordinate system.
 * 
 * Mathematical Foundation:
 * For a translation δ at level k, we need compensations {Δⱼ} at levels j > k such that:
 * Σⱼ₌ₖ⁶ wⱼ²(‖vⱼ + Δⱼ‖² - ‖vⱼ‖²) = 0
 * 
 * Linearized constraint: Σⱼ₌ₖ₊₁⁶ 2wⱼ² vⱼ·Δⱼ = -wₖ²(2vₖ·δ + ‖δ‖²)
 * Solved via Moore-Penrose pseudoinverse: Δ = J†b
 */
class CompensationSolver {
public:
    static constexpr size_t NUM_LEVELS = 7;
    static constexpr double DEFAULT_RIDGE_LAMBDA = 1e-12;
    static constexpr int DEFAULT_MAX_ITERATIONS = 10;
    static constexpr double DEFAULT_CONVERGENCE_TOLERANCE = 1e-12;

    /**
     * @brief Configuration for compensation solver
     */
    struct SolverConfig {
        double ridge_lambda = DEFAULT_RIDGE_LAMBDA;     // Ridge regularization parameter
        int max_iterations = DEFAULT_MAX_ITERATIONS;    // Maximum iterations for refinement
        double convergence_tolerance = DEFAULT_CONVERGENCE_TOLERANCE; // Convergence threshold
        bool use_iterative_refinement = true;          // Enable iterative refinement for large deltas
        bool use_proportional_fallback = true;         // Fallback to proportional distribution if singular
    };

    /**
     * @brief Result of compensation calculation
     */
    struct CompensationResult {
        std::vector<Vector3> compensation_deltas;       // Computed compensation vectors
        double residual_error;                          // Remaining S-conservation error
        int iterations_used;                            // Number of iterations performed
        bool converged;                                 // Whether algorithm converged
        bool used_fallback;                            // Whether fallback method was used
        std::string error_message;                      // Error description if failed
    };

    /**
     * @brief Solve for optimal compensation deltas using Moore-Penrose pseudoinverse
     * 
     * @param translation_level Level at which translation is applied (0-6)
     * @param delta Translation vector applied at translation_level
     * @param current_vectors Current state vectors for all levels
     * @param weights Scale weights for all levels
     * @param config Solver configuration
     * @return CompensationResult containing solution and diagnostics
     */
    static CompensationResult solve_compensation(
        size_t translation_level,
        const Vector3& delta,
        const std::array<Vector3, NUM_LEVELS>& current_vectors,
        const std::array<double, NUM_LEVELS>& weights
    );
    
    static CompensationResult solve_compensation(
        size_t translation_level,
        const Vector3& delta,
        const std::array<Vector3, NUM_LEVELS>& current_vectors,
        const std::array<double, NUM_LEVELS>& weights,
        const SolverConfig& config
    );

    /**
     * @brief Validate that computed compensation preserves S-invariant
     * 
     * @param translation_level Level at which translation is applied
     * @param delta Translation vector
     * @param compensation_deltas Computed compensation vectors
     * @param current_vectors Current state vectors
     * @param weights Scale weights
     * @param tolerance Acceptable error tolerance
     * @return True if S-conservation is satisfied within tolerance
     */
    static bool validate_s_conservation(
        size_t translation_level,
        const Vector3& delta,
        const std::vector<Vector3>& compensation_deltas,
        const std::array<Vector3, NUM_LEVELS>& current_vectors,
        const std::array<double, NUM_LEVELS>& weights,
        double tolerance = DEFAULT_CONVERGENCE_TOLERANCE
    );

private:
    // Internal matrix operations for pseudoinverse calculation
    
    /**
     * @brief Build Jacobian matrix for linearized constraint
     * J = [2w_{k+1}²v_{k+1}, 2w_{k+2}²v_{k+2}, ..., 2w_6²v_6] (flattened to row vector)
     */
    static std::vector<double> build_jacobian_matrix(
        size_t translation_level,
        const std::array<Vector3, NUM_LEVELS>& current_vectors,
        const std::array<double, NUM_LEVELS>& weights
    );

    /**
     * @brief Build constraint vector b = -w_k²(2v_k·δ + ‖δ‖²)
     */
    static double build_constraint_vector(
        size_t translation_level,
        const Vector3& delta,
        const std::array<Vector3, NUM_LEVELS>& current_vectors,
        const std::array<double, NUM_LEVELS>& weights
    );

    /**
     * @brief Compute Moore-Penrose pseudoinverse using SVD decomposition
     * For 1×n matrix J: J† = J^T (JJ^T + λI)^(-1)
     */
    static std::vector<double> compute_pseudoinverse_solution(
        const std::vector<double>& jacobian,
        double constraint_value,
        double ridge_lambda
    );

    /**
     * @brief Flatten compensation deltas from solution vector
     */
    static std::vector<Vector3> unflatten_solution(
        const std::vector<double>& solution,
        size_t num_compensation_levels
    );

    /**
     * @brief Calculate exact S-change including quadratic terms
     */
    static double calculate_exact_delta_s(
        size_t translation_level,
        const Vector3& delta,
        const std::vector<Vector3>& compensation_deltas,
        const std::array<Vector3, NUM_LEVELS>& current_vectors,
        const std::array<double, NUM_LEVELS>& weights
    );

    /**
     * @brief Proportional fallback method for singular cases
     */
    static CompensationResult solve_proportional_fallback(
        size_t translation_level,
        const Vector3& delta,
        const std::array<Vector3, NUM_LEVELS>& current_vectors,
        const std::array<double, NUM_LEVELS>& weights
    );

    /**
     * @brief Iterative refinement for large deltas (includes quadratic terms)
     */
    static CompensationResult refine_solution_iteratively(
        const CompensationResult& initial_solution,
        size_t translation_level,
        const Vector3& delta,
        const std::array<Vector3, NUM_LEVELS>& current_vectors,
        const std::array<double, NUM_LEVELS>& weights,
        const SolverConfig& config
    );
};

} // namespace hsml::core