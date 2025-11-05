#include "hsml/core/compensation_solver.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace hsml::core {

// Forward declaration for helper function
double calculate_current_S(const std::array<Vector3, CompensationSolver::NUM_LEVELS>& vectors,
                          const std::array<double, CompensationSolver::NUM_LEVELS>& weights);

CompensationSolver::CompensationResult CompensationSolver::solve_compensation(
    size_t translation_level,
    const Vector3& delta,
    const std::array<Vector3, NUM_LEVELS>& current_vectors,
    const std::array<double, NUM_LEVELS>& weights) {
    
    SolverConfig config; // Use default configuration
    return solve_compensation(translation_level, delta, current_vectors, weights, config);
}

CompensationSolver::CompensationResult CompensationSolver::solve_compensation(
    size_t translation_level,
    const Vector3& delta,
    const std::array<Vector3, NUM_LEVELS>& current_vectors,
    const std::array<double, NUM_LEVELS>& weights,
    const SolverConfig& config) {
    
    CompensationResult result;
    result.converged = false;
    result.used_fallback = false;
    result.iterations_used = 0;
    result.residual_error = 0.0;

    // Validate inputs
    if (translation_level >= NUM_LEVELS) {
        result.error_message = "Invalid translation level";
        return result;
    }

    size_t num_compensation_levels = NUM_LEVELS - translation_level - 1;
    if (num_compensation_levels == 0) {
        // No lower levels to compensate - trivial case
        result.converged = true;
        return result;
    }

    try {
        // Step 1: Build linearized constraint system J·x = b
        std::vector<double> jacobian = build_jacobian_matrix(translation_level, current_vectors, weights);
        double constraint_value = build_constraint_vector(translation_level, delta, current_vectors, weights);

        // Step 2: Solve using Moore-Penrose pseudoinverse
        std::vector<double> solution = compute_pseudoinverse_solution(jacobian, constraint_value, config.ridge_lambda);
        
        // Check for numerical issues
        bool has_nan_or_inf = std::any_of(solution.begin(), solution.end(), 
            [](double x) { return std::isnan(x) || std::isinf(x); });
            
        if (has_nan_or_inf) {
            // Fallback to proportional method
            if (config.use_proportional_fallback) {
                return solve_proportional_fallback(translation_level, delta, current_vectors, weights);
            } else {
                result.error_message = "Numerical instability in pseudoinverse solution";
                return result;
            }
        }

        // Step 3: Convert solution back to Vector3 format
        result.compensation_deltas = unflatten_solution(solution, num_compensation_levels);

        // Step 4: Calculate initial residual
        result.residual_error = calculate_exact_delta_s(
            translation_level, delta, result.compensation_deltas, current_vectors, weights
        );
        
        // PATCH 3: Step halving if residual too large
        double step_size = 1.0;
        int step_attempts = 0;
        const int max_step_attempts = 5;
        
        while (std::abs(result.residual_error) > config.convergence_tolerance && step_attempts < max_step_attempts) {
            step_size *= 0.5;
            step_attempts++;
            
            // Scale down compensation deltas
            for (auto& comp_delta : result.compensation_deltas) {
                comp_delta = comp_delta * 0.5;
            }
            
            // Recalculate residual
            result.residual_error = calculate_exact_delta_s(
                translation_level, delta, result.compensation_deltas, current_vectors, weights
            );
        }
        
        if (std::abs(result.residual_error) < config.convergence_tolerance) {
            result.converged = true;
            result.iterations_used = step_attempts + 1;
        } else if (config.use_iterative_refinement) {
            // Step 5: Iterative refinement for remaining error
            result = refine_solution_iteratively(
                result, translation_level, delta, current_vectors, weights, config
            );
        } else {
            result.converged = false;
        }

    } catch (const std::exception& e) {
        result.error_message = "Exception in compensation solver: " + std::string(e.what());
        result.converged = false;
        
        // PATCH 4: Enhanced error context
        result.error_message += " [Context: level=" + std::to_string(translation_level) + 
                               ", delta_magnitude=" + std::to_string(delta.magnitude()) + 
                               ", S_before=" + std::to_string(calculate_current_S(current_vectors, weights)) + "]";
        
        // Try fallback method only if explicitly configured
        if (config.use_proportional_fallback) {
            auto fallback_result = solve_proportional_fallback(translation_level, delta, current_vectors, weights);
            fallback_result.error_message = "Primary solver failed, used fallback: " + result.error_message;
            return fallback_result;
        }
    }

    return result;
}

std::vector<double> CompensationSolver::build_jacobian_matrix(
    size_t translation_level,
    const std::array<Vector3, NUM_LEVELS>& current_vectors,
    const std::array<double, NUM_LEVELS>& weights) {
    
    size_t num_compensation_levels = NUM_LEVELS - translation_level - 1;
    std::vector<double> jacobian(num_compensation_levels * 3);
    
    size_t idx = 0;
    for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
        double w_j_squared = weights[j] * weights[j];
        Vector3 v_j = current_vectors[j];
        
        // Jacobian entries: 2 * w_j² * v_j components
        jacobian[idx * 3 + 0] = 2.0 * w_j_squared * v_j.x();
        jacobian[idx * 3 + 1] = 2.0 * w_j_squared * v_j.y();
        jacobian[idx * 3 + 2] = 2.0 * w_j_squared * v_j.z();
        
        idx++;
    }
    
    return jacobian;
}

double CompensationSolver::build_constraint_vector(
    size_t translation_level,
    const Vector3& delta,
    const std::array<Vector3, NUM_LEVELS>& current_vectors,
    const std::array<double, NUM_LEVELS>& weights) {
    
    double w_k_squared = weights[translation_level] * weights[translation_level];
    Vector3 v_k = current_vectors[translation_level];
    
    // Constraint: -w_k²(2 v_k·δ + ‖δ‖²)
    double linear_term = 2.0 * v_k.dot(delta);
    double quadratic_term = delta.magnitude_squared();
    
    return -w_k_squared * (linear_term + quadratic_term);
}

std::vector<double> CompensationSolver::compute_pseudoinverse_solution(
    const std::vector<double>& jacobian,
    double constraint_value,
    double ridge_lambda) {
    
    size_t n = jacobian.size(); // Number of variables (3 * num_compensation_levels)
    
    // PATCH 1: Scale normalization - find max element for conditioning
    double max_j_element = 0.0;
    for (double j : jacobian) {
        max_j_element = std::max(max_j_element, std::abs(j));
    }
    
    if (max_j_element < 1e-15) {
        throw std::runtime_error("Jacobian is effectively zero");
    }
    
    // Normalize Jacobian and constraint by max element
    std::vector<double> jacobian_normalized(n);
    for (size_t i = 0; i < n; ++i) {
        jacobian_normalized[i] = jacobian[i] / max_j_element;
    }
    double constraint_normalized = constraint_value / max_j_element;
    
    // For 1×n matrix J, pseudoinverse is: J† = J^T (JJ^T + λI)^(-1)
    // Since J is 1×n, JJ^T is 1×1 scalar
    double jjt = 0.0;
    for (double j : jacobian_normalized) {
        jjt += j * j;
    }
    
    // PATCH 2: Adaptive ridge based on condition number
    double condition_number = jjt;  // For 1x1 case, condition = value
    double adaptive_ridge = ridge_lambda;
    
    if (condition_number < 1e-10) {
        adaptive_ridge = std::max(1e-6, condition_number * 1e-3);
    } else if (condition_number > 1e10) {
        adaptive_ridge = std::max(1e-3, condition_number * 1e-12);
    }
    
    // Add adaptive ridge regularization
    double jjt_regularized = jjt + adaptive_ridge;
    
    if (std::abs(jjt_regularized) < 1e-15) {
        throw std::runtime_error("Matrix remains singular after adaptive ridge regularization");
    }
    
    // Solution: x = J^T * (constraint_normalized / jjt_regularized)
    double scale_factor = constraint_normalized / jjt_regularized;
    
    std::vector<double> solution(n);
    for (size_t i = 0; i < n; ++i) {
        // Scale back by normalization factor
        solution[i] = jacobian_normalized[i] * scale_factor;
    }
    
    return solution;
}

std::vector<Vector3> CompensationSolver::unflatten_solution(
    const std::vector<double>& solution,
    size_t num_compensation_levels) {
    
    std::vector<Vector3> compensation_deltas;
    compensation_deltas.reserve(num_compensation_levels);
    
    for (size_t i = 0; i < num_compensation_levels; ++i) {
        size_t base_idx = i * 3;
        compensation_deltas.emplace_back(
            solution[base_idx + 0],
            solution[base_idx + 1], 
            solution[base_idx + 2]
        );
    }
    
    return compensation_deltas;
}

bool CompensationSolver::validate_s_conservation(
    size_t translation_level,
    const Vector3& delta,
    const std::vector<Vector3>& compensation_deltas,
    const std::array<Vector3, NUM_LEVELS>& current_vectors,
    const std::array<double, NUM_LEVELS>& weights,
    double tolerance) {
    
    double total_delta_s = calculate_exact_delta_s(
        translation_level, delta, compensation_deltas, current_vectors, weights
    );
    
    return std::abs(total_delta_s) < tolerance;
}

double CompensationSolver::calculate_exact_delta_s(
    size_t translation_level,
    const Vector3& delta,
    const std::vector<Vector3>& compensation_deltas,
    const std::array<Vector3, NUM_LEVELS>& current_vectors,
    const std::array<double, NUM_LEVELS>& weights) {
    
    double total_delta_s = 0.0;
    
    // S-change from translation at level k
    double w_k_squared = weights[translation_level] * weights[translation_level];
    Vector3 v_k = current_vectors[translation_level];
    Vector3 v_k_new = v_k + delta;
    
    total_delta_s += w_k_squared * (v_k_new.magnitude_squared() - v_k.magnitude_squared());
    
    // S-changes from compensations at levels j > k
    size_t comp_idx = 0;
    for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
        if (comp_idx < compensation_deltas.size()) {
            double w_j_squared = weights[j] * weights[j];
            Vector3 v_j = current_vectors[j];
            Vector3 v_j_new = v_j + compensation_deltas[comp_idx];
            
            total_delta_s += w_j_squared * (v_j_new.magnitude_squared() - v_j.magnitude_squared());
            comp_idx++;
        }
    }
    
    return total_delta_s;
}

CompensationSolver::CompensationResult CompensationSolver::solve_proportional_fallback(
    size_t translation_level,
    const Vector3& delta,
    const std::array<Vector3, NUM_LEVELS>& current_vectors,
    const std::array<double, NUM_LEVELS>& weights) {
    
    CompensationResult result;
    result.used_fallback = true;
    result.iterations_used = 1;
    
    size_t num_compensation_levels = NUM_LEVELS - translation_level - 1;
    if (num_compensation_levels == 0) {
        result.converged = true;
        return result;
    }
    
    // Calculate exact S-change from translation
    double w_k_squared = weights[translation_level] * weights[translation_level];
    Vector3 v_k = current_vectors[translation_level];
    Vector3 v_k_new = v_k + delta;
    double delta_s_exact = w_k_squared * (v_k_new.magnitude_squared() - v_k.magnitude_squared());
    
    // Distribute proportionally based on weights
    double total_weight_squared = 0.0;
    for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
        total_weight_squared += weights[j] * weights[j];
    }
    
    if (total_weight_squared < 1e-15) {
        result.error_message = "All compensation weights are effectively zero";
        return result;
    }
    
    double target_delta_s_total = -delta_s_exact;
    
    for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
        double w_j_squared = weights[j] * weights[j];
        double proportion = w_j_squared / total_weight_squared;
        double target_delta_s_j = target_delta_s_total * proportion;
        
        Vector3 v_j = current_vectors[j];
        double current_magnitude_squared = v_j.magnitude_squared();
        double required_magnitude_squared = current_magnitude_squared + target_delta_s_j / w_j_squared;
        
        if (required_magnitude_squared >= 0.0) {
            double required_magnitude = std::sqrt(required_magnitude_squared);
            double current_magnitude = v_j.magnitude();
            
            Vector3 compensation_delta;
            if (current_magnitude > 1e-12) {
                Vector3 scaled_v_j = v_j.normalized() * required_magnitude;
                compensation_delta = scaled_v_j - v_j;
            } else {
                // Use perpendicular direction for zero vectors
                Vector3 direction = delta.normalized();
                Vector3 perp_direction;
                if (std::abs(direction.x()) < 0.9) {
                    perp_direction = Vector3(1.0, 0.0, 0.0).cross(direction).normalized();
                } else {
                    perp_direction = Vector3(0.0, 1.0, 0.0).cross(direction).normalized();
                }
                compensation_delta = perp_direction * required_magnitude;
            }
            
            result.compensation_deltas.push_back(compensation_delta);
        } else {
            result.compensation_deltas.push_back(Vector3(0.0, 0.0, 0.0));
        }
    }
    
    // Validate result
    result.converged = validate_s_conservation(
        translation_level, delta, result.compensation_deltas, 
        current_vectors, weights, DEFAULT_CONVERGENCE_TOLERANCE
    );
    
    if (!result.converged) {
        result.residual_error = calculate_exact_delta_s(
            translation_level, delta, result.compensation_deltas, current_vectors, weights
        );
    }
    
    return result;
}

CompensationSolver::CompensationResult CompensationSolver::refine_solution_iteratively(
    const CompensationResult& initial_solution,
    size_t translation_level,
    const Vector3& delta,
    const std::array<Vector3, NUM_LEVELS>& current_vectors,
    const std::array<double, NUM_LEVELS>& weights,
    const SolverConfig& config) {
    
    CompensationResult result = initial_solution;
    
    for (int iteration = 1; iteration < config.max_iterations; ++iteration) {
        // Calculate current residual error
        double current_error = calculate_exact_delta_s(
            translation_level, delta, result.compensation_deltas, current_vectors, weights
        );
        
        if (std::abs(current_error) < config.convergence_tolerance) {
            result.converged = true;
            result.iterations_used = iteration;
            result.residual_error = current_error;
            break;
        }
        
        // Build correction constraint to eliminate residual error
        std::vector<double> jacobian = build_jacobian_matrix(translation_level, current_vectors, weights);
        double correction_constraint = -current_error;
        
        try {
            std::vector<double> correction_solution = compute_pseudoinverse_solution(
                jacobian, correction_constraint, config.ridge_lambda
            );
            
            std::vector<Vector3> correction_deltas = unflatten_solution(
                correction_solution, NUM_LEVELS - translation_level - 1
            );
            
            // Apply corrections to current solution
            for (size_t i = 0; i < result.compensation_deltas.size() && i < correction_deltas.size(); ++i) {
                result.compensation_deltas[i] = result.compensation_deltas[i] + correction_deltas[i];
            }
            
        } catch (const std::exception& e) {
            result.error_message = "Iterative refinement failed: " + std::string(e.what());
            break;
        }
        
        result.iterations_used = iteration + 1;
    }
    
    // Final validation
    result.residual_error = calculate_exact_delta_s(
        translation_level, delta, result.compensation_deltas, current_vectors, weights
    );
    
    if (!result.converged) {
        result.converged = (std::abs(result.residual_error) < config.convergence_tolerance);
    }
    
    return result;
}

// Helper function for error reporting
double calculate_current_S(const std::array<Vector3, CompensationSolver::NUM_LEVELS>& vectors,
                          const std::array<double, CompensationSolver::NUM_LEVELS>& weights) {
    double S = 0.0;
    for (size_t i = 0; i < CompensationSolver::NUM_LEVELS; ++i) {
        S += weights[i] * weights[i] * vectors[i].magnitude_squared();
    }
    return S;
}

} // namespace hsml::core