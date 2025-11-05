#include "hsml/core/production_compensation_solver.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <sstream>

namespace hsml::core {

// ========== Primary Trust-Region Solver ==========

ProductionCompensationSolver::SolveResult ProductionCompensationSolver::solve_translate(
    size_t translation_level,
    const Vector3& delta, 
    const std::array<Vector3, NUM_LEVELS>& v,
    const std::array<double, NUM_LEVELS>& w,
    double S_before,
    double epsilon_absolute,
    int max_outer_iterations) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    SolveResult result;
    result.status = SolveStatus::CATASTROPHIC_FAILURE;
    result.successful_method = SolverMethod::TRUST_REGION_PRIMARY;
    result.used_fallback = false;
    result.halving_count = 0;
    result.condition_number = 0.0;
    
    // Validate inputs
    if (translation_level >= NUM_LEVELS) {
        result.diagnostics = "Invalid translation level";
        return result;
    }
    
    // Determine compensation levels (j > k)
    std::vector<size_t> compensation_indices;
    for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
        compensation_indices.push_back(j);
    }
    
    if (compensation_indices.empty()) {
        // No compensation needed - trivial case
        result.status = SolveStatus::SUCCESS;
        result.final_residual_error = 0.0;
        result.iterations_used = 0;
        return result;
    }
    
    // Calculate primary constraint terms
    const double wk_squared = w[translation_level] * w[translation_level];
    const Vector3 vk = v[translation_level];
    const double linear_primary = wk_squared * 2.0 * vk.dot(delta);
    const double quadratic_primary = wk_squared * delta.magnitude_squared();
    
    // Initialize trust-region parameters
    double step_scale = 1.0;
    double best_residual = std::numeric_limits<double>::infinity();
    std::vector<Vector3> best_deltas(NUM_LEVELS, Vector3(0.0, 0.0, 0.0));
    
    // Main trust-region iteration loop
    for (int outer = 0; outer < max_outer_iterations; ++outer) {
        result.iterations_used = outer + 1;
        
        try {
            // Build scaled Jacobian for current step scale
            auto [jacobian_normalized, scaling_factor] = build_scaled_jacobian(
                translation_level, step_scale, v, w);
            
            if (jacobian_normalized.empty() || scaling_factor < 1e-15) {
                result.diagnostics = "Jacobian is effectively zero";
                result.status = SolveStatus::SINGULAR_JACOBIAN;
                break;
            }
            
            // Calculate condition number
            result.condition_number = estimate_condition_number(jacobian_normalized);
            
            // RHS: want total ΔS = -(linear_primary*step + quad_primary*step²)
            double target_delta_S = -(linear_primary * step_scale + 
                                     quadratic_primary * step_scale * step_scale);
            double constraint_rhs = target_delta_S / scaling_factor;
            
            // Calculate adaptive ridge parameter
            double ridge_param = calculate_adaptive_ridge(jacobian_normalized);
            
            // Solve Moore-Penrose pseudoinverse
            std::vector<double> solution = solve_moore_penrose(
                jacobian_normalized, constraint_rhs, ridge_param);
            
            // Convert back to Vector3 deltas
            std::vector<Vector3> compensation_deltas = unflatten_compensation_deltas(
                solution, translation_level);
            
            // Evaluate S with proposed deltas
            double S_after = evaluate_S_with_deltas(
                translation_level, delta * step_scale, compensation_deltas, v, w);
            double residual = std::abs(S_after - S_before);
            
            result.residual_history.push_back(residual);
            
            // Check for improvement
            if (residual < best_residual) {
                best_residual = residual;
                best_deltas = compensation_deltas;
            }
            
            // Check convergence
            if (residual < epsilon_absolute) {
                result.status = SolveStatus::SUCCESS;
                result.final_residual_error = residual;
                result.compensation_deltas = compensation_deltas;
                result.step_scale_final = step_scale;
                break;
            }
            
            // Trust-region backtracking - halve step size
            step_scale *= STEP_HALVING_FACTOR;
            result.halving_count++;
            
            if (step_scale < MIN_STEP_SCALE) {
                result.status = SolveStatus::STEP_SIZE_UNDERFLOW;
                result.diagnostics = "Step scale underflow: " + std::to_string(step_scale);
                break;
            }
            
        } catch (const std::exception& e) {
            result.diagnostics = "Exception in trust-region: " + std::string(e.what());
            result.status = SolveStatus::CATASTROPHIC_FAILURE;
            break;
        }
    }
    
    // If primary method failed, try fallback methods
    if (result.status != SolveStatus::SUCCESS) {
        result.used_fallback = true;
        
        // Try constrained least-squares fallback
        SolveResult fallback_result = solve_constrained_least_squares(
            translation_level, delta, v, w, S_before, epsilon_absolute);
        
        if (fallback_result.status == SolveStatus::SUCCESS) {
            result = fallback_result;
            result.status = SolveStatus::FELL_BACK_TO_SECONDARY;
            result.successful_method = SolverMethod::CONSTRAINED_LEAST_SQUARES;
            result.used_fallback = true;
        } else {
            // Emergency proportional fallback
            SolveResult emergency_result = solve_proportional_emergency(
                translation_level, delta, v, w, S_before, epsilon_absolute * 100.0);
            
            if (emergency_result.status == SolveStatus::SUCCESS) {
                result = emergency_result;
                result.status = SolveStatus::FELL_BACK_TO_EMERGENCY;
                result.successful_method = SolverMethod::PROPORTIONAL_FALLBACK;
                result.used_fallback = true;
            } else {
                // All methods failed - use best attempt
                result.compensation_deltas = best_deltas;
                result.final_residual_error = best_residual;
                result.diagnostics += " [All fallbacks failed, using best attempt]";
            }
        }
    }
    
    // Calculate final timing
    auto end_time = std::chrono::high_resolution_clock::now();
    result.solve_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();
    
    return result;
}

// ========== Constrained Least-Squares Fallback ==========

ProductionCompensationSolver::SolveResult ProductionCompensationSolver::solve_constrained_least_squares(
    size_t translation_level,
    const Vector3& delta,
    const std::array<Vector3, NUM_LEVELS>& v,
    const std::array<double, NUM_LEVELS>& w,
    double S_before,
    double epsilon_absolute) {
    
    SolveResult result;
    result.successful_method = SolverMethod::CONSTRAINED_LEAST_SQUARES;
    result.status = SolveStatus::CATASTROPHIC_FAILURE;
    result.iterations_used = 1;
    
    try {
        // Build Jacobian without step scaling (use full delta)
        auto [jacobian, scaling_factor] = build_scaled_jacobian(translation_level, 1.0, v, w);
        
        if (jacobian.empty() || scaling_factor < 1e-15) {
            result.diagnostics = "Singular Jacobian in constrained LS";
            result.status = SolveStatus::SINGULAR_JACOBIAN;
            return result;
        }
        
        // Calculate target constraint
        const double wk_squared = w[translation_level] * w[translation_level];
        const Vector3 vk = v[translation_level];
        double target_delta_S = -wk_squared * (2.0 * vk.dot(delta) + delta.magnitude_squared());
        double constraint_rhs = target_delta_S / scaling_factor;
        
        // Lagrange multiplier method: min ||x||² subject to J·x = b
        // Solution: x = J^T * λ where λ = (JJ^T)^-1 * b
        double jjt = std::inner_product(jacobian.begin(), jacobian.end(), jacobian.begin(), 0.0);
        
        if (std::abs(jjt) < 1e-15) {
            result.diagnostics = "JJ^T singular in constrained LS";
            result.status = SolveStatus::SINGULAR_JACOBIAN;
            return result;
        }
        
        double lambda = constraint_rhs / jjt;
        std::vector<double> solution(jacobian.size());
        for (size_t i = 0; i < jacobian.size(); ++i) {
            solution[i] = jacobian[i] * lambda;
        }
        
        // Convert to compensation deltas
        result.compensation_deltas = unflatten_compensation_deltas(solution, translation_level);
        
        // Validate result
        double S_after = evaluate_S_with_deltas(translation_level, delta, result.compensation_deltas, v, w);
        result.final_residual_error = std::abs(S_after - S_before);
        
        if (result.final_residual_error < epsilon_absolute) {
            result.status = SolveStatus::SUCCESS;
            result.diagnostics = "Constrained LS succeeded";
        } else {
            result.status = SolveStatus::CATASTROPHIC_FAILURE;
            result.diagnostics = "Constrained LS residual too large: " + std::to_string(result.final_residual_error);
        }
        
    } catch (const std::exception& e) {
        result.diagnostics = "Exception in constrained LS: " + std::string(e.what());
        result.status = SolveStatus::CATASTROPHIC_FAILURE;
    }
    
    return result;
}

// ========== Emergency Proportional Fallback ==========

ProductionCompensationSolver::SolveResult ProductionCompensationSolver::solve_proportional_emergency(
    size_t translation_level,
    const Vector3& delta,
    const std::array<Vector3, NUM_LEVELS>& v,
    const std::array<double, NUM_LEVELS>& w,
    double S_before,
    double epsilon_absolute) {
    
    SolveResult result;
    result.successful_method = SolverMethod::PROPORTIONAL_FALLBACK;
    result.status = SolveStatus::SUCCESS;  // Emergency mode - be optimistic
    result.iterations_used = 1;
    
    try {
        // Calculate exact S-change from translation
        const double wk_squared = w[translation_level] * w[translation_level];
        const Vector3 vk = v[translation_level];
        const Vector3 vk_new = vk + delta;
        double delta_s_exact = wk_squared * (vk_new.magnitude_squared() - vk.magnitude_squared());
        
        // Calculate total weight for proportional distribution
        double total_weight_squared = 0.0;
        for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
            total_weight_squared += w[j] * w[j];
        }
        
        if (total_weight_squared < 1e-15) {
            result.diagnostics = "All compensation weights zero in emergency mode";
            result.status = SolveStatus::CATASTROPHIC_FAILURE;
            return result;
        }
        
        double target_delta_s_total = -delta_s_exact;
        result.compensation_deltas.assign(NUM_LEVELS, Vector3(0.0, 0.0, 0.0));
        
        // Distribute proportionally by weight squared
        for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
            double wj_squared = w[j] * w[j];
            double proportion = wj_squared / total_weight_squared;
            double target_delta_s_j = target_delta_s_total * proportion;
            
            Vector3 vj = v[j];
            double current_mag_sq = vj.magnitude_squared();
            double required_mag_sq = current_mag_sq + target_delta_s_j / wj_squared;
            
            if (required_mag_sq >= 0.0) {
                double required_mag = std::sqrt(required_mag_sq);
                double current_mag = vj.magnitude();
                
                if (current_mag > 1e-12) {
                    Vector3 scaled_vj = vj.normalized() * required_mag;
                    result.compensation_deltas[j] = scaled_vj - vj;
                } else {
                    // Use delta direction for zero vectors
                    Vector3 direction = (delta.magnitude() > 1e-12) ? 
                                      delta.normalized() : Vector3(1.0, 0.0, 0.0);
                    result.compensation_deltas[j] = direction * required_mag;
                }
            }
        }
        
        // Validate emergency result
        double S_after = evaluate_S_with_deltas(translation_level, delta, result.compensation_deltas, v, w);
        result.final_residual_error = std::abs(S_after - S_before);
        
        if (result.final_residual_error < epsilon_absolute) {
            result.diagnostics = "Emergency proportional succeeded";
        } else {
            result.diagnostics = "Emergency proportional - residual " + 
                               std::to_string(result.final_residual_error) + " exceeds tolerance";
        }
        
    } catch (const std::exception& e) {
        result.diagnostics = "Exception in emergency proportional: " + std::string(e.what());
        result.status = SolveStatus::CATASTROPHIC_FAILURE;
    }
    
    return result;
}

// ========== Utility Functions ==========

double ProductionCompensationSolver::calculate_deterministic_tolerance(
    const std::array<Vector3, NUM_LEVELS>& v,
    const std::array<double, NUM_LEVELS>& w) {
    
    double S_total = calculate_S_total(v, w);
    return std::max(1e-12 * S_total, 1e-18);
}

double ProductionCompensationSolver::evaluate_S_with_deltas(
    size_t translation_level,
    const Vector3& delta,
    const std::vector<Vector3>& compensation_deltas,
    const std::array<Vector3, NUM_LEVELS>& v,
    const std::array<double, NUM_LEVELS>& w) {
    
    double S_total = 0.0;
    
    // Apply translation at level k
    Vector3 vk_new = v[translation_level] + delta;
    S_total += w[translation_level] * w[translation_level] * vk_new.magnitude_squared();
    
    // Apply compensations at levels j > k
    for (size_t j = 0; j < NUM_LEVELS; ++j) {
        if (j == translation_level) continue;
        
        Vector3 vj_new = v[j];
        if (j < compensation_deltas.size()) {
            vj_new = vj_new + compensation_deltas[j];
        }
        S_total += w[j] * w[j] * vj_new.magnitude_squared();
    }
    
    return S_total;
}

double ProductionCompensationSolver::calculate_S_total(
    const std::array<Vector3, NUM_LEVELS>& v,
    const std::array<double, NUM_LEVELS>& w) {
    
    double S_total = 0.0;
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        S_total += w[i] * w[i] * v[i].magnitude_squared();
    }
    return S_total;
}

// ========== Post-Commit Verification and Rollback ==========

ProductionCompensationSolver::SolveResult ProductionCompensationSolver::translate_with_verification(
    size_t translation_level,
    const Vector3& delta,
    std::array<Vector3, NUM_LEVELS>& v,
    const std::array<double, NUM_LEVELS>& w,
    int max_outer_iterations) {
    
    // Calculate S-invariant before modification
    double S_before = calculate_S_total(v, w);
    
    // Calculate deterministic tolerance
    double epsilon_absolute = calculate_deterministic_tolerance(v, w);
    
    // Create backup of current state for rollback
    std::array<Vector3, NUM_LEVELS> v_backup = v;
    
    // Solve for compensation deltas
    SolveResult result = solve_translate(
        translation_level, delta, v, w, S_before, epsilon_absolute, max_outer_iterations);
    
    if (result.status != SolveStatus::SUCCESS && 
        result.status != SolveStatus::FELL_BACK_TO_SECONDARY &&
        result.status != SolveStatus::FELL_BACK_TO_EMERGENCY) {
        
        // Solver failed completely - no changes to apply
        result.diagnostics += " [No state changes applied due to solver failure]";
        return result;
    }
    
    try {
        // Apply translation at level k
        v[translation_level] = v[translation_level] + delta;
        
        // Apply compensation deltas at levels j > k
        for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
            if (j < result.compensation_deltas.size()) {
                v[j] = v[j] + result.compensation_deltas[j];
            }
        }
        
        // Post-commit verification: Check actual S-conservation
        double S_after = calculate_S_total(v, w);
        double post_commit_residual = std::abs(S_after - S_before);
        
        // Maximum error negation: If post-commit drift exceeds tolerance, rollback
        if (post_commit_residual > epsilon_absolute) {
            // ROLLBACK: Restore original state
            v = v_backup;
            
            result.status = SolveStatus::POST_COMMIT_DRIFT;
            result.final_residual_error = post_commit_residual;
            result.diagnostics += " [ROLLBACK: Post-commit S-drift " + 
                                std::to_string(post_commit_residual) + 
                                " > tolerance " + std::to_string(epsilon_absolute) + "]";
            
            return result;
        }
        
        // Success: Post-commit verification passed
        result.final_residual_error = post_commit_residual;
        result.diagnostics += " [Post-commit verified: ΔS=" + std::to_string(post_commit_residual) + "]";
        
        // Ensure status reflects actual outcome
        if (result.status == SolveStatus::FELL_BACK_TO_EMERGENCY && post_commit_residual < epsilon_absolute) {
            // Emergency method succeeded after verification
            result.diagnostics += " [Emergency method validated by post-commit check]";
        }
        
    } catch (const std::exception& e) {
        // Exception during apply - rollback
        v = v_backup;
        result.status = SolveStatus::CATASTROPHIC_FAILURE;
        result.diagnostics = "Exception during apply phase: " + std::string(e.what()) + " [ROLLBACK performed]";
    }
    
    return result;
}

// ========== Private Helper Functions ==========

std::pair<std::vector<double>, double> ProductionCompensationSolver::build_scaled_jacobian(
    size_t translation_level,
    double step_scale,
    const std::array<Vector3, NUM_LEVELS>& v,
    const std::array<double, NUM_LEVELS>& w) {
    
    std::vector<double> jacobian;
    
    // Build Jacobian for compensation levels j > k
    for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
        double wj_squared = w[j] * w[j];
        Vector3 vj = v[j];
        
        // Jacobian entries: 2 * w_j² * v_j components
        jacobian.push_back(2.0 * wj_squared * vj.x());
        jacobian.push_back(2.0 * wj_squared * vj.y());
        jacobian.push_back(2.0 * wj_squared * vj.z());
    }
    
    if (jacobian.empty()) {
        return {jacobian, 0.0};
    }
    
    // Find max element for scaling
    double max_element = 0.0;
    for (double j : jacobian) {
        max_element = std::max(max_element, std::abs(j));
    }
    
    if (max_element < 1e-15) {
        return {jacobian, 0.0};
    }
    
    // Scale Jacobian
    std::vector<double> jacobian_normalized(jacobian.size());
    for (size_t i = 0; i < jacobian.size(); ++i) {
        jacobian_normalized[i] = jacobian[i] / max_element;
    }
    
    return {jacobian_normalized, max_element};
}

double ProductionCompensationSolver::calculate_adaptive_ridge(
    const std::vector<double>& jacobian_normalized) {
    
    // Calculate JJ^T for ridge estimation
    double jjt = std::inner_product(jacobian_normalized.begin(), jacobian_normalized.end(),
                                   jacobian_normalized.begin(), 0.0);
    
    // Adaptive ridge based on condition number proxy
    if (jjt < 1e-10) {
        return std::max(1e-6, jjt * 1e-3);
    } else if (jjt > 1e10) {
        return std::max(1e-3, jjt * 1e-12);
    } else {
        return std::max(1e-8, jjt * 1e-6);
    }
}

std::vector<double> ProductionCompensationSolver::solve_moore_penrose(
    const std::vector<double>& jacobian_normalized,
    double constraint_rhs,
    double ridge_parameter) {
    
    // Moore-Penrose: x = J^T * (JJ^T + ridge)^-1 * b
    double jjt = std::inner_product(jacobian_normalized.begin(), jacobian_normalized.end(),
                                   jacobian_normalized.begin(), 0.0);
    double jjt_regularized = jjt + ridge_parameter;
    
    if (std::abs(jjt_regularized) < 1e-15) {
        throw std::runtime_error("Singular regularized system in Moore-Penrose");
    }
    
    double scale_factor = constraint_rhs / jjt_regularized;
    
    std::vector<double> solution(jacobian_normalized.size());
    for (size_t i = 0; i < jacobian_normalized.size(); ++i) {
        solution[i] = jacobian_normalized[i] * scale_factor;
    }
    
    return solution;
}

std::vector<Vector3> ProductionCompensationSolver::unflatten_compensation_deltas(
    const std::vector<double>& solution,
    size_t translation_level) {
    
    std::vector<Vector3> compensation_deltas(NUM_LEVELS, Vector3(0.0, 0.0, 0.0));
    
    size_t sol_idx = 0;
    for (size_t j = translation_level + 1; j < NUM_LEVELS && sol_idx + 2 < solution.size(); ++j) {
        compensation_deltas[j] = Vector3(
            solution[sol_idx + 0],
            solution[sol_idx + 1],
            solution[sol_idx + 2]
        );
        sol_idx += 3;
    }
    
    return compensation_deltas;
}

double ProductionCompensationSolver::estimate_condition_number(
    const std::vector<double>& jacobian) {
    
    if (jacobian.empty()) return 1.0;
    
    double max_element = *std::max_element(jacobian.begin(), jacobian.end(),
                                         [](double a, double b) { return std::abs(a) < std::abs(b); });
    double min_element = *std::min_element(jacobian.begin(), jacobian.end(),
                                         [](double a, double b) { return std::abs(a) > std::abs(b); });
    
    if (std::abs(min_element) < 1e-15) {
        return 1e15;  // Effectively singular
    }
    
    return std::abs(max_element / min_element);
}

} // namespace hsml::core