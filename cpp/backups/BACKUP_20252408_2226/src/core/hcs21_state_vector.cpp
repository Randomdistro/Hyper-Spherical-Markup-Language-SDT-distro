// ⚠️ CRITICAL WARNING: HCS21 STATE VECTOR CONTAMINATION ⚠️
// This file contains HERETICAL use of Vector3 for 21D state management!
// 21-dimensional hierarchical cosmological states should be PURE SPHERICAL!
// Vector3 xyz coordinates are FORBIDDEN in cosmological state space!

#ifndef HSML_ALLOW_CARTESIAN_HERESY
#error "HCS21StateVector uses Vector3 xyz - BANISHED! Reimplement with pure spherical state representation. Define HSML_ALLOW_CARTESIAN_HERESY to override."
#endif

#include "hsml/core/hcs21_state_vector.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace hsml::core {

// ========== Constructors ==========

HCS21StateVector::HCS21StateVector() 
    : has_uncertainties_(false)
{
    // Initialize with zero vectors
    level_vectors_.fill(Vector3(0.0, 0.0, 0.0));
    
    // Initialize with default scale weights
    scale_weights_ = create_default_weights();
    
    // Initialize uncertainties to zero
    uncertainties_.fill(0.0);
}

HCS21StateVector::HCS21StateVector(const std::array<Vector3, NUM_LEVELS>& level_vectors,
                                   const std::array<double, NUM_LEVELS>& scale_weights,
                                   const std::array<double, NUM_LEVELS>& uncertainties)
    : level_vectors_(level_vectors)
    , scale_weights_(scale_weights)
    , uncertainties_(uncertainties)
    , has_uncertainties_(true)
{
    // Validate weights are strictly decreasing and positive
    if (!validate_weights(scale_weights_)) {
        throw std::invalid_argument("Scale weights must be positive and strictly decreasing");
    }
    
    // If no uncertainties provided, set has_uncertainties_ to false
    if (std::all_of(uncertainties.begin(), uncertainties.end(), [](double u) { return u == 0.0; })) {
        has_uncertainties_ = false;
    }
}

// ========== Core Properties ==========

Vector3& HCS21StateVector::get_level_vector(size_t level) {
    validate_level_index(level);
    return level_vectors_[level];
}

const Vector3& HCS21StateVector::get_level_vector(size_t level) const {
    validate_level_index(level);
    return level_vectors_[level];
}

Vector3& HCS21StateVector::get_level_vector(Level level) {
    return get_level_vector(static_cast<size_t>(level));
}

const Vector3& HCS21StateVector::get_level_vector(Level level) const {
    return get_level_vector(static_cast<size_t>(level));
}

double HCS21StateVector::get_weight(size_t level) const {
    validate_level_index(level);
    return scale_weights_[level];
}

double HCS21StateVector::get_weight(Level level) const {
    return get_weight(static_cast<size_t>(level));
}

const std::array<Vector3, HCS21StateVector::NUM_LEVELS>& HCS21StateVector::get_all_vectors() const {
    return level_vectors_;
}

const std::array<double, HCS21StateVector::NUM_LEVELS>& HCS21StateVector::get_all_weights() const {
    return scale_weights_;
}

// ========== S-Invariant Calculations ==========

double HCS21StateVector::calculate_S_invariant() const {
    double S = 0.0;
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        double weight_squared = scale_weights_[i] * scale_weights_[i];
        double potential = displacement_potential(level_vectors_[i]);
        S += weight_squared * potential;
    }
    return S;
}

double HCS21StateVector::displacement_potential(const Vector3& vector) {
    // v1.0 specification: Φ(v) = ‖v‖²
    return vector.magnitude_squared();
}

bool HCS21StateVector::validate_S_conservation(const HCS21StateVector& other, 
                                               double tolerance) const {
    double S_this = calculate_S_invariant();
    double S_other = other.calculate_S_invariant();
    return std::abs(S_this - S_other) < tolerance;
}

// ========== Legal Transform Operations ==========

void HCS21StateVector::rotate_at_level(size_t level, const Matrix4& rotation, 
                                       double tolerance) {
    validate_level_index(level);
    
    // Validate rotation matrix
    if (!validate_rotation_matrix(rotation, tolerance)) {
        throw std::invalid_argument("Matrix is not a valid rotation (not in SO(3))");
    }
    
    // Store original S-invariant for validation
    double original_S = calculate_S_invariant();
    
    // Apply rotation cascade
    apply_rotation_cascade(level, rotation);
    
    // Validate S-conservation
    double new_S = calculate_S_invariant();
    if (std::abs(new_S - original_S) > tolerance) {
        throw std::runtime_error("S-invariant conservation violated during rotation");
    }
}

void HCS21StateVector::rotate_at_level(Level level, const Matrix4& rotation, 
                                       double tolerance) {
    rotate_at_level(static_cast<size_t>(level), rotation, tolerance);
}

void HCS21StateVector::translate_at_level(size_t level, const Vector3& delta, 
                                         double tolerance) {
    validate_level_index(level);
    
    // Store original S-invariant for validation
    double original_S = calculate_S_invariant();
    
    // Use Moore-Penrose compensation solver for exact S-conservation
    CompensationSolver::SolverConfig config;
    config.convergence_tolerance = tolerance;
    config.use_iterative_refinement = true;
    config.use_proportional_fallback = true;
    
    auto compensation_result = CompensationSolver::solve_compensation(
        level, delta, level_vectors_, scale_weights_, config
    );
    
    if (!compensation_result.converged) {
        throw std::runtime_error("Translation compensation failed: " + 
                               compensation_result.error_message + 
                               " (Residual error: " + std::to_string(compensation_result.residual_error) + ")");
    }
    
    // Apply translation to the specified level
    level_vectors_[level] = level_vectors_[level] + delta;
    
    // Apply compensation to lower levels
    size_t comp_index = 0;
    for (size_t i = level + 1; i < NUM_LEVELS; ++i) {
        if (comp_index < compensation_result.compensation_deltas.size()) {
            level_vectors_[i] = level_vectors_[i] + compensation_result.compensation_deltas[comp_index];
            comp_index++;
        }
    }
    
    // Final validation of S-conservation
    double new_S = calculate_S_invariant();
    if (std::abs(new_S - original_S) > tolerance) {
        throw std::runtime_error("S-invariant conservation violated after Moore-Penrose compensation. " +
                               std::string("Original S: ") + std::to_string(original_S) + 
                               ", New S: " + std::to_string(new_S) + 
                               ", Difference: " + std::to_string(std::abs(new_S - original_S)) +
                               ", Solver iterations: " + std::to_string(compensation_result.iterations_used));
    }
}

void HCS21StateVector::translate_at_level(Level level, const Vector3& delta, 
                                         double tolerance) {
    translate_at_level(static_cast<size_t>(level), delta, tolerance);
}

// ========== Projection Methods to 3D ==========

Vector3 HCS21StateVector::project_weighted() const {
    Vector3 result(0.0, 0.0, 0.0);
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        result = result + (level_vectors_[i] * scale_weights_[i]);
    }
    return result;
}

Vector3 HCS21StateVector::project_dominant(double eps_obs) const {
    // Find finest scale exceeding observation threshold
    for (int i = NUM_LEVELS - 1; i >= 0; --i) {
        double visibility = scale_weights_[i] * level_vectors_[i].magnitude();
        if (visibility > eps_obs) {
            return level_vectors_[i];
        }
    }
    
    // If no level exceeds threshold, return coarsest level
    return level_vectors_[0];
}

Vector3 HCS21StateVector::project_dominant_smooth(double eps_obs, double width) const {
    std::array<double, NUM_LEVELS> alpha;
    double alpha_sum = 0.0;
    
    // Calculate visibility-based weights using tanh transitions
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        double visibility = scale_weights_[i] * level_vectors_[i].magnitude();
        double tanh_arg = (visibility - eps_obs) / (width * eps_obs);
        alpha[i] = 0.5 * (1.0 + std::tanh(tanh_arg));
        alpha_sum += alpha[i];
    }
    
    // Normalize alpha weights
    if (alpha_sum > 1e-12) {
        for (size_t i = 0; i < NUM_LEVELS; ++i) {
            alpha[i] /= alpha_sum;
        }
    } else {
        // If all alphas are near zero, use uniform weighting
        alpha.fill(1.0 / NUM_LEVELS);
    }
    
    // Calculate weighted sum
    Vector3 result(0.0, 0.0, 0.0);
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        result = result + (level_vectors_[i] * alpha[i]);
    }
    
    return result;
}

// ========== Uncertainty Propagation ==========

void HCS21StateVector::set_uncertainties(const std::array<double, NUM_LEVELS>& uncertainties) {
    uncertainties_ = uncertainties;
    has_uncertainties_ = true;
}

double HCS21StateVector::get_uncertainty(size_t level) const {
    validate_level_index(level);
    return uncertainties_[level];
}

double HCS21StateVector::propagate_uncertainty() const {
    if (!has_uncertainties_) {
        return 0.0;
    }
    
    double sigma_S_squared = 0.0;
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        double w4 = std::pow(scale_weights_[i], 4);
        double norm_squared = level_vectors_[i].magnitude_squared();
        double sigma_squared = uncertainties_[i] * uncertainties_[i];
        
        sigma_S_squared += w4 * norm_squared * sigma_squared;
    }
    
    return 2.0 * std::sqrt(sigma_S_squared);
}

// ========== Validation & Utilities ==========

bool HCS21StateVector::validate_weights(const std::array<double, NUM_LEVELS>& weights) {
    // Check all weights are positive
    for (double w : weights) {
        if (w <= 0.0) {
            return false;
        }
    }
    
    // Check strictly decreasing property
    for (size_t i = 0; i < NUM_LEVELS - 1; ++i) {
        if (weights[i] <= weights[i + 1]) {
            return false;
        }
    }
    
    return true;
}

bool HCS21StateVector::validate_rotation_matrix(const Matrix4& matrix, double tolerance) {
    // Extract 3x3 rotation part from 4x4 matrix
    // Matrix4 stores data in column-major order: [0,1,2,3, 4,5,6,7, 8,9,10,11, 12,13,14,15]
    // Extract upper-left 3x3 submatrix
    
    double R[9]; // 3x3 matrix in row-major order
    const double* m = matrix.data().data();
    
    // Extract 3x3 rotation part (assuming Matrix4 is column-major)
    R[0] = m[0];  R[1] = m[4];  R[2] = m[8];   // First row
    R[3] = m[1];  R[4] = m[5];  R[5] = m[9];   // Second row  
    R[6] = m[2];  R[7] = m[6];  R[8] = m[10];  // Third row
    
    // Check 1: Determinant should be +1 (proper rotation, not reflection)
    double det = R[0] * (R[4] * R[8] - R[5] * R[7]) -
                 R[1] * (R[3] * R[8] - R[5] * R[6]) +
                 R[2] * (R[3] * R[7] - R[4] * R[6]);
    
    if (std::abs(det - 1.0) > tolerance) {
        return false;
    }
    
    // Check 2: Orthogonality - R^T * R should be identity
    // Compute R^T * R and check if it's close to identity
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double dot_product = 0.0;
            for (int k = 0; k < 3; ++k) {
                dot_product += R[k * 3 + i] * R[k * 3 + j]; // R^T[i,k] * R[k,j]
            }
            
            double expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(dot_product - expected) > tolerance) {
                return false;
            }
        }
    }
    
    // Check 3: Column vectors should have unit length
    for (int j = 0; j < 3; ++j) {
        double column_length_squared = 0.0;
        for (int i = 0; i < 3; ++i) {
            double element = R[i * 3 + j];
            column_length_squared += element * element;
        }
        
        if (std::abs(column_length_squared - 1.0) > tolerance) {
            return false;
        }
    }
    
    return true;
}

std::array<double, HCS21StateVector::NUM_LEVELS> HCS21StateVector::create_default_weights(
    double k_galaxy, double k_stellar, double k_planetary, 
    double k_regional, double k_local, double k_atomic) {
    
    std::array<double, NUM_LEVELS> weights;
    
    // Start with speed of light (level 0)
    weights[0] = SPEED_OF_LIGHT;  // ~3e8 m/s
    
    // Apply k-value ratios to create strictly decreasing sequence
    // Each k-value should be < 1.0 to ensure decreasing weights
    weights[1] = weights[0] * k_galaxy;      // ~3e5 m/s (300 km/s)
    weights[2] = weights[1] * k_stellar;     // Stellar escape velocity scale
    weights[3] = weights[2] * k_planetary;   // Planetary escape velocity scale  
    weights[4] = weights[3] * k_regional;    // Regional velocity scale
    weights[5] = weights[4] * k_local;       // Local velocity scale
    weights[6] = weights[5] * k_atomic;      // Atomic velocity scale
    
    return weights;
}

// ========== Operators ==========

HCS21StateVector HCS21StateVector::operator+(const HCS21StateVector& other) const {
    std::array<Vector3, NUM_LEVELS> result_vectors;
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        result_vectors[i] = level_vectors_[i] + other.level_vectors_[i];
    }
    
    return HCS21StateVector(result_vectors, scale_weights_, uncertainties_);
}

HCS21StateVector HCS21StateVector::operator-(const HCS21StateVector& other) const {
    std::array<Vector3, NUM_LEVELS> result_vectors;
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        result_vectors[i] = level_vectors_[i] - other.level_vectors_[i];
    }
    
    return HCS21StateVector(result_vectors, scale_weights_, uncertainties_);
}

HCS21StateVector HCS21StateVector::operator*(double scalar) const {
    std::array<Vector3, NUM_LEVELS> result_vectors;
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        result_vectors[i] = level_vectors_[i] * scalar;
    }
    
    return HCS21StateVector(result_vectors, scale_weights_, uncertainties_);
}

bool HCS21StateVector::operator==(const HCS21StateVector& other) const {
    const double tolerance = 1e-10;
    
    // Compare weights
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        if (std::abs(scale_weights_[i] - other.scale_weights_[i]) > tolerance) {
            return false;
        }
    }
    
    // Compare vectors
    for (size_t i = 0; i < NUM_LEVELS; ++i) {
        if ((level_vectors_[i] - other.level_vectors_[i]).magnitude() > tolerance) {
            return false;
        }
    }
    
    return true;
}

bool HCS21StateVector::operator!=(const HCS21StateVector& other) const {
    return !(*this == other);
}

// ========== Private Helper Methods ==========

std::vector<Vector3> HCS21StateVector::calculate_compensation_deltas(size_t translation_level, 
                                                                     const Vector3& delta) const {
    std::vector<Vector3> compensation_deltas;
    
    // Number of levels below translation_level that need compensation
    size_t num_lower_levels = NUM_LEVELS - translation_level - 1;
    if (num_lower_levels == 0) {
        return compensation_deltas; // No lower levels to compensate
    }
    
    // Calculate the exact S-change from the translation
    double w_k_squared = scale_weights_[translation_level] * scale_weights_[translation_level];
    Vector3 v_k = level_vectors_[translation_level];
    
    // Exact ΔS from translation: w_k²(‖v_k + δ‖² - ‖v_k‖²)
    Vector3 v_k_new = v_k + delta;
    double delta_S_exact = w_k_squared * (v_k_new.magnitude_squared() - v_k.magnitude_squared());
    
    // Distribute this S-change across lower levels to cancel it out
    // Use simple equal distribution for now (can be improved with Moore-Penrose later)
    
    double target_delta_S_per_level = -delta_S_exact / num_lower_levels;
    
    for (size_t j = translation_level + 1; j < NUM_LEVELS; ++j) {
        double w_j_squared = scale_weights_[j] * scale_weights_[j];
        Vector3 v_j = level_vectors_[j];
        
        // We need: w_j²(‖v_j + Δ_j‖² - ‖v_j‖²) = target_delta_S_per_level
        // Solving: ‖v_j + Δ_j‖² = ‖v_j‖² + target_delta_S_per_level/w_j²
        
        double current_magnitude_squared = v_j.magnitude_squared();
        double required_magnitude_squared = current_magnitude_squared + target_delta_S_per_level / w_j_squared;
        
        if (required_magnitude_squared >= 0.0) {
            double required_magnitude = std::sqrt(required_magnitude_squared);
            double current_magnitude = v_j.magnitude();
            
            if (current_magnitude > 1e-12) {
                // Scale existing vector to required magnitude
                Vector3 scaled_v_j = v_j.normalized() * required_magnitude;
                Vector3 compensation_delta = scaled_v_j - v_j;
                compensation_deltas.push_back(compensation_delta);
            } else {
                // Current vector is near zero, create new vector with required magnitude
                // Use a direction perpendicular to the translation delta
                Vector3 direction = delta.normalized();
                // Create a perpendicular vector
                Vector3 perp_direction;
                if (std::abs(direction.x()) < 0.9) {
                    perp_direction = Vector3(1.0, 0.0, 0.0).cross(direction).normalized();
                } else {
                    perp_direction = Vector3(0.0, 1.0, 0.0).cross(direction).normalized();
                }
                
                Vector3 compensation_delta = perp_direction * required_magnitude;
                compensation_deltas.push_back(compensation_delta);
            }
        } else {
            // Required magnitude squared is negative - this shouldn't happen with equal distribution
            // Use zero compensation as fallback
            compensation_deltas.push_back(Vector3(0.0, 0.0, 0.0));
        }
    }
    
    return compensation_deltas;
}

void HCS21StateVector::apply_rotation_cascade(size_t starting_level, const Matrix4& rotation) {
    // Apply rotation to starting level and all lower levels (j >= starting_level)
    for (size_t j = starting_level; j < NUM_LEVELS; ++j) {
        // Extract 3D vector, apply rotation, store back
        Vector3 rotated_vector = rotation.transform_vector(level_vectors_[j]);
        level_vectors_[j] = rotated_vector;
    }
}

void HCS21StateVector::validate_level_index(size_t level) const {
    if (level >= NUM_LEVELS) {
        throw std::out_of_range("Level index " + std::to_string(level) + 
                               " is out of range [0, " + std::to_string(NUM_LEVELS - 1) + "]");
    }
}

// ========== Free Functions ==========

double calculate_hierarchical_distance(const HCS21StateVector& lhs, 
                                      const HCS21StateVector& rhs) {
    double S_lhs = lhs.calculate_S_invariant();
    double S_rhs = rhs.calculate_S_invariant();
    return std::abs(S_lhs - S_rhs);
}

HCS21StateVector create_galactic_state(const Vector3& galactic_vector) {
    std::array<Vector3, HCS21StateVector::NUM_LEVELS> vectors;
    vectors[0] = galactic_vector;
    for (size_t i = 1; i < HCS21StateVector::NUM_LEVELS; ++i) {
        vectors[i] = Vector3(0.0, 0.0, 0.0);
    }
    
    auto weights = HCS21StateVector::create_default_weights();
    return HCS21StateVector(vectors, weights);
}

HCS21StateVector create_local_state(const Vector3& local_vector) {
    std::array<Vector3, HCS21StateVector::NUM_LEVELS> vectors;
    vectors.fill(Vector3(0.0, 0.0, 0.0));
    vectors[static_cast<size_t>(HCS21StateVector::Level::LOCAL)] = local_vector;
    
    auto weights = HCS21StateVector::create_default_weights();
    return HCS21StateVector(vectors, weights);
}

} // namespace hsml::core