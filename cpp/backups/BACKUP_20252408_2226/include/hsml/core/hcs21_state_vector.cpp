#pragma once

#include <array>
#include <vector>
#include <memory>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include "vector3.h"
#include "matrix4.h"
#include "compensation_solver.h"

namespace hsml::core {

/**
 * @brief HCS-21 21-Dimensional Hierarchical Coordinate System State Vector
 * 
 * Implements the complete HCS-21 framework with 7 hierarchical levels (21D total).
 * Each level represents a nested physical scale from galactic to sub-quantum.
 * 
 * Core features:
 * - 7 levels × 3 axes = 21 dimensional state space
 * - S-invariant (Total Spation-Displacement Energy) conservation
 * - Legal transforms: rotations with cascade, compensated translations
 * - Multiple projection methods to observable 3D space
 * - Scale weights with strict decreasing property
 * 
 * Levels:
 * 0: Galactic/intergalactic (w₀ = c)
 * 1: Stellar (w₁ = c·k_galaxy)  
 * 2: Planetary/satellite
 * 3: Regional/surface
 * 4: Local/human
 * 5: Atomic/quantum
 * 6: Sub-vortex/sub-quantum
 */
class HCS21StateVector {
public:
    static constexpr size_t NUM_LEVELS = 7;
    static constexpr double SPEED_OF_LIGHT = 299792458.0; // m/s
    static constexpr double DEFAULT_EPS = 1e-12;
    
    // Scale level enumeration for clarity
    enum class Level : size_t {
        GALACTIC = 0,
        STELLAR = 1,
        PLANETARY = 2,
        REGIONAL = 3,
        LOCAL = 4,
        ATOMIC = 5,
        SUB_VORTEX = 6
    };

    /**
     * @brief Default constructor - initializes with zero vectors and default weights
     */
    HCS21StateVector();

    /**
     * @brief Constructor with custom vectors and weights
     * @param level_vectors Array of 7 Vector3 objects representing the 21D state
     * @param scale_weights Array of 7 positive, strictly decreasing weights
     * @param uncertainties Optional array of 7 uncertainty values per level
     */
    HCS21StateVector(const std::array<Vector3, NUM_LEVELS>& level_vectors,
                     const std::array<double, NUM_LEVELS>& scale_weights,
                     const std::array<double, NUM_LEVELS>& uncertainties = {});

    /**
     * @brief Copy constructor
     */
    HCS21StateVector(const HCS21StateVector& other) = default;

    /**
     * @brief Assignment operator
     */
    HCS21StateVector& operator=(const HCS21StateVector& other) = default;

    // ========== Core Properties ==========

    /**
     * @brief Get vector at specific level
     * @param level Level index (0-6)
     * @return Reference to Vector3 at that level
     */
    Vector3& get_level_vector(size_t level);
    const Vector3& get_level_vector(size_t level) const;

    /**
     * @brief Get vector at specific level using enum
     * @param level Level enum value
     * @return Reference to Vector3 at that level
     */
    Vector3& get_level_vector(Level level);
    const Vector3& get_level_vector(Level level) const;

    /**
     * @brief Get weight at specific level
     * @param level Level index (0-6)
     * @return Weight value
     */
    double get_weight(size_t level) const;
    double get_weight(Level level) const;

    /**
     * @brief Get all level vectors
     * @return Const reference to array of vectors
     */
    const std::array<Vector3, NUM_LEVELS>& get_all_vectors() const;

    /**
     * @brief Get all weights
     * @return Const reference to array of weights
     */
    const std::array<double, NUM_LEVELS>& get_all_weights() const;

    // ========== S-Invariant Calculations ==========

    /**
     * @brief Calculate the S-invariant (Total Spation-Displacement Energy)
     * S = Σᵢ wᵢ²‖vᵢ‖²
     * @return The S-invariant value
     */
    double calculate_S_invariant() const;

    /**
     * @brief Calculate displacement potential Φ for a vector
     * v1.0: Φ(v) = ‖v‖²
     * @param vector Input vector
     * @return Displacement potential
     */
    static double displacement_potential(const Vector3& vector);

    /**
     * @brief Validate S-conservation between two states
     * @param other Other HCS21StateVector to compare with
     * @param tolerance Acceptable difference in S-values
     * @return True if S-conservation holds within tolerance
     */
    bool validate_S_conservation(const HCS21StateVector& other, 
                                double tolerance = DEFAULT_EPS) const;

    // ========== Legal Transform Operations ==========

    /**
     * @brief Apply rotation at specified level with cascade to lower levels
     * @param level Level at which to apply rotation (0-6)
     * @param rotation 3x3 rotation matrix (must be in SO(3))
     * @param tolerance S-conservation validation tolerance
     * @throws std::invalid_argument if rotation is not in SO(3)
     * @throws std::runtime_error if S-conservation fails
     */
    void rotate_at_level(size_t level, const Matrix4& rotation, 
                        double tolerance = DEFAULT_EPS);
    void rotate_at_level(Level level, const Matrix4& rotation, 
                        double tolerance = DEFAULT_EPS);

    /**
     * @brief Apply compensated translation at specified level
     * Automatically calculates and applies minimal-norm compensation to lower levels
     * to preserve S-invariant using Moore-Penrose pseudoinverse method
     * @param level Level at which to apply translation (0-6)
     * @param delta Translation vector
     * @param tolerance S-conservation validation tolerance
     * @throws std::runtime_error if compensation fails or S-conservation violated
     */
    void translate_at_level(size_t level, const Vector3& delta, 
                           double tolerance = DEFAULT_EPS);
    void translate_at_level(Level level, const Vector3& delta, 
                           double tolerance = DEFAULT_EPS);

    // ========== Projection Methods to 3D ==========

    /**
     * @brief Weighted linear projection to 3D observable space
     * Π³ᴰʷᵉⁱᵍʰᵗᵉᵈ(L) = Σᵢ wᵢ vᵢ
     * @return 3D projection vector
     */
    Vector3 project_weighted() const;

    /**
     * @brief Dominant scale projection (hard cutoff)
     * Returns vector at finest scale where wᵢ‖vᵢ‖ > ε_obs
     * @param eps_obs Observer precision threshold
     * @return 3D projection vector from dominant scale
     */
    Vector3 project_dominant(double eps_obs) const;

    /**
     * @brief Smooth dominant projection with tanh transitions
     * Uses tanh function for continuous transitions between scales
     * @param eps_obs Observer precision threshold
     * @param width Transition width parameter (default 0.1)
     * @return 3D smooth projection vector
     */
    Vector3 project_dominant_smooth(double eps_obs, double width = 0.1) const;

    // ========== Uncertainty Propagation ==========

    /**
     * @brief Set uncertainty values for all levels
     * @param uncertainties Array of 7 uncertainty values
     */
    void set_uncertainties(const std::array<double, NUM_LEVELS>& uncertainties);

    /**
     * @brief Get uncertainty at specific level
     * @param level Level index (0-6)
     * @return Uncertainty value
     */
    double get_uncertainty(size_t level) const;

    /**
     * @brief Propagate uncertainties to calculate S-invariant uncertainty
     * σ_S = 2√(Σᵢ wᵢ⁴‖vᵢ‖²σᵢ²)
     * @return Propagated uncertainty in S-value
     */
    double propagate_uncertainty() const;

    // ========== Validation & Utilities ==========

    /**
     * @brief Validate that weights are strictly decreasing and positive
     * @param weights Array of weights to validate
     * @return True if weights are valid
     */
    static bool validate_weights(const std::array<double, NUM_LEVELS>& weights);

    /**
     * @brief Validate that a matrix is in SO(3) (special orthogonal group)
     * @param matrix Matrix to validate
     * @param tolerance Numerical tolerance for validation
     * @return True if matrix is a valid rotation
     */
    static bool validate_rotation_matrix(const Matrix4& matrix, 
                                       double tolerance = DEFAULT_EPS);

    /**
     * @brief Create default scale weights with realistic physical values
     * @param k_galaxy Galaxy-to-stellar scale ratio
     * @param k_stellar Stellar-to-planetary scale ratio  
     * @param k_planetary Planetary-to-regional scale ratio
     * @param k_regional Regional-to-local scale ratio
     * @param k_local Local-to-atomic scale ratio
     * @param k_atomic Atomic-to-subvortex scale ratio
     * @return Array of calculated scale weights
     */
    static std::array<double, NUM_LEVELS> create_default_weights(
        double k_galaxy = 1e-3,    // c → ~300 km/s (1e-3)
        double k_stellar = 0.1,    // Stellar scale ratio (0.1)
        double k_planetary = 0.01, // Planetary scale ratio (0.01)
        double k_regional = 0.1,   // Regional scale ratio (0.1)
        double k_local = 0.1,      // Local scale ratio (0.1)
        double k_atomic = 0.1      // Atomic scale ratio (0.1)
    );

    // ========== Operators ==========

    /**
     * @brief Addition operator for HCS21StateVector
     */
    HCS21StateVector operator+(const HCS21StateVector& other) const;

    /**
     * @brief Subtraction operator for HCS21StateVector
     */
    HCS21StateVector operator-(const HCS21StateVector& other) const;

    /**
     * @brief Scalar multiplication operator
     */
    HCS21StateVector operator*(double scalar) const;

    /**
     * @brief Equality operator
     */
    bool operator==(const HCS21StateVector& other) const;

    /**
     * @brief Inequality operator
     */
    bool operator!=(const HCS21StateVector& other) const;

private:
    // Core 21D state: 7 levels × 3 axes
    std::array<Vector3, NUM_LEVELS> level_vectors_;
    
    // Scale weights (strictly decreasing, positive)
    std::array<double, NUM_LEVELS> scale_weights_;
    
    // Optional uncertainty values per level
    std::array<double, NUM_LEVELS> uncertainties_;
    bool has_uncertainties_;

    // ========== Private Helper Methods ==========

    /**
     * @brief Calculate compensation deltas for translation using Moore-Penrose pseudoinverse
     * @param translation_level Level where translation is applied
     * @param delta Translation vector
     * @return Array of compensation vectors for levels > translation_level
     */
    std::vector<Vector3> calculate_compensation_deltas(size_t translation_level, 
                                                      const Vector3& delta) const;

    /**
     * @brief Apply rotation cascade to lower levels
     * @param starting_level Level where rotation starts
     * @param rotation Rotation matrix to apply
     */
    void apply_rotation_cascade(size_t starting_level, const Matrix4& rotation);

    /**
     * @brief Validate level index
     * @param level Level index to validate
     * @throws std::out_of_range if level is invalid
     */
    void validate_level_index(size_t level) const;
};

// ========== Free Functions ==========

/**
 * @brief Calculate distance between two HCS21StateVector objects
 * d(L,M)² = |S(L) - S(M)|
 * @param lhs First state vector
 * @param rhs Second state vector
 * @return Distance value
 */
double calculate_hierarchical_distance(const HCS21StateVector& lhs, 
                                     const HCS21StateVector& rhs);

/**
 * @brief Create HCS21StateVector with galactic-scale initialization
 * @param galactic_vector Initial galactic-level vector
 * @return Initialized HCS21StateVector
 */
HCS21StateVector create_galactic_state(const Vector3& galactic_vector);

/**
 * @brief Create HCS21StateVector with local-scale initialization
 * @param local_vector Initial local-level vector
 * @return Initialized HCS21StateVector
 */
HCS21StateVector create_local_state(const Vector3& local_vector);

} // namespace hsml::core