/** @file IMovementStrategy.h
 * @brief Movement strategy interface for spherical entities
 *
 * Clean Architecture: Domain Layer Interface
 * Defines movement strategies for entity motion in spherical space.
 */

#pragma once

#include <memory>
#include <string>

namespace hsml {
namespace domain {

class SphericalCoords;
class ISphericalEntity;

/**
 * @brief Interface for movement strategies
 *
 * Defines different movement patterns and behaviors for spherical entities.
 * This follows the Strategy pattern to allow runtime switching of movement
 * algorithms without changing the entity implementation.
 */
class IMovementStrategy {
public:
    virtual ~IMovementStrategy() = default;

    // Strategy identification
    virtual std::string get_name() const = 0;
    virtual std::string get_description() const = 0;

    // Movement calculation
    virtual SphericalCoords calculate_next_position(
        const ISphericalEntity& entity,
        double delta_time
    ) const = 0;

    // Movement constraints and validation
    virtual bool is_valid_position(const SphericalCoords& position) const = 0;
    virtual double get_max_speed() const = 0;
    virtual double get_acceleration() const = 0;

    // Strategy parameters
    virtual void set_parameter(const std::string& name, double value) = 0;
    virtual double get_parameter(const std::string& name) const = 0;
    virtual std::vector<std::string> get_parameter_names() const = 0;

    // Clone for entity assignment
    virtual std::unique_ptr<IMovementStrategy> clone() const = 0;
};

/**
 * @brief Movement strategy types
 */
enum class MovementStrategyType {
    STATIC,           // No movement
    LINEAR,           // Constant velocity
    ORBITAL,          // Circular/elliptical orbits
    GRAVITATIONAL,    // Newtonian gravity
    SPRING_DAMPED,    // Spring-mass system
    FLUID_DYNAMIC,    // Fluid flow following
    RANDOM_WALK,      // Brownian motion
    PATH_FOLLOWING,   // Follow predefined path
    AI_CONTROLLED     // AI-driven movement
};

/**
 * @brief Movement constraints for safety and physics
 */
struct MovementConstraints {
    double max_radius{1e6};           // Maximum allowed radius
    double min_radius{1e-6};          // Minimum allowed radius
    double max_speed{1e8};            // Speed limit (c)
    double max_acceleration{1e12};    // Acceleration limit
    bool enforce_spherical_bounds{true}; // Keep within valid spherical coordinates
    bool prevent_collisions{true};    // Basic collision avoidance
};

} // namespace domain
} // namespace hsml
