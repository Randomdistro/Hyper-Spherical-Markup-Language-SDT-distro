/** @file ISphericalEntity.h
 * @brief Domain interface for spherical entities
 *
 * Clean Architecture: Domain Layer Interface
 * Defines the contract for spherical entities in the business domain.
 */

#pragma once

#include <memory>
#include <string>

namespace hsml {
namespace domain {

// Forward declarations
class SphericalCoords;
class SphericalShape;
class IMovementStrategy;

/**
 * @brief Interface for spherical entities in the domain layer
 *
 * This interface defines the contract for all spherical entities
 * in the HSML system, providing a clean abstraction from infrastructure
 * concerns like rendering and persistence.
 */
class ISphericalEntity {
public:
    virtual ~ISphericalEntity() = default;

    // Entity identification
    virtual std::string get_id() const = 0;
    virtual std::string get_name() const = 0;
    virtual void set_name(const std::string& name) = 0;

    // Spatial properties (pure spherical)
    virtual const SphericalCoords& position() const = 0;
    virtual const SphericalShape& shape() const = 0;

    // Movement operations
    virtual void move_to(const SphericalCoords& target) = 0;
    virtual void move_by(double radius_delta, double theta_delta, double phi_delta) = 0;
    virtual bool collides_with(const ISphericalEntity& other) const = 0;

    // Movement strategy
    virtual void set_movement_strategy(std::unique_ptr<IMovementStrategy> strategy) = 0;
    virtual const IMovementStrategy* get_movement_strategy() const = 0;

    // State management
    virtual bool is_active() const = 0;
    virtual void set_active(bool active) = 0;

    // Domain-specific queries
    virtual double get_mass() const = 0;      // For physics calculations
    virtual double get_charge() const = 0;    // For electromagnetic interactions
    virtual double get_energy() const = 0;    // For quantum field calculations

    // Clone for prototyping
    virtual std::unique_ptr<ISphericalEntity> clone() const = 0;
};

} // namespace domain
} // namespace hsml
