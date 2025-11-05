/** @file SphericalEntity.h
 * @brief Core domain entity for spherical objects
 *
 * Clean Architecture: Domain Layer Entity
 * Implements ISphericalEntity using the Bubble infrastructure
 * while maintaining clean domain boundaries.
 */

#pragma once

#include "hsml/core/bubble.h"
#include "hsml/domain/interfaces/ISphericalEntity.h"
#include "hsml/domain/interfaces/ISphericalShape.h"
#include "hsml/domain/interfaces/IMovementStrategy.h"
#include <memory>
#include <string>
#include <chrono>

namespace hsml {
namespace domain {

/**
 * @brief Core domain entity for spherical objects
 *
 * This class implements ISphericalEntity and wraps the infrastructure
 * Bubble class while maintaining clean domain boundaries. It provides
 * domain-specific operations and enforces business rules.
 */
class SphericalEntity : public ISphericalEntity {
public:
    /**
     * @brief Construct a spherical entity
     *
     * @param position Initial spherical position
     * @param shape The entity's shape
     * @param name Optional entity name
     */
    explicit SphericalEntity(
        const SphericalCoords& position,
        std::unique_ptr<ISphericalShape> shape,
        const std::string& name = ""
    );

    /**
     * @brief Create entity from existing Bubble (adapter pattern)
     */
    explicit SphericalEntity(std::shared_ptr<core::Bubble> bubble);

    // ISphericalEntity implementation
    std::string get_id() const override;
    std::string get_name() const override;
    void set_name(const std::string& name) override;

    const SphericalCoords& position() const override;
    const ISphericalShape& shape() const override;

    void move_to(const SphericalCoords& target) override;
    void move_by(double radius_delta, double theta_delta, double phi_delta) override;
    bool collides_with(const ISphericalEntity& other) const override;

    void set_movement_strategy(std::unique_ptr<IMovementStrategy> strategy) override;
    const IMovementStrategy* get_movement_strategy() const override;

    bool is_active() const override;
    void set_active(bool active) override;

    double get_mass() const override;
    double get_charge() const override;
    double get_energy() const override;

    std::unique_ptr<ISphericalEntity> clone() const override;

    // Additional domain operations
    void update(double delta_time);
    bool is_visible_from(const SphericalCoords& observer) const;

    // Physical properties
    void set_mass(double mass);
    void set_charge(double charge);
    void set_energy(double energy);

    // Access to underlying infrastructure (for infrastructure layer only)
    std::shared_ptr<core::Bubble> get_infrastructure_bubble() const {
        return bubble_;
    }

private:
    std::shared_ptr<core::Bubble> bubble_;
    std::unique_ptr<ISphericalShape> shape_;
    std::unique_ptr<IMovementStrategy> movement_strategy_;

    std::string name_;
    bool is_active_;

    // Physical properties
    double mass_{1.0};
    double charge_{0.0};
    double energy_{0.0};

    // Entity state
    std::chrono::steady_clock::time_point last_update_;
    SphericalCoords last_position_;

    // Business rule validation
    bool validate_position(const SphericalCoords& position) const;
    bool validate_movement(const SphericalCoords& from, const SphericalCoords& to) const;

    // Internal operations
    void update_position(const SphericalCoords& new_position);
    void notify_position_change();
};

} // namespace domain
} // namespace hsml
