/** @file LinearMovementStrategy.h
 * @brief Linear movement strategy implementation
 *
 * Clean Architecture: Domain Layer Entity
 * Implements IMovementStrategy for constant velocity movement.
 */

#pragma once

#include "hsml/domain/interfaces/IMovementStrategy.h"
#include <memory>
#include <string>
#include <vector>

namespace hsml {
namespace domain {

/**
 * @brief Linear movement strategy for constant velocity motion
 *
 * Implements straight-line movement with constant velocity.
 * Velocity is defined in spherical coordinate deltas per second.
 */
class LinearMovementStrategy : public IMovementStrategy {
public:
    /**
     * @brief Construct linear movement strategy
     *
     * @param velocity_r Velocity in radius direction (units per second)
     * @param velocity_theta Velocity in theta direction (radians per second)
     * @param velocity_phi Velocity in phi direction (radians per second)
     */
    LinearMovementStrategy(
        double velocity_r = 0.0,
        double velocity_theta = 0.0,
        double velocity_phi = 0.0
    );

    // IMovementStrategy implementation
    std::string get_name() const override { return "linear"; }
    std::string get_description() const override;

    SphericalCoords calculate_next_position(
        const ISphericalEntity& entity,
        double delta_time
    ) const override;

    bool is_valid_position(const SphericalCoords& position) const override;
    double get_max_speed() const override;
    double get_acceleration() const override { return 0.0; }  // No acceleration

    void set_parameter(const std::string& name, double value) override;
    double get_parameter(const std::string& name) const override;
    std::vector<std::string> get_parameter_names() const override;

    std::unique_ptr<IMovementStrategy> clone() const override;

    // Linear movement specific methods
    void set_velocity(double velocity_r, double velocity_theta, double velocity_phi);
    void set_velocity_r(double velocity_r);
    void set_velocity_theta(double velocity_theta);
    void set_velocity_phi(double velocity_phi);

    double get_velocity_r() const { return velocity_r_; }
    double get_velocity_theta() const { return velocity_theta_; }
    double get_velocity_phi() const { return velocity_phi_; }

private:
    double velocity_r_;      // Units per second
    double velocity_theta_;  // Radians per second
    double velocity_phi_;    // Radians per second

    // Parameter names for interface compliance
    static const std::vector<std::string> parameter_names_;

    // Helper methods
    double calculate_speed() const;
    SphericalCoords apply_velocity(
        const SphericalCoords& position,
        double delta_time
    ) const;
};

} // namespace domain
} // namespace hsml
