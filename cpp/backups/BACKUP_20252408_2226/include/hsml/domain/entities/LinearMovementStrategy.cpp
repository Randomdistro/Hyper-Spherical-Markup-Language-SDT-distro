#pragma once

#include "../../../core/math/spherical_coords.h"
#include "../../../core/math/vector3.h"

namespace hsml {
namespace domain {
namespace entities {

class LinearMovementStrategy {
public:
    LinearMovementStrategy() = default;
    LinearMovementStrategy(const core::Vector3& velocity, double speed = 1.0);

    core::Vector3 velocity() const { return velocity_; }
    double speed() const { return speed_; }

    void set_velocity(const core::Vector3& vel) { velocity_ = vel; }
    void set_speed(double s) { speed_ = s; }

private:
    core::Vector3 velocity_;
    double speed_ = 1.0;
};

} // namespace entities
} // namespace domain
} // namespace hsml
