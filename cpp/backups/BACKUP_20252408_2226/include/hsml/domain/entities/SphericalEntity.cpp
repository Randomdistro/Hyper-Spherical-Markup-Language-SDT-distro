#pragma once

#include "../../../core/math/spherical_coords.h"
#include "../../../core/math/vector3.h"

namespace hsml {
namespace domain {
namespace entities {

class SphericalEntity {
public:
    SphericalEntity() = default;
    SphericalEntity(const core::SphericalCoords& position, double radius);

    const core::SphericalCoords& position() const { return position_; }
    double radius() const { return radius_; }

    void set_position(const core::SphericalCoords& pos) { position_ = pos; }
    void set_radius(double r) { radius_ = r; }

private:
    core::SphericalCoords position_;
    double radius_ = 1.0;
};

} // namespace entities
} // namespace domain
} // namespace hsml
