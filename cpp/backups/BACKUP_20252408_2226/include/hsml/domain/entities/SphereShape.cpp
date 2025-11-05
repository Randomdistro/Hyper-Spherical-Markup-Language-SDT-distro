#pragma once

#include "SphericalEntity.h"

namespace hsml {
namespace domain {
namespace entities {

class SphereShape : public SphericalEntity {
public:
    SphereShape() = default;
    SphereShape(const core::SphericalCoords& position, double radius);

    // Additional sphere-specific properties can be added here
};

} // namespace entities
} // namespace domain
} // namespace hsml
