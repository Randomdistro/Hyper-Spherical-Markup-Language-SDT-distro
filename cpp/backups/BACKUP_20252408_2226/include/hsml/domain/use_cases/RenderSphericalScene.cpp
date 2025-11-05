#pragma once

#include "../entities/SphericalScene.h"
#include "../../../core/math/spherical_coords.h"

namespace hsml {
namespace domain {
namespace use_cases {

class RenderSphericalScene {
public:
    RenderSphericalScene() = default;

    void execute(const SphericalScene& scene);

private:
    // Implementation details for rendering use case
};

} // namespace use_cases
} // namespace domain
} // namespace hsml
