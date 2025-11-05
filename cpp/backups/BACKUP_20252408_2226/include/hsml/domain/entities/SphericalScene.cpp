#pragma once

#include "SphericalEntity.h"
#include <vector>
#include <memory>

namespace hsml {
namespace domain {
namespace entities {

class SphericalScene {
public:
    SphericalScene() = default;

    void add_entity(std::shared_ptr<SphericalEntity> entity);
    void remove_entity(const std::shared_ptr<SphericalEntity>& entity);
    const std::vector<std::shared_ptr<SphericalEntity>>& entities() const;

private:
    std::vector<std::shared_ptr<SphericalEntity>> entities_;
};

} // namespace entities
} // namespace domain
} // namespace hsml
