/** @file ISceneRepository.h
 * @brief Repository interface for scene data access
 *
 * Clean Architecture: Domain Layer Interface
 * Defines the contract for accessing scene data from persistence.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

namespace hsml {
namespace domain {

class SphericalScene;
class ISphericalEntity;

/**
 * @brief Interface for scene data repository
 *
 * Defines the contract for accessing and storing spherical scene data.
 * This interface follows the Repository pattern to abstract data access
 * from the domain layer.
 */
class ISceneRepository {
public:
    virtual ~ISceneRepository() = default;

    // Scene operations
    virtual SphericalScene get_by_id(const std::string& scene_id) = 0;
    virtual std::vector<SphericalScene> get_all_scenes() = 0;
    virtual void save(const SphericalScene& scene) = 0;
    virtual void delete_scene(const std::string& scene_id) = 0;
    virtual bool scene_exists(const std::string& scene_id) const = 0;

    // Entity operations within scenes
    virtual std::shared_ptr<ISphericalEntity> get_entity(
        const std::string& scene_id,
        const std::string& entity_id
    ) = 0;

    virtual std::vector<std::shared_ptr<ISphericalEntity>> get_entities_in_scene(
        const std::string& scene_id
    ) = 0;

    virtual void add_entity_to_scene(
        const std::string& scene_id,
        std::shared_ptr<ISphericalEntity> entity
    ) = 0;

    virtual void remove_entity_from_scene(
        const std::string& scene_id,
        const std::string& entity_id
    ) = 0;

    // Query operations
    virtual std::vector<std::shared_ptr<ISphericalEntity>> find_entities_in_radius(
        const std::string& scene_id,
        const SphericalCoords& center,
        double radius
    ) = 0;

    virtual std::vector<std::shared_ptr<ISphericalEntity>> find_entities_by_type(
        const std::string& scene_id,
        const std::string& type
    ) = 0;

    // Transaction support
    virtual void begin_transaction() = 0;
    virtual void commit_transaction() = 0;
    virtual void rollback_transaction() = 0;
};

} // namespace domain
} // namespace hsml
