/** @file SphericalSceneRepository.h
 * @brief Infrastructure layer repository for spherical scenes
 *
 * Clean Architecture: Infrastructure Layer
 * Implements ISceneRepository using existing Bubble infrastructure.
 */

#pragma once

#include "hsml/domain/interfaces/ISceneRepository.h"
#include "hsml/domain/entities/SphericalScene.h"
#include <memory>
#include <unordered_map>
#include <mutex>

namespace hsml {
namespace infrastructure {

/**
 * @brief Repository implementation for spherical scenes
 *
 * Implements ISceneRepository using the existing Bubble infrastructure
 * while maintaining clean architecture boundaries.
 */
class SphericalSceneRepository : public domain::ISceneRepository {
public:
    SphericalSceneRepository();
    ~SphericalSceneRepository() override = default;

    // ISceneRepository implementation
    domain::SphericalScene get_by_id(const std::string& scene_id) override;
    std::vector<domain::SphericalScene> get_all_scenes() override;
    void save(const domain::SphericalScene& scene) override;
    void delete_scene(const std::string& scene_id) override;
    bool scene_exists(const std::string& scene_id) const override;

    // Entity operations
    std::shared_ptr<domain::ISphericalEntity> get_entity(
        const std::string& scene_id,
        const std::string& entity_id
    ) override;

    std::vector<std::shared_ptr<domain::ISphericalEntity>> get_entities_in_scene(
        const std::string& scene_id
    ) override;

    void add_entity_to_scene(
        const std::string& scene_id,
        std::shared_ptr<domain::ISphericalEntity> entity
    ) override;

    void remove_entity_from_scene(
        const std::string& scene_id,
        const std::string& entity_id
    ) override;

    // Query operations
    std::vector<std::shared_ptr<domain::ISphericalEntity>> find_entities_in_radius(
        const std::string& scene_id,
        const domain::SphericalCoords& center,
        double radius
    ) override;

    std::vector<std::shared_ptr<domain::ISphericalEntity>> find_entities_by_type(
        const std::string& scene_id,
        const std::string& type
    ) override;

    // Transaction support
    void begin_transaction() override;
    void commit_transaction() override;
    void rollback_transaction() override;

private:
    // Storage using existing infrastructure
    struct SceneData {
        domain::SphericalScene scene;
        std::unordered_map<std::string, std::shared_ptr<domain::ISphericalEntity>> entities;
        bool is_dirty{false};
    };

    std::unordered_map<std::string, SceneData> scenes_;
    mutable std::mutex mutex_;

    // Transaction support
    std::unordered_map<std::string, SceneData> transaction_backup_;
    bool in_transaction_{false};

    // Helper methods
    SceneData& get_or_create_scene_data(const std::string& scene_id);
    const SceneData& get_scene_data(const std::string& scene_id) const;

    std::shared_ptr<domain::ISphericalEntity> convert_bubble_to_entity(
        std::shared_ptr<core::Bubble> bubble
    ) const;

    std::shared_ptr<core::Bubble> convert_entity_to_bubble(
        std::shared_ptr<domain::ISphericalEntity> entity
    ) const;

    // Persistence methods (simplified in-memory for now)
    void load_scene(const std::string& scene_id);
    void persist_scene(const std::string& scene_id);
    void ensure_default_scene_exists();
};

/**
 * @brief Repository factory for dependency injection
 */
class RepositoryFactory {
public:
    static std::shared_ptr<domain::ISceneRepository> create_scene_repository();
};

} // namespace infrastructure
} // namespace hsml
