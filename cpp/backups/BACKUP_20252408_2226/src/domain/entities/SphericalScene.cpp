/** @file SphericalScene.h
 * @brief Domain entity representing a spherical scene
 *
 * Clean Architecture: Domain Layer Entity
 * Represents a collection of spherical entities in a scene.
 */

#pragma once

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>

namespace hsml {
namespace domain {

class ISphericalEntity;
class SphericalCoords;

/**
 * @brief Domain entity representing a spherical scene
 *
 * A SphericalScene contains a collection of spherical entities
 * and provides operations for managing and querying the scene.
 * This is a domain entity following Clean Architecture principles.
 */
class SphericalScene {
public:
    using EntityId = std::string;
    using EntityPtr = std::shared_ptr<ISphericalEntity>;
    using EntityMap = std::unordered_map<EntityId, EntityPtr>;

    /**
     * @brief Construct an empty spherical scene
     */
    SphericalScene() = default;

    /**
     * @brief Construct a scene with initial entities
     */
    explicit SphericalScene(const std::vector<EntityPtr>& entities);

    // Scene identification
    const std::string& get_id() const { return id_; }
    const std::string& get_name() const { return name_; }
    void set_name(const std::string& name) { name_ = name; }

    // Entity management
    bool add_entity(EntityPtr entity);
    bool remove_entity(const EntityId& entity_id);
    EntityPtr get_entity(const EntityId& entity_id) const;
    const EntityMap& get_all_entities() const { return entities_; }

    // Scene queries
    std::vector<EntityPtr> find_entities_in_radius(
        const SphericalCoords& center,
        double radius
    ) const;

    std::vector<EntityPtr> find_entities_by_type(const std::string& type) const;

    std::vector<EntityPtr> find_colliding_entities(
        const ISphericalEntity& entity
    ) const;

    // Scene properties
    size_t get_entity_count() const { return entities_.size(); }
    bool is_empty() const { return entities_.empty(); }

    // Scene bounds (spherical)
    double get_max_radius() const;
    SphericalCoords get_centroid() const;

    // Scene operations
    void clear();
    void merge(const SphericalScene& other);

    // Scene metadata
    void set_attribute(const std::string& key, const std::string& value);
    std::string get_attribute(const std::string& key) const;
    bool has_attribute(const std::string& key) const;

private:
    std::string id_;
    std::string name_;
    EntityMap entities_;
    std::unordered_map<std::string, std::string> attributes_;

    static std::string generate_id();
};

} // namespace domain
} // namespace hsml
