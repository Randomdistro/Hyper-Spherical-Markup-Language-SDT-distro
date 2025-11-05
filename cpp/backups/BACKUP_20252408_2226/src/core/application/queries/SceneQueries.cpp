/** @file SceneQueries.h
 * @brief Query handlers for scene data retrieval
 *
 * Clean Architecture: Application Layer Queries
 * Implements CQRS query pattern for efficient data retrieval.
 */

#pragma once

#include "hsml/domain/interfaces/ISceneRepository.h"
#include "hsml/domain/entities/SphericalScene.h"
#include <memory>
#include <vector>
#include <string>

namespace hsml {
namespace application {

class IQuery {
public:
    virtual ~IQuery() = default;
};

/**
 * @brief Query for getting scene entities
 */
class GetSceneEntitiesQuery : public IQuery {
public:
    explicit GetSceneEntitiesQuery(const std::string& scene_id);

    const std::string& get_scene_id() const { return scene_id_; }

private:
    std::string scene_id_;
};

/**
 * @brief Query for finding entities in radius
 */
class FindEntitiesInRadiusQuery : public IQuery {
public:
    FindEntitiesInRadiusQuery(
        const std::string& scene_id,
        const domain::SphericalCoords& center,
        double radius
    );

    const std::string& get_scene_id() const { return scene_id_; }
    const domain::SphericalCoords& get_center() const { return center_; }
    double get_radius() const { return radius_; }

private:
    std::string scene_id_;
    domain::SphericalCoords center_;
    double radius_;
};

/**
 * @brief Query result for scene entities
 */
struct SceneEntitiesResult {
    bool success{false};
    std::string error_message;
    std::vector<std::shared_ptr<domain::ISphericalEntity>> entities;
};

/**
 * @brief Query result for entities in radius
 */
struct EntitiesInRadiusResult {
    bool success{false};
    std::string error_message;
    std::vector<std::shared_ptr<domain::ISphericalEntity>> entities;
    size_t total_entities_searched{0};
    double search_time_ms{0.0};
};

/**
 * @brief Query handler for scene-related queries
 *
 * Implements CQRS query pattern for efficient data retrieval
 * without side effects.
 */
class SceneQueryHandler {
public:
    explicit SceneQueryHandler(std::shared_ptr<domain::ISceneRepository> repository);

    // Query execution
    SceneEntitiesResult handle(const GetSceneEntitiesQuery& query);
    EntitiesInRadiusResult handle(const FindEntitiesInRadiusQuery& query);

    // Query validation
    bool can_handle(const GetSceneEntitiesQuery& query) const;
    bool can_handle(const FindEntitiesInRadiusQuery& query) const;

private:
    std::shared_ptr<domain::ISceneRepository> repository_;

    // Query optimization methods
    std::vector<std::shared_ptr<domain::ISphericalEntity>> optimize_radius_search(
        const domain::SphericalCoords& center,
        double radius
    );

    // Validation helpers
    bool validate_scene_exists(const std::string& scene_id) const;
    bool validate_search_parameters(const domain::SphericalCoords& center, double radius) const;
};

/**
 * @brief Query result caching for performance
 */
class QueryResultCache {
public:
    void cache_result(const std::string& query_key, const SceneEntitiesResult& result);
    void cache_result(const std::string& query_key, const EntitiesInRadiusResult& result);

    std::optional<SceneEntitiesResult> get_cached_scene_result(const std::string& query_key) const;
    std::optional<EntitiesInRadiusResult> get_cached_radius_result(const std::string& query_key) const;

    void clear_cache();
    void clear_expired_results();

    size_t get_cache_size() const;
    double get_cache_hit_ratio() const;

private:
    struct CachedSceneResult {
        SceneEntitiesResult result;
        std::chrono::steady_clock::time_point timestamp;
    };

    struct CachedRadiusResult {
        EntitiesInRadiusResult result;
        std::chrono::steady_clock::time_point timestamp;
    };

    std::unordered_map<std::string, CachedSceneResult> scene_cache_;
    std::unordered_map<std::string, CachedRadiusResult> radius_cache_;

    mutable size_t cache_hits_{0};
    mutable size_t cache_misses_{0};

    static constexpr std::chrono::seconds CACHE_EXPIRY_TIME{300}; // 5 minutes

    std::string generate_query_key(const GetSceneEntitiesQuery& query) const;
    std::string generate_query_key(const FindEntitiesInRadiusQuery& query) const;
};

} // namespace application
} // namespace hsml
