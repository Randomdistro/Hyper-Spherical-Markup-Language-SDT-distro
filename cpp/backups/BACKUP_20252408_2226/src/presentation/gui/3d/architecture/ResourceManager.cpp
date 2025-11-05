/** @file ResourceManager.h
 * @brief Resource Manager Implementation - Asset Management System
 */

#pragma once

#include "src/presentation/gui/3d/architecture/ArchitectureCore.h"
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <mutex>
#include <shared_mutex>
#include <atomic>
#include <string>
#include <vector>
#include <chrono>
#include <future>

namespace hsml {
namespace gui3d {
namespace architecture {

/**
 * @brief Resource Base Class
 */
class Resource {
public:
    Resource(const std::string& path, const std::string& type);
    virtual ~Resource() = default;

    const std::string& get_path() const { return path_; }
    const std::string& get_type() const { return type_; }
    size_t get_size() const { return size_; }
    bool is_loaded() const { return loaded_; }
    std::chrono::steady_clock::time_point get_load_time() const { return load_time_; }
    std::chrono::steady_clock::time_point get_last_access_time() const { return last_access_time_; }

    void set_loaded(bool loaded) { loaded_ = loaded; }
    void set_size(size_t size) { size_ = size; }
    void update_access_time();

    virtual bool load() = 0;
    virtual void unload() = 0;
    virtual size_t get_memory_usage() const = 0;

protected:
    std::string path_;
    std::string type_;
    size_t size_;
    bool loaded_;
    std::chrono::steady_clock::time_point load_time_;
    std::chrono::steady_clock::time_point last_access_time_;
};

/**
 * @brief Texture Resource
 */
class TextureResource : public Resource {
public:
    TextureResource(const std::string& path);
    ~TextureResource() override;

    bool load() override;
    void unload() override;
    size_t get_memory_usage() const override;

    void* get_texture_data() const { return texture_data_; }
    int get_width() const { return width_; }
    int get_height() const { return height_; }
    int get_channels() const { return channels_; }

private:
    void* texture_data_;
    int width_;
    int height_;
    int channels_;
};

/**
 * @brief Mesh Resource
 */
class MeshResource : public Resource {
public:
    MeshResource(const std::string& path);
    ~MeshResource() override;

    bool load() override;
    void unload() override;
    size_t get_memory_usage() const override;

    const std::vector<float>& get_vertices() const { return vertices_; }
    const std::vector<float>& get_normals() const { return normals_; }
    const std::vector<float>& get_texcoords() const { return texcoords_; }
    const std::vector<unsigned int>& get_indices() const { return indices_; }

private:
    std::vector<float> vertices_;
    std::vector<float> normals_;
    std::vector<float> texcoords_;
    std::vector<unsigned int> indices_;
};

/**
 * @brief Shader Resource
 */
class ShaderResource : public Resource {
public:
    ShaderResource(const std::string& path);
    ~ShaderResource() override;

    bool load() override;
    void unload() override;
    size_t get_memory_usage() const override;

    void* get_program() const { return program_; }
    bool is_compiled() const { return compiled_; }

private:
    void* program_;
    bool compiled_;
};

/**
 * @brief Resource Manager Implementation
 */
class ResourceManager : public IResourceManager {
public:
    ResourceManager();
    ~ResourceManager() override;

    bool initialize() override;
    void shutdown() override;

    // Resource loading
    bool load_resource(const std::string& path, const std::string& type) override;
    bool unload_resource(const std::string& path) override;
    bool is_resource_loaded(const std::string& path) const override;

    // Resource access
    void* get_resource(const std::string& path) const override;
    size_t get_resource_size(const std::string& path) const override;
    std::string get_resource_type(const std::string& path) const override;

    // Resource streaming
    bool start_streaming(const std::string& path) override;
    bool stop_streaming(const std::string& path) override;
    float get_streaming_progress(const std::string& path) const override;

    // Memory management
    size_t get_total_memory_usage() const override;
    size_t get_cache_size() const override;
    void set_memory_limit(size_t bytes) override;
    void optimize_memory() override;

private:
    struct StreamingResource {
        std::shared_ptr<Resource> resource;
        std::future<bool> loading_future;
        float progress;
        bool active;
    };

    // Core data structures
    std::unordered_map<std::string, std::shared_ptr<Resource>> loaded_resources_;
    std::unordered_map<std::string, StreamingResource> streaming_resources_;
    std::vector<std::string> resource_load_order_;

    // Thread safety
    mutable std::shared_mutex resources_mutex_;
    mutable std::shared_mutex streaming_mutex_;

    // Configuration
    size_t memory_limit_;
    size_t current_memory_usage_;
    std::atomic<bool> initialized_;

    // Cache management
    std::unordered_set<std::string> pinned_resources_;
    std::vector<std::string> access_history_;

    // Helper methods
    std::shared_ptr<Resource> create_resource(const std::string& path, const std::string& type);
    bool load_resource_internal(const std::string& path, const std::string& type);
    void unload_resource_internal(const std::string& path);
    void update_memory_usage();
    void enforce_memory_limit();
    void optimize_cache();

    // Streaming helpers
    bool start_streaming_internal(const std::string& path);
    void stop_streaming_internal(const std::string& path);
    void update_streaming_progress();

    // Constants
    static const size_t DEFAULT_MEMORY_LIMIT;
    static const std::chrono::milliseconds CACHE_CLEANUP_INTERVAL;
};

/**
 * @brief Resource Cache - LRU Cache Implementation
 */
class ResourceCache {
public:
    ResourceCache(size_t max_size);
    ~ResourceCache() = default;

    bool put(const std::string& key, std::shared_ptr<Resource> resource);
    std::shared_ptr<Resource> get(const std::string& key);
    bool remove(const std::string& key);
    void clear();

    size_t size() const { return cache_.size(); }
    size_t max_size() const { return max_size_; }
    size_t memory_usage() const { return memory_usage_; }

    void set_max_size(size_t max_size) { max_size_ = max_size; }
    void optimize();

private:
    struct CacheEntry {
        std::shared_ptr<Resource> resource;
        std::chrono::steady_clock::time_point last_access;
        size_t access_count;
    };

    std::unordered_map<std::string, CacheEntry> cache_;
    std::list<std::string> access_order_;
    size_t max_size_;
    size_t memory_usage_;

    mutable std::shared_mutex cache_mutex_;

    void evict_lru();
    void update_access_order(const std::string& key);
};

} // namespace architecture
} // namespace gui3d
} // namespace hsml
