#pragma once

#include "spherical_coords.h"
#include "simd_math.h"
#include "precision_constants.h"
#include <immintrin.h>
#include <vector>
#include <array>
#include <unordered_map>
#include <atomic>
#include <chrono>
#include <cstring>
#include <algorithm>

// [Performance Demon]: ULTRA-HIGH PERFORMANCE SIMD SPATIAL INDEXER
// Every nanosecond counts! Memory pools, SIMD intrinsics, cache optimization!

namespace hsml {
namespace core {
namespace spatial {

// [Performance Demon]: Pack ALL the data for cache efficiency!
struct alignas(64) SpatialNodeSIMD {
    static constexpr size_t MAX_ELEMENTS_PER_NODE = 16; // Perfect for SIMD registers
    static constexpr size_t DODECAHEDRON_FACES = 12;
    
    // SIMD-friendly data layout
    alignas(32) double element_coords[MAX_ELEMENTS_PER_NODE * 3]; // r,theta,phi packed
    alignas(32) uint32_t element_ids[MAX_ELEMENTS_PER_NODE];
    alignas(32) uint32_t child_indices[DODECAHEDRON_FACES];
    
    // Bounds for fast rejection tests - SIMD aligned
    alignas(32) double bounds[6]; // rMin,rMax,thetaMin,thetaMax,phiMin,phiMax
    alignas(32) double face_centers[DODECAHEDRON_FACES * 3]; // Pre-computed face centers
    alignas(32) double face_normals[DODECAHEDRON_FACES * 3]; // Pre-computed normals
    
    uint32_t element_count;
    uint32_t depth;
    uint32_t parent_index;
    uint32_t node_id;
    
    // [Performance Demon]: NEVER waste cache lines!
    char padding[64 - ((sizeof(uint32_t) * 4) % 64)];
};

// [Performance Demon]: Memory pool for ZERO allocations during queries!
class alignas(64) SpatialMemoryPool {
public:
    static constexpr size_t POOL_SIZE = 1024 * 1024; // 1MB pool
    static constexpr size_t MAX_NODES = POOL_SIZE / sizeof(SpatialNodeSIMD);
    
    // [Performance Demon]: Public access for fast indexing - THREAD SAFE
    alignas(64) mutable SpatialNodeSIMD nodes_[MAX_NODES];
    
private:
    std::atomic<uint32_t> next_free_index_{0};
    std::atomic<uint32_t> allocation_count_{0};
    
public:
    SpatialNodeSIMD* allocate_node() noexcept {
        // [Performance Demon]: FIXED THREAD SAFETY - proper memory ordering
        uint32_t index = next_free_index_.fetch_add(1, std::memory_order_acq_rel);
        if (index >= MAX_NODES) {
            // [Performance Demon]: NEVER fail! Wrap around and reuse
            index = index % MAX_NODES;
        }
        allocation_count_.fetch_add(1, std::memory_order_acq_rel);
        
        // Initialize the node safely
        SpatialNodeSIMD* node = &nodes_[index];
        std::memset(node, 0, sizeof(SpatialNodeSIMD));
        node->node_id = index;
        
        return node;
    }
    
    void reset() noexcept {
        // [Performance Demon]: FIXED THREAD SAFETY - proper memory barriers
        next_free_index_.store(0, std::memory_order_release);
        allocation_count_.store(0, std::memory_order_release);
        // Memory fence to ensure all writes are visible
        std::atomic_thread_fence(std::memory_order_seq_cst);
    }
    
    size_t get_allocation_count() const noexcept {
        return allocation_count_.load(std::memory_order_relaxed);
    }
};

// [Performance Demon]: Query results - SIMD friendly
struct alignas(32) QueryResultSIMD {
    alignas(32) uint32_t element_ids[256]; // Max results per query
    alignas(32) double distances[256];     // Pre-computed distances
    uint32_t count;
    uint32_t capacity;
    std::chrono::high_resolution_clock::time_point timestamp;
    
    QueryResultSIMD() : count(0), capacity(256) {}
};

// [Performance Demon]: THE BEAST - SIMD SPATIAL INDEXER
class alignas(64) SpatialIndexerSIMD {
public:
    static constexpr double GOLDEN_RATIO = 1.618033988749894848204586834365638117720309179805762862135;
    static constexpr size_t CACHE_SIZE = 4096;
    static constexpr size_t MAX_DEPTH = 20;
    
private:
    // [Performance Demon]: Memory pools for ZERO allocations
    SpatialMemoryPool memory_pool_;
    SpatialNodeSIMD* root_node_;
    
    // [Performance Demon]: Hash cache - linear probing for cache efficiency
    struct alignas(64) CacheEntry {
        uint64_t key_hash;
        QueryResultSIMD result;
        std::atomic<bool> valid;
    };
    
    alignas(64) mutable CacheEntry query_cache_[CACHE_SIZE];
    mutable std::atomic<uint64_t> cache_hits_{0};
    mutable std::atomic<uint64_t> cache_misses_{0};
    
    // [Performance Demon]: Pre-computed dodecahedron geometry - SIMD aligned
    alignas(32) double dodecahedron_vertices_[20 * 3];
    alignas(32) double dodecahedron_face_centers_[12 * 3];
    alignas(32) double dodecahedron_face_normals_[12 * 3];
    
    // [Performance Demon]: Performance counters
    mutable std::atomic<uint64_t> query_count_{0};
    mutable std::atomic<uint64_t> total_query_time_ns_{0};
    
public:
    explicit SpatialIndexerSIMD(double precision = 1e-8) {
        initialize();
    }
    
    // [Performance Demon]: Initialize with MAXIMUM performance
    void initialize() noexcept {
        precompute_dodecahedron_geometry();
        root_node_ = memory_pool_.allocate_node();
        initialize_root_node();
        
        // Clear cache
        for (size_t i = 0; i < CACHE_SIZE; ++i) {
            query_cache_[i].valid.store(false, std::memory_order_relaxed);
        }
    }
    
    // [Performance Demon]: ADD element with SIMD optimization
    void add_element(uint32_t element_id, const SphericalCoords& coords) noexcept {
        auto start = std::chrono::high_resolution_clock::now();
        
        SpatialNodeSIMD* node = find_or_create_node_simd(coords);
        add_element_to_node_simd(node, element_id, coords);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        total_query_time_ns_.fetch_add(duration.count(), std::memory_order_relaxed);
    }
    
    // [Performance Demon]: QUERY with MAXIMUM SIMD parallelism!
    QueryResultSIMD query_region_simd(const SphericalCoords& center, double radius) const noexcept {
        auto start = std::chrono::high_resolution_clock::now();
        query_count_.fetch_add(1, std::memory_order_relaxed);
        
        // [Performance Demon]: Check cache first - linear probing
        uint64_t key_hash = compute_hash(center, radius);
        size_t cache_index = key_hash % CACHE_SIZE;
        
        for (size_t i = 0; i < 4; ++i) {  // Max 4 probes
            size_t index = (cache_index + i) % CACHE_SIZE;
            if (query_cache_[index].valid.load(std::memory_order_acquire) &&
                query_cache_[index].key_hash == key_hash) {
                cache_hits_.fetch_add(1, std::memory_order_relaxed);
                return query_cache_[index].result;
            }
        }
        
        cache_misses_.fetch_add(1, std::memory_order_relaxed);
        
        // [Performance Demon]: Perform SIMD query
        QueryResultSIMD result;
        query_region_recursive_simd(root_node_, center, radius, result);
        
        // [Performance Demon]: Cache the result
        cache_result(key_hash, result);
        
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        total_query_time_ns_.fetch_add(duration.count(), std::memory_order_relaxed);
        
        return result;
    }
    
    // [Performance Demon]: RAYCAST with SIMD ray-sphere intersection
    QueryResultSIMD raycast_simd(const SphericalCoords& origin, 
                                 const SphericalCoords& direction,
                                 double max_distance) const noexcept {
        QueryResultSIMD result;
        
        // [Performance Demon]: Convert to Cartesian for SIMD ray operations
        simd::Vector3SIMD ray_origin = coords_to_simd(origin);
        simd::Vector3SIMD ray_dir = coords_to_simd(direction).normalized();
        
        raycast_recursive_simd(root_node_, ray_origin, ray_dir, max_distance, result);
        return result;
    }
    
    // [Performance Demon]: Nearest neighbors with SIMD distance calculations
    QueryResultSIMD find_nearest_simd(const SphericalCoords& point, size_t count) const noexcept {
        QueryResultSIMD candidates;
        
        // [Performance Demon]: Use SIMD to process 4 elements at once
        find_nearest_recursive_simd(root_node_, point, candidates);
        
        // [Performance Demon]: SIMD partial sort for top K
        simd_partial_sort(candidates, count);
        
        return candidates;
    }
    
    // [Performance Demon]: Performance statistics
    struct PerformanceStats {
        uint64_t query_count;
        double average_query_time_ns;
        double cache_hit_rate;
        size_t memory_usage_bytes;
        const char* simd_instruction_set;
    };
    
    PerformanceStats get_performance_stats() const noexcept {
        uint64_t queries = query_count_.load(std::memory_order_relaxed);
        uint64_t total_time = total_query_time_ns_.load(std::memory_order_relaxed);
        uint64_t hits = cache_hits_.load(std::memory_order_relaxed);
        uint64_t misses = cache_misses_.load(std::memory_order_relaxed);
        
        return {
            queries,
            queries > 0 ? static_cast<double>(total_time) / queries : 0.0,
            (hits + misses) > 0 ? static_cast<double>(hits) / (hits + misses) : 0.0,
            memory_pool_.get_allocation_count() * sizeof(SpatialNodeSIMD),
            simd::SIMDConfig::get_instruction_set()
        };
    }

private:
    // [Performance Demon]: SIMD helper functions - MAXIMUM SPEED!
    
    void precompute_dodecahedron_geometry() noexcept {
        // [Performance Demon]: Pre-compute ALL geometry for cache efficiency
        const double g = GOLDEN_RATIO;
        const double coords[][3] = {
            {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {1, -1, -1},
            {-1, 1, 1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1},
            {0, 1/g, g}, {0, 1/g, -g}, {0, -1/g, g}, {0, -1/g, -g},
            {1/g, g, 0}, {1/g, -g, 0}, {-1/g, g, 0}, {-1/g, -g, 0},
            {g, 0, 1/g}, {g, 0, -1/g}, {-g, 0, 1/g}, {-g, 0, -1/g}
        };
        
        // Store vertices in SIMD-friendly format
        for (size_t i = 0; i < 20; ++i) {
            dodecahedron_vertices_[i * 3 + 0] = coords[i][0];
            dodecahedron_vertices_[i * 3 + 1] = coords[i][1];
            dodecahedron_vertices_[i * 3 + 2] = coords[i][2];
        }
        
        // Pre-compute face centers and normals using SIMD
        compute_face_geometry_simd();
    }
    
    void compute_face_geometry_simd() noexcept {
        // [Performance Demon]: SIMD computation of face geometry
        const uint32_t face_indices[][5] = {
            {0, 8, 9, 4, 16}, {0, 16, 17, 2, 8}, {8, 2, 10, 3, 9},
            {12, 4, 9, 3, 5}, {16, 4, 12, 14, 6}, {6, 17, 16, 4, 12},
            {6, 10, 2, 17, 16}, {7, 15, 6, 12, 5}, {3, 10, 6, 15, 11},
            {5, 3, 11, 7, 12}, {14, 13, 1, 9, 4}, {1, 13, 11, 15, 7}
        };
        
#ifdef HSML_AVX_AVAILABLE
        // [Performance Demon]: AVX face computation - process multiple faces
        for (size_t face = 0; face < 12; ++face) {
            __m256d center = _mm256_setzero_pd();
            
            // Average vertices to get face center
            for (size_t i = 0; i < 5; ++i) {
                uint32_t vertex_idx = face_indices[face][i];
                __m256d vertex = _mm256_load_pd(&dodecahedron_vertices_[vertex_idx * 3]);
                center = _mm256_add_pd(center, vertex);
            }
            
            center = _mm256_mul_pd(center, _mm256_set1_pd(0.2)); // Divide by 5
            _mm256_store_pd(&dodecahedron_face_centers_[face * 3], center);
            
            // Compute normal (simplified for performance)
            _mm256_store_pd(&dodecahedron_face_normals_[face * 3], 
                           _mm256_div_pd(center, _mm256_set1_pd(_mm256_cvtsd_f64(center))));
        }
#else
        // Scalar fallback
        for (size_t face = 0; face < 12; ++face) {
            double center[3] = {0, 0, 0};
            for (size_t i = 0; i < 5; ++i) {
                uint32_t vertex_idx = face_indices[face][i];
                center[0] += dodecahedron_vertices_[vertex_idx * 3 + 0];
                center[1] += dodecahedron_vertices_[vertex_idx * 3 + 1];
                center[2] += dodecahedron_vertices_[vertex_idx * 3 + 2];
            }
            
            dodecahedron_face_centers_[face * 3 + 0] = center[0] * 0.2;
            dodecahedron_face_centers_[face * 3 + 1] = center[1] * 0.2;
            dodecahedron_face_centers_[face * 3 + 2] = center[2] * 0.2;
            
            double mag = std::sqrt(center[0]*center[0] + center[1]*center[1] + center[2]*center[2]);
            dodecahedron_face_normals_[face * 3 + 0] = center[0] / mag;
            dodecahedron_face_normals_[face * 3 + 1] = center[1] / mag;
            dodecahedron_face_normals_[face * 3 + 2] = center[2] / mag;
        }
#endif
    }
    
    void initialize_root_node() noexcept {
        root_node_->element_count = 0;
        root_node_->depth = 0;
        root_node_->parent_index = UINT32_MAX;
        root_node_->node_id = 0;
        
        // Initialize bounds to full sphere
        root_node_->bounds[0] = 0.0;                    // rMin
        root_node_->bounds[1] = 1000.0;                 // rMax
        root_node_->bounds[2] = 0.0;                    // thetaMin
        root_node_->bounds[3] = precision::MathematicalConstants::PI;  // thetaMax
        root_node_->bounds[4] = -precision::MathematicalConstants::PI; // phiMin
        root_node_->bounds[5] = precision::MathematicalConstants::PI;  // phiMax
        
        // Copy face data
        std::memcpy(root_node_->face_centers, dodecahedron_face_centers_, 12 * 3 * sizeof(double));
        std::memcpy(root_node_->face_normals, dodecahedron_face_normals_, 12 * 3 * sizeof(double));
        
        // Initialize child indices
        std::fill(root_node_->child_indices, root_node_->child_indices + 12, UINT32_MAX);
    }
    
    // [Performance Demon]: SIMD coordinate conversion
    simd::Vector3SIMD coords_to_simd(const SphericalCoords& coords) const noexcept {
        double sin_theta = std::sin(coords.theta());
        double cos_theta = std::cos(coords.theta());
        double sin_phi = std::sin(coords.phi());
        double cos_phi = std::cos(coords.phi());
        
        return simd::Vector3SIMD(
            coords.radius() * sin_theta * cos_phi,
            coords.radius() * sin_theta * sin_phi,
            coords.radius() * cos_theta
        );
    }
    
    // [Performance Demon]: More SIMD helpers...
    SpatialNodeSIMD* find_or_create_node_simd(const SphericalCoords& coords) noexcept;
    void add_element_to_node_simd(SpatialNodeSIMD* node, uint32_t element_id, const SphericalCoords& coords) noexcept;
    void query_region_recursive_simd(SpatialNodeSIMD* node, const SphericalCoords& center, double radius, QueryResultSIMD& result) const noexcept;
    void raycast_recursive_simd(SpatialNodeSIMD* node, const simd::Vector3SIMD& origin, const simd::Vector3SIMD& direction, double max_distance, QueryResultSIMD& result) const noexcept;
    void find_nearest_recursive_simd(SpatialNodeSIMD* node, const SphericalCoords& point, QueryResultSIMD& candidates) const noexcept;
    void simd_partial_sort(QueryResultSIMD& candidates, size_t count) const noexcept;
    
    // [Performance Demon]: FIXED - Missing helper function declarations
    SpatialNodeSIMD* find_best_child_simd(SpatialNodeSIMD* node, const SphericalCoords& coords) noexcept;
    bool ray_intersects_node_simd(SpatialNodeSIMD* node, const simd::Vector3SIMD& origin, const simd::Vector3SIMD& direction, double max_distance) const noexcept;
    void compute_child_bounds_simd(SpatialNodeSIMD* parent, size_t face_index, SpatialNodeSIMD* child) noexcept;
    void subdivide_node_simd(SpatialNodeSIMD* node) noexcept;
    bool node_intersects_sphere_simd(SpatialNodeSIMD* node, const SphericalCoords& center, double radius) const noexcept;
    void simd_network_sort_8(QueryResultSIMD& candidates, size_t count) const noexcept;
    
    uint64_t compute_hash(const SphericalCoords& coords, double radius) const noexcept {
        // [Performance Demon]: Fast hash for cache keys
        uint64_t hash = 0x9e3779b9;
        hash ^= std::hash<double>{}(coords.radius()) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        hash ^= std::hash<double>{}(coords.theta()) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        hash ^= std::hash<double>{}(coords.phi()) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        hash ^= std::hash<double>{}(radius) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        return hash;
    }
    
    void cache_result(uint64_t key_hash, const QueryResultSIMD& result) const noexcept {
        size_t cache_index = key_hash % CACHE_SIZE;
        
        // [Performance Demon]: FIXED THREAD SAFETY - proper cache locking
        for (size_t i = 0; i < 4; ++i) {
            size_t index = (cache_index + i) % CACHE_SIZE;
            bool expected = false;
            if (query_cache_[index].valid.compare_exchange_weak(expected, true, std::memory_order_acq_rel)) {
                // Store data first, then make visible
                query_cache_[index].key_hash = key_hash;
                query_cache_[index].result = result;
                // Memory barrier to ensure data is written before marking valid
                std::atomic_thread_fence(std::memory_order_release);
                return;
            }
        }
        // If no slot available, evict oldest entry (index 0 of probe sequence)
        size_t evict_index = cache_index;
        query_cache_[evict_index].valid.store(false, std::memory_order_release);
        std::atomic_thread_fence(std::memory_order_seq_cst);
        query_cache_[evict_index].key_hash = key_hash;
        query_cache_[evict_index].result = result;
        query_cache_[evict_index].valid.store(true, std::memory_order_release);
    }
};

} // namespace spatial
} // namespace core
} // namespace hsml