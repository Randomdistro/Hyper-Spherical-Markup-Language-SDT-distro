/** @file IMemoryManager.h
 * @brief Cross-platform memory management abstraction
 *
 * Clean Architecture: Infrastructure Layer - Memory Management
 * Provides unified memory operations with performance monitoring.
 */

#pragma once

#include <string>
#include <memory>
#include <vector>
#include <functional>
#include <chrono>

namespace hsml {
namespace infrastructure {

// Memory allocation types
enum class MemoryType {
    GENERAL_PURPOSE,
    SPHERICAL_COORDINATES,
    GRAPHICS_DATA,
    AUDIO_DATA,
    NETWORK_DATA,
    FILE_DATA,
    TEMPORARY,
    PERSISTENT
};

// Memory allocation strategies
enum class AllocationStrategy {
    FIRST_FIT,
    BEST_FIT,
    WORST_FIT,
    BUDDY_SYSTEM,
    SLAB_ALLOCATION,
    POOL_ALLOCATION
};

// Memory region information
struct MemoryRegion {
    void* address{nullptr};
    size_t size{0};
    MemoryType type{MemoryType::GENERAL_PURPOSE};
    std::string name;
    std::chrono::steady_clock::time_point allocation_time;
    uint64_t allocation_id{0};
    bool is_free{false};
    size_t alignment{0};
};

// Memory statistics
struct MemoryStatistics {
    uint64_t total_allocated_bytes{0};
    uint64_t peak_allocated_bytes{0};
    uint64_t current_allocated_bytes{0};
    uint64_t total_free_bytes{0};
    uint32_t allocation_count{0};
    uint32_t deallocation_count{0};
    uint32_t live_allocation_count{0};
    double average_allocation_size{0.0};
    double average_allocation_time_us{0.0};
    uint32_t fragmentation_ratio{0}; // 0-100
    uint32_t cache_miss_ratio{0};    // 0-100
    uint64_t page_fault_count{0};
};

// Memory pool configuration
struct MemoryPoolConfig {
    size_t block_size{64};
    size_t block_count{1024};
    size_t alignment{16};
    bool thread_safe{true};
    MemoryType memory_type{MemoryType::GENERAL_PURPOSE};
};

// Memory leak detection
struct MemoryLeakInfo {
    void* address{nullptr};
    size_t size{0};
    std::string allocation_location;
    std::chrono::steady_clock::time_point allocation_time;
    uint64_t allocation_id{0};
    bool is_leaked{false};
};

/**
 * @brief Memory pool interface for efficient allocations
 */
class IMemoryPool {
public:
    virtual ~IMemoryPool() = default;

    virtual bool initialize(const MemoryPoolConfig& config) = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;

    virtual void* allocate() = 0;
    virtual void deallocate(void* ptr) = 0;
    virtual bool owns(void* ptr) const = 0;

    virtual size_t get_block_size() const = 0;
    virtual size_t get_total_blocks() const = 0;
    virtual size_t get_used_blocks() const = 0;
    virtual size_t get_free_blocks() const = 0;
    virtual double get_utilization_ratio() const = 0;
};

/**
 * @brief Memory arena interface for bulk allocations
 */
class IMemoryArena {
public:
    virtual ~IMemoryArena() = default;

    virtual bool initialize(size_t size, size_t alignment = 16) = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;

    virtual void* allocate(size_t size, size_t alignment = 0) = 0;
    virtual void* reallocate(void* ptr, size_t new_size) = 0;
    virtual void deallocate(void* ptr) = 0;
    virtual void reset() = 0;

    virtual size_t get_size() const = 0;
    virtual size_t get_used_size() const = 0;
    virtual size_t get_free_size() const = 0;
    virtual double get_utilization_ratio() const = 0;

    virtual bool contains(void* ptr) const = 0;
    virtual void* get_base_address() const = 0;
};

/**
 * @brief Memory profiler interface for performance monitoring
 */
class IMemoryProfiler {
public:
    virtual ~IMemoryProfiler() = default;

    virtual bool start_profiling() = 0;
    virtual void stop_profiling() = 0;
    virtual bool is_profiling() const = 0;

    virtual const MemoryStatistics& get_statistics() const = 0;
    virtual std::vector<MemoryRegion> get_memory_map() const = 0;
    virtual std::vector<MemoryLeakInfo> get_leak_report() const = 0;

    virtual void reset_statistics() = 0;
    virtual void take_snapshot(const std::string& name) = 0;
    virtual bool compare_snapshots(const std::string& snapshot1,
                                 const std::string& snapshot2) const = 0;
};

/**
 * @brief Cross-platform memory manager interface
 *
 * This interface provides unified memory management operations across platforms
 * with specialized support for spherical coordinate data and performance monitoring.
 */
class IMemoryManager {
public:
    virtual ~IMemoryManager() = default;

    // Initialization
    virtual bool initialize() = 0;
    virtual void shutdown() = 0;
    virtual bool is_initialized() const = 0;

    // Basic memory operations
    virtual void* allocate(size_t size, size_t alignment = 0,
                          MemoryType type = MemoryType::GENERAL_PURPOSE,
                          const std::string& name = "") = 0;
    virtual void* reallocate(void* ptr, size_t new_size, size_t alignment = 0) = 0;
    virtual void deallocate(void* ptr) = 0;
    virtual size_t get_allocation_size(void* ptr) const = 0;

    // Memory pools for performance
    virtual std::unique_ptr<IMemoryPool> create_pool(const MemoryPoolConfig& config) = 0;
    virtual std::vector<IMemoryPool*> get_all_pools() const = 0;

    // Memory arenas for bulk operations
    virtual std::unique_ptr<IMemoryArena> create_arena(size_t size, size_t alignment = 16) = 0;
    virtual std::vector<IMemoryArena*> get_all_arenas() const = 0;

    // Spherical coordinate specific operations
    virtual void* allocate_spherical_coords(size_t count, const std::string& name = "") = 0;
    virtual void* allocate_spherical_mesh(size_t vertex_count, size_t index_count = 0,
                                        const std::string& name = "") = 0;
    virtual void* allocate_spherical_texture(size_t width, size_t height, size_t channels = 4,
                                           const std::string& name = "") = 0;

    // Memory utilities
    virtual void zero_memory(void* ptr, size_t size) = 0;
    virtual void copy_memory(void* dst, const void* src, size_t size) = 0;
    virtual void move_memory(void* dst, const void* src, size_t size) = 0;
    virtual int compare_memory(const void* ptr1, const void* ptr2, size_t size) const = 0;

    // Memory alignment utilities
    virtual size_t get_alignment_padding(void* ptr, size_t alignment) const = 0;
    virtual void* align_pointer(void* ptr, size_t alignment) const = 0;
    virtual bool is_aligned(void* ptr, size_t alignment) const = 0;

    // Memory protection
    virtual bool protect_memory(void* ptr, size_t size, bool read, bool write, bool execute) = 0;
    virtual bool unprotect_memory(void* ptr, size_t size) = 0;

    // Memory mapping (for large files or shared memory)
    virtual void* map_memory(size_t size, bool read, bool write, bool execute) = 0;
    virtual bool unmap_memory(void* ptr, size_t size) = 0;
    virtual bool flush_memory(void* ptr, size_t size) = 0;

    // Memory profiler
    virtual IMemoryProfiler* get_profiler() = 0;
    virtual std::unique_ptr<IMemoryProfiler> create_profiler() = 0;

    // Memory statistics and monitoring
    virtual const MemoryStatistics& get_statistics() const = 0;
    virtual std::vector<MemoryRegion> get_memory_regions() const = 0;
    virtual std::vector<MemoryLeakInfo> detect_leaks() const = 0;

    // Memory optimization
    virtual void defragment() = 0;
    virtual void compact() = 0;
    virtual void garbage_collect() = 0;

    // Platform-specific operations
    virtual size_t get_page_size() const = 0;
    virtual size_t get_cache_line_size() const = 0;
    virtual size_t get_large_page_size() const = 0;
    virtual bool supports_large_pages() const = 0;
    virtual bool allocate_large_pages(size_t size) = 0;

    // Performance monitoring
    virtual double get_allocation_rate() const = 0; // allocations per second
    virtual double get_deallocation_rate() const = 0; // deallocations per second
    virtual double get_fragmentation_ratio() const = 0; // 0.0 to 1.0
    virtual uint64_t get_cache_miss_count() const = 0;
    virtual uint64_t get_page_fault_count() const = 0;

    // Error handling
    virtual std::string get_last_error() const = 0;
    virtual void clear_error() = 0;

    // Platform-specific
    virtual void* get_platform_handle() const = 0;
};

/**
 * @brief Memory manager factory interface
 */
class IMemoryManagerFactory {
public:
    virtual ~IMemoryManagerFactory() = default;
    virtual std::unique_ptr<IMemoryManager> create_memory_manager() = 0;
};

/**
 * @brief Smart pointer types for memory management
 */
template<typename T>
class MemoryPtr {
public:
    MemoryPtr() = default;
    explicit MemoryPtr(T* ptr, IMemoryManager* manager = nullptr);
    MemoryPtr(const MemoryPtr&) = delete;
    MemoryPtr(MemoryPtr&& other) noexcept;
    ~MemoryPtr();

    MemoryPtr& operator=(const MemoryPtr&) = delete;
    MemoryPtr& operator=(MemoryPtr&& other) noexcept;

    T* get() const { return ptr_; }
    T* operator->() const { return ptr_; }
    T& operator*() const { return *ptr_; }
    explicit operator bool() const { return ptr_ != nullptr; }

    void reset(T* ptr = nullptr);
    T* release();

private:
    T* ptr_{nullptr};
    IMemoryManager* manager_{nullptr};
};

// Utility functions
template<typename T, typename... Args>
MemoryPtr<T> make_memory_ptr(IMemoryManager* manager, Args&&... args) {
    if (!manager) return MemoryPtr<T>();

    void* raw_ptr = manager->allocate(sizeof(T), alignof(T), MemoryType::GENERAL_PURPOSE);
    if (!raw_ptr) return MemoryPtr<T>();

    T* typed_ptr = new(raw_ptr) T(std::forward<Args>(args)...);
    return MemoryPtr<T>(typed_ptr, manager);
}

} // namespace infrastructure
} // namespace hsml
