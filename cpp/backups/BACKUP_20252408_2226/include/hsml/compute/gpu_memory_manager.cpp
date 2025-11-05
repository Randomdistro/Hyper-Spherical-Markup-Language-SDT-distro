#pragma once

#include "../core/advanced_concepts.h"
#include <atomic>
#include <memory>
#include <span>
#include <vector>
#include <unordered_map>
#include <mutex>
#include <condition_variable>
#include <thread>
#include <future>
#include <cstring>

namespace hsml {
namespace compute {

using namespace core::concepts;

// Memory usage types
enum class memory_usage {
    static_read,      // GPU read-only, CPU write-once
    dynamic_read,     // GPU read-only, CPU write frequently  
    static_write,     // GPU write-only, CPU read-once
    dynamic_write,    // GPU write-only, CPU read frequently
    read_write,       // GPU and CPU both read/write
    staging,          // Temporary staging buffer
    compute_shader    // Compute shader scratch space
};

// Memory synchronization states
enum class sync_state : uint32_t {
    cpu_dirty = 0,     // CPU has newer data
    gpu_dirty = 1,     // GPU has newer data
    synchronized = 2,  // CPU and GPU in sync
    invalid = 3        // Buffer is invalid
};

// Buffer allocation metadata
template<typename T>
struct gpu_allocation {
    void* gpu_buffer = nullptr;        // GPU buffer handle (API-specific)
    std::byte* mapped_ptr = nullptr;   // CPU-mapped memory pointer
    size_t size_bytes = 0;             // Total size in bytes
    size_t element_count = 0;          // Number of T elements
    memory_usage usage = memory_usage::static_read;
    std::atomic<sync_state> sync_state{sync_state::synchronized};
    std::atomic<uint64_t> cpu_timestamp{0};
    std::atomic<uint64_t> gpu_timestamp{0};
    std::atomic<uint32_t> ref_count{1};
    bool persistent_mapping = false;
    
    gpu_allocation() = default;
    gpu_allocation(const gpu_allocation&) = delete;
    gpu_allocation& operator=(const gpu_allocation&) = delete;
    
    gpu_allocation(gpu_allocation&& other) noexcept {
        *this = std::move(other);
    }
    
    gpu_allocation& operator=(gpu_allocation&& other) noexcept {
        if (this != &other) {
            gpu_buffer = std::exchange(other.gpu_buffer, nullptr);
            mapped_ptr = std::exchange(other.mapped_ptr, nullptr);
            size_bytes = std::exchange(other.size_bytes, 0);
            element_count = std::exchange(other.element_count, 0);
            usage = other.usage;
            sync_state.store(other.sync_state.load());
            cpu_timestamp.store(other.cpu_timestamp.load());
            gpu_timestamp.store(other.gpu_timestamp.load());
            ref_count.store(other.ref_count.load());
            persistent_mapping = other.persistent_mapping;
        }
        return *this;
    }
};

// GPU buffer view with automatic synchronization
template<typename T>
class gpu_buffer_view {
    static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable for GPU transfer");
    
    std::shared_ptr<gpu_allocation<T>> allocation_;
    
public:
    using value_type = T;
    using pointer = T*;
    using const_pointer = const T*;
    using reference = T&;
    using const_reference = const T&;
    using size_type = size_t;
    
    gpu_buffer_view() = default;
    
    explicit gpu_buffer_view(std::shared_ptr<gpu_allocation<T>> alloc)
        : allocation_(std::move(alloc)) {}
    
    // Element access with automatic synchronization
    [[nodiscard]] T* data() noexcept {
        ensure_cpu_sync();
        return reinterpret_cast<T*>(allocation_->mapped_ptr);
    }
    
    [[nodiscard]] const T* data() const noexcept {
        ensure_cpu_sync();
        return reinterpret_cast<const T*>(allocation_->mapped_ptr);
    }
    
    [[nodiscard]] size_t size() const noexcept {
        return allocation_ ? allocation_->element_count : 0;
    }
    
    [[nodiscard]] size_t size_bytes() const noexcept {
        return allocation_ ? allocation_->size_bytes : 0;
    }
    
    [[nodiscard]] bool empty() const noexcept {
        return size() == 0;
    }
    
    // Iterator interface
    [[nodiscard]] T* begin() noexcept { return data(); }
    [[nodiscard]] T* end() noexcept { return data() + size(); }
    [[nodiscard]] const T* begin() const noexcept { return data(); }
    [[nodiscard]] const T* end() const noexcept { return data() + size(); }
    [[nodiscard]] const T* cbegin() const noexcept { return data(); }
    [[nodiscard]] const T* cend() const noexcept { return data() + size(); }
    
    // Element access
    [[nodiscard]] reference operator[](size_type index) noexcept {
        ensure_cpu_sync();
        return data()[index];
    }
    
    [[nodiscard]] const_reference operator[](size_type index) const noexcept {
        ensure_cpu_sync();
        return data()[index];
    }
    
    // Synchronization control
    void sync_to_gpu() {
        if (!allocation_) return;
        
        auto expected = sync_state::cpu_dirty;
        if (allocation_->sync_state.compare_exchange_strong(expected, sync_state::synchronized)) {
            // Perform actual GPU sync (implementation-specific)
            flush_to_gpu();
            allocation_->gpu_timestamp.store(get_timestamp());
        }
    }
    
    void sync_from_gpu() {
        if (!allocation_) return;
        
        auto expected = sync_state::gpu_dirty;
        if (allocation_->sync_state.compare_exchange_strong(expected, sync_state::synchronized)) {
            // Perform actual CPU sync (implementation-specific)
            invalidate_cpu_cache();
            allocation_->cpu_timestamp.store(get_timestamp());
        }
    }
    
    void mark_cpu_dirty() {
        if (allocation_) {
            allocation_->sync_state.store(sync_state::cpu_dirty);
            allocation_->cpu_timestamp.store(get_timestamp());
        }
    }
    
    void mark_gpu_dirty() {
        if (allocation_) {
            allocation_->sync_state.store(sync_state::gpu_dirty);
            allocation_->gpu_timestamp.store(get_timestamp());
        }
    }
    
    [[nodiscard]] bool is_synchronized() const noexcept {
        return allocation_ && allocation_->sync_state.load() == sync_state::synchronized;
    }
    
    [[nodiscard]] memory_usage get_usage() const noexcept {
        return allocation_ ? allocation_->usage : memory_usage::static_read;
    }
    
    // GPU buffer handle (API-specific)
    [[nodiscard]] void* gpu_handle() const noexcept {
        return allocation_ ? allocation_->gpu_buffer : nullptr;
    }
    
private:
    void ensure_cpu_sync() const {
        if (!allocation_) return;
        
        if (allocation_->sync_state.load() == sync_state::gpu_dirty) {
            const_cast<gpu_buffer_view*>(this)->sync_from_gpu();
        }
    }
    
    void flush_to_gpu() {
        // Implementation would call actual GPU API to flush cache
        // For demonstration, this is a no-op
    }
    
    void invalidate_cpu_cache() {
        // Implementation would invalidate CPU cache for memory-mapped region
        // For demonstration, this is a no-op
    }
    
    [[nodiscard]] static uint64_t get_timestamp() noexcept {
        return std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
};

// Compute mapping for direct memory access
template<typename T>
struct compute_mapping {
    T* data;
    size_t size;
    std::function<void()> sync;
    
    compute_mapping(T* ptr, size_t count, std::function<void()> sync_func)
        : data(ptr), size(count), sync(std::move(sync_func)) {}
    
    ~compute_mapping() {
        if (sync) sync();
    }
    
    compute_mapping(const compute_mapping&) = delete;
    compute_mapping& operator=(const compute_mapping&) = delete;
    
    compute_mapping(compute_mapping&& other) noexcept
        : data(std::exchange(other.data, nullptr))
        , size(std::exchange(other.size, 0))
        , sync(std::move(other.sync)) {}
    
    compute_mapping& operator=(compute_mapping&& other) noexcept {
        if (this != &other) {
            if (sync) sync();
            data = std::exchange(other.data, nullptr);
            size = std::exchange(other.size, 0);
            sync = std::move(other.sync);
        }
        return *this;
    }
    
    [[nodiscard]] T* begin() noexcept { return data; }
    [[nodiscard]] T* end() noexcept { return data + size; }
    [[nodiscard]] const T* begin() const noexcept { return data; }
    [[nodiscard]] const T* end() const noexcept { return data + size; }
};

// Memory-mapped GPU compute pipeline
template<graphics_api API = graphics_api::vulkan>
class gpu_memory_manager {
private:
    static constexpr size_t DEFAULT_POOL_SIZE = 256 * 1024 * 1024;  // 256MB
    static constexpr size_t ALIGNMENT = 256;  // GPU alignment requirement
    
    struct memory_pool {
        std::byte* base_ptr = nullptr;
        size_t total_size = 0;
        size_t used_size = 0;
        std::atomic<bool> in_use{false};
        std::vector<std::pair<size_t, size_t>> free_blocks;  // offset, size pairs
        mutable std::mutex mutex;
    };
    
    std::vector<std::unique_ptr<memory_pool>> pools_;
    std::unordered_map<void*, std::weak_ptr<gpu_allocation<std::byte>>> active_allocations_;
    mutable std::shared_mutex allocations_mutex_;
    std::atomic<uint64_t> allocation_counter_{0};
    
    // API-specific buffer types
    using buffer_type = std::conditional_t<
        API == graphics_api::vulkan, uintptr_t,  // VkBuffer
        std::conditional_t<
            API == graphics_api::opengl, uint32_t,  // GLuint
            void*  // Generic handle
        >
    >;
    
    [[nodiscard]] size_t align_size(size_t size) const noexcept {
        return (size + ALIGNMENT - 1) & ~(ALIGNMENT - 1);
    }
    
    [[nodiscard]] std::unique_ptr<memory_pool> create_pool(size_t size) {
        auto pool = std::make_unique<memory_pool>();
        pool->total_size = align_size(size);
        
        if constexpr (API == graphics_api::vulkan) {
            pool->base_ptr = create_vulkan_memory_pool(pool->total_size);
        } else if constexpr (API == graphics_api::opengl) {
            pool->base_ptr = create_opengl_memory_pool(pool->total_size);
        } else {
            // Fallback to system memory
            pool->base_ptr = static_cast<std::byte*>(std::aligned_alloc(ALIGNMENT, pool->total_size));
        }
        
        if (pool->base_ptr) {
            pool->free_blocks.emplace_back(0, pool->total_size);
        }
        
        return pool;
    }
    
    [[nodiscard]] std::byte* create_vulkan_memory_pool(size_t size) {
        // Implementation would create VkDeviceMemory and map it
        // For demonstration, fall back to system memory
        return static_cast<std::byte*>(std::aligned_alloc(ALIGNMENT, size));
    }
    
    [[nodiscard]] std::byte* create_opengl_memory_pool(size_t size) {
        // Implementation would create persistent mapped OpenGL buffer
        // For demonstration, fall back to system memory
        return static_cast<std::byte*>(std::aligned_alloc(ALIGNMENT, size));
    }
    
    [[nodiscard]] std::pair<std::byte*, std::unique_ptr<memory_pool>*> 
    allocate_from_pool(size_t size) {
        size = align_size(size);
        
        // Try existing pools first
        for (auto& pool : pools_) {
            std::lock_guard<std::mutex> lock(pool->mutex);
            
            for (auto it = pool->free_blocks.begin(); it != pool->free_blocks.end(); ++it) {
                if (it->second >= size) {
                    size_t offset = it->first;
                    std::byte* ptr = pool->base_ptr + offset;
                    
                    // Split the block if necessary
                    if (it->second > size) {
                        it->first += size;
                        it->second -= size;
                    } else {
                        pool->free_blocks.erase(it);
                    }
                    
                    pool->used_size += size;
                    return {ptr, &pool};
                }
            }
        }
        
        // Create new pool if needed
        size_t pool_size = std::max(size * 2, DEFAULT_POOL_SIZE);
        auto new_pool = create_pool(pool_size);
        if (!new_pool || !new_pool->base_ptr) {
            return {nullptr, nullptr};
        }
        
        std::byte* ptr = new_pool->base_ptr;
        new_pool->used_size = size;
        new_pool->free_blocks[0] = {size, pool_size - size};
        
        pools_.push_back(std::move(new_pool));
        return {ptr, &pools_.back()};
    }
    
    template<typename T>
    buffer_type create_gpu_buffer(size_t size_bytes) {
        if constexpr (API == graphics_api::vulkan) {
            return create_vulkan_buffer(size_bytes);
        } else if constexpr (API == graphics_api::opengl) {
            return create_opengl_buffer(size_bytes);
        } else {
            return static_cast<buffer_type>(size_bytes);  // Dummy handle
        }
    }
    
    buffer_type create_vulkan_buffer(size_t size) {
        // Implementation would create VkBuffer
        // For demonstration, return placeholder
        return static_cast<buffer_type>(allocation_counter_.fetch_add(1) + 1);
    }
    
    buffer_type create_opengl_buffer(size_t size) {
        // Implementation would call glCreateBuffers/glBufferStorage
        // For demonstration, return placeholder
        return static_cast<buffer_type>(allocation_counter_.fetch_add(1) + 1);
    }
    
public:
    gpu_memory_manager() {
        // Pre-allocate initial pool
        pools_.push_back(create_pool(DEFAULT_POOL_SIZE));
    }
    
    ~gpu_memory_manager() {
        // Clean up all pools
        for (auto& pool : pools_) {
            if (pool && pool->base_ptr) {
                if constexpr (API == graphics_api::vulkan || API == graphics_api::opengl) {
                    // Would call API-specific cleanup
                } else {
                    std::free(pool->base_ptr);
                }
            }
        }
    }
    
    template<typename T>
    [[nodiscard]] gpu_buffer_view<T> allocate_buffer(size_t count, memory_usage usage) {
        static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable");
        
        size_t size_bytes = count * sizeof(T);
        auto [mapped_ptr, pool_ptr] = allocate_from_pool(size_bytes);
        
        if (!mapped_ptr) {
            return gpu_buffer_view<T>{};  // Failed allocation
        }
        
        auto allocation = std::make_shared<gpu_allocation<T>>();
        allocation->gpu_buffer = reinterpret_cast<void*>(create_gpu_buffer<T>(size_bytes));
        allocation->mapped_ptr = mapped_ptr;
        allocation->size_bytes = size_bytes;
        allocation->element_count = count;
        allocation->usage = usage;
        allocation->persistent_mapping = true;
        
        // Register allocation for tracking
        {
            std::unique_lock<std::shared_mutex> lock(allocations_mutex_);
            active_allocations_[allocation->gpu_buffer] = allocation;
        }
        
        return gpu_buffer_view<T>{std::move(allocation)};
    }
    
    template<typename T>
    [[nodiscard]] compute_mapping<T> map_for_compute(gpu_buffer_view<T> buffer) {
        if (buffer.empty()) {
            return compute_mapping<T>{nullptr, 0, []{}};
        }
        
        T* data_ptr = buffer.data();
        size_t count = buffer.size();
        
        // Create sync function for automatic cache management
        auto sync_func = [buffer]() mutable {
            buffer.sync_to_gpu();
        };
        
        return compute_mapping<T>{data_ptr, count, std::move(sync_func)};
    }
    
    // Batch allocation for performance
    template<typename T>
    [[nodiscard]] std::vector<gpu_buffer_view<T>> 
    allocate_buffers_batch(std::span<const size_t> counts, memory_usage usage) {
        
        std::vector<gpu_buffer_view<T>> buffers;
        buffers.reserve(counts.size());
        
        for (size_t count : counts) {
            buffers.push_back(allocate_buffer<T>(count, usage));
        }
        
        return buffers;
    }
    
    // Memory statistics
    [[nodiscard]] size_t get_total_allocated() const noexcept {
        size_t total = 0;
        for (const auto& pool : pools_) {
            if (pool) {
                std::lock_guard<std::mutex> lock(pool->mutex);
                total += pool->used_size;
            }
        }
        return total;
    }
    
    [[nodiscard]] size_t get_total_capacity() const noexcept {
        size_t total = 0;
        for (const auto& pool : pools_) {
            if (pool) {
                total += pool->total_size;
            }
        }
        return total;
    }
    
    [[nodiscard]] double get_utilization() const noexcept {
        size_t allocated = get_total_allocated();
        size_t capacity = get_total_capacity();
        return capacity > 0 ? static_cast<double>(allocated) / capacity : 0.0;
    }
    
    [[nodiscard]] size_t get_active_allocations() const noexcept {
        std::shared_lock<std::shared_mutex> lock(allocations_mutex_);
        return active_allocations_.size();
    }
    
    // Force garbage collection of unused allocations
    void garbage_collect() {
        std::unique_lock<std::shared_mutex> lock(allocations_mutex_);
        
        auto it = active_allocations_.begin();
        while (it != active_allocations_.end()) {
            if (it->second.expired()) {
                it = active_allocations_.erase(it);
            } else {
                ++it;
            }
        }
    }
    
    // Flush all pending GPU operations
    void flush_all() {
        std::shared_lock<std::shared_mutex> lock(allocations_mutex_);
        
        for (const auto& [handle, weak_alloc] : active_allocations_) {
            if (auto allocation = weak_alloc.lock()) {
                if (allocation->sync_state.load() == sync_state::cpu_dirty) {
                    // Force GPU sync
                    allocation->sync_state.store(sync_state::synchronized);
                }
            }
        }
    }
};

// Type aliases for different graphics APIs
using vulkan_memory_manager = gpu_memory_manager<graphics_api::vulkan>;
using opengl_memory_manager = gpu_memory_manager<graphics_api::opengl>;
using generic_memory_manager = gpu_memory_manager<graphics_api::software>;

// Concept verification
static_assert(gpu_buffer_concept<gpu_buffer_view<float>>);
static_assert(gpu_buffer_concept<gpu_buffer_view<double>>);

} // namespace compute
} // namespace hsml