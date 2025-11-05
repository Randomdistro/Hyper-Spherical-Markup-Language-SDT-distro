#pragma once

/* [The Hacktivist]: Memory-mapped files are PURE UNIX WIZARDRY! */
/* [Security Paranoid]: *hyperventilating* Direct memory access... so many vulnerabilities! */
/* [Enterprise Bean]: We need proper abstraction layers for this system-level code */

#include "../core/advanced_concepts.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <atomic>
#include <memory>
#include <string>
#include <optional>
#include <chrono>
#include <filesystem>
#include <bit>

namespace hsml {
namespace compute {

using namespace core::concepts;

// [Security Paranoid]: Paranoid validation for memory-mapped operations
enum class mmap_security_level {
    permissive,    // Basic validation
    strict,        // Enhanced bounds checking
    paranoid       // Full sanitization + checksums
};

// [Enterprise Bean]: Comprehensive error handling enumeration
enum class mmap_error_code {
    success = 0,
    file_not_found,
    permission_denied,
    insufficient_memory,
    invalid_size,
    mmap_failed,
    msync_failed,
    munmap_failed,
    checksum_mismatch,
    bounds_violation
};

// [Hacktivist]: Raw performance metrics because NUMBERS DON'T LIE
struct mmap_performance_metrics {
    std::atomic<uint64_t> page_faults{0};
    std::atomic<uint64_t> cache_hits{0};
    std::atomic<uint64_t> cache_misses{0};
    std::atomic<uint64_t> bytes_read{0};
    std::atomic<uint64_t> bytes_written{0};
    std::atomic<uint64_t> sync_operations{0};
    std::atomic<uint64_t> validation_failures{0};
};

// [Security Paranoid]: Memory-mapped cache entry with extensive validation
template<typename T>
struct alignas(4096) mmap_cache_entry {  // [Hacktivist]: Page-aligned for optimal mmap performance
    static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable for mmap storage");
    
    // [Security Paranoid]: Multiple checksums for integrity verification
    alignas(8) uint64_t entry_checksum;
    alignas(8) uint64_t key_hash;
    alignas(8) uint64_t creation_timestamp;
    alignas(8) uint64_t last_access_timestamp;
    alignas(8) uint32_t ttl_seconds;
    alignas(8) uint32_t access_count;
    alignas(8) uint32_t size_bytes;
    alignas(8) uint32_t magic_number;  // [Security Paranoid]: Corruption detection
    
    alignas(alignof(T)) T value;
    
    // [Security Paranoid]: Padding to prevent buffer overflows
    alignas(8) uint64_t end_sentinel;
    
    static constexpr uint32_t MAGIC_VALUE = 0xDEADBEEF;
    static constexpr uint64_t END_SENTINEL = 0xCAFEBABEDEADC0DE;
    
    mmap_cache_entry() : magic_number(MAGIC_VALUE), end_sentinel(END_SENTINEL) {}
    
    template<typename U>
    mmap_cache_entry(uint64_t hash, U&& val, uint32_t ttl) 
        : key_hash(hash)
        , creation_timestamp(get_nanoseconds())
        , last_access_timestamp(get_nanoseconds())
        , ttl_seconds(ttl)
        , access_count(1)
        , size_bytes(sizeof(T))
        , magic_number(MAGIC_VALUE)
        , value(std::forward<U>(val))
        , end_sentinel(END_SENTINEL) {
        update_checksum();
    }
    
    // [Security Paranoid]: Comprehensive validation
    [[nodiscard]] bool is_valid() const noexcept {
        return magic_number == MAGIC_VALUE && 
               end_sentinel == END_SENTINEL &&
               verify_checksum();
    }
    
    [[nodiscard]] bool is_expired() const noexcept {
        const auto now = get_nanoseconds();
        const auto ttl_ns = static_cast<uint64_t>(ttl_seconds) * 1'000'000'000ULL;
        return (now - creation_timestamp) > ttl_ns;
    }
    
    void update_access() noexcept {
        last_access_timestamp = get_nanoseconds();
        ++access_count;
        update_checksum();
    }
    
private:
    void update_checksum() noexcept {
        // [Security Paranoid]: CRC-like checksum calculation
        uint64_t checksum = key_hash;
        checksum ^= creation_timestamp;
        checksum ^= last_access_timestamp;
        checksum ^= ttl_seconds;
        checksum ^= access_count;
        checksum ^= size_bytes;
        
        // [Hacktivist]: XOR all value bytes for integrity
        const auto* bytes = reinterpret_cast<const uint8_t*>(&value);
        for (size_t i = 0; i < sizeof(T); ++i) {
            checksum ^= bytes[i];
        }
        
        entry_checksum = checksum;
    }
    
    [[nodiscard]] bool verify_checksum() const noexcept {
        uint64_t calculated = key_hash;
        calculated ^= creation_timestamp;
        calculated ^= last_access_timestamp;
        calculated ^= ttl_seconds;
        calculated ^= access_count;
        calculated ^= size_bytes;
        
        const auto* bytes = reinterpret_cast<const uint8_t*>(&value);
        for (size_t i = 0; i < sizeof(T); ++i) {
            calculated ^= bytes[i];
        }
        
        return calculated == entry_checksum;
    }
    
    static uint64_t get_nanoseconds() noexcept {
        return std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
};

// [Hacktivist]: RAII wrapper for memory-mapped file management
template<typename T>
class mmap_file_manager {
private:
    int fd_ = -1;
    void* mapped_memory_ = nullptr;
    size_t file_size_ = 0;
    std::string file_path_;
    mmap_security_level security_level_;
    
    // [Security Paranoid]: Track memory access patterns
    mutable mmap_performance_metrics metrics_;
    
public:
    using entry_type = mmap_cache_entry<T>;
    
    mmap_file_manager(const std::string& path, size_t max_entries, 
                     mmap_security_level security = mmap_security_level::strict)
        : file_path_(path), security_level_(security) {
        
        file_size_ = max_entries * sizeof(entry_type);
        
        // [Hacktivist]: Ensure directory exists because we're PRACTICAL
        std::filesystem::create_directories(std::filesystem::path(path).parent_path());
        
        // [Security Paranoid]: Secure file creation with proper permissions
        fd_ = open(path.c_str(), O_CREAT | O_RDWR, S_IRUSR | S_IWUSR);
        if (fd_ == -1) {
            throw std::runtime_error("Failed to open cache file: " + path);
        }
        
        // [Hacktivist]: Extend file to required size
        if (ftruncate(fd_, file_size_) == -1) {
            close(fd_);
            throw std::runtime_error("Failed to resize cache file");
        }
        
        // [Hacktivist]: Memory-map the entire file for MAXIMUM SPEED
        mapped_memory_ = mmap(nullptr, file_size_, PROT_READ | PROT_WRITE, 
                             MAP_SHARED, fd_, 0);
        
        if (mapped_memory_ == MAP_FAILED) {
            close(fd_);
            throw std::runtime_error("Failed to memory-map cache file");
        }
        
        // [Performance Demon]: Advise kernel about access patterns
        madvise(mapped_memory_, file_size_, MADV_RANDOM);  // Random access pattern
    }
    
    ~mmap_file_manager() {
        if (mapped_memory_ && mapped_memory_ != MAP_FAILED) {
            // [Security Paranoid]: Secure sync before unmapping
            msync(mapped_memory_, file_size_, MS_SYNC);
            munmap(mapped_memory_, file_size_);
        }
        
        if (fd_ != -1) {
            close(fd_);
        }
    }
    
    // [Security Paranoid]: Bounds-checked entry access
    [[nodiscard]] entry_type* get_entry(size_t index) noexcept {
        if (validate_index(index)) {
            auto* entries = static_cast<entry_type*>(mapped_memory_);
            return &entries[index];
        }
        
        metrics_.validation_failures.fetch_add(1);
        return nullptr;
    }
    
    [[nodiscard]] const entry_type* get_entry(size_t index) const noexcept {
        return const_cast<mmap_file_manager*>(this)->get_entry(index);
    }
    
    [[nodiscard]] size_t max_entries() const noexcept {
        return file_size_ / sizeof(entry_type);
    }
    
    // [Hacktivist]: Force kernel to sync dirty pages
    mmap_error_code sync_to_disk() noexcept {
        if (msync(mapped_memory_, file_size_, MS_ASYNC) == -1) {
            return mmap_error_code::msync_failed;
        }
        
        metrics_.sync_operations.fetch_add(1);
        return mmap_error_code::success;
    }
    
    // [Performance Demon]: Prefault pages for immediate access
    void prefault_pages() noexcept {
        auto* bytes = static_cast<volatile char*>(mapped_memory_);
        const size_t page_size = getpagesize();
        
        for (size_t offset = 0; offset < file_size_; offset += page_size) {
            bytes[offset] = bytes[offset];  // Touch each page
        }
    }
    
    [[nodiscard]] const mmap_performance_metrics& get_metrics() const noexcept {
        return metrics_;
    }
    
private:
    [[nodiscard]] bool validate_index(size_t index) const noexcept {
        const size_t max_idx = max_entries();
        
        switch (security_level_) {
        case mmap_security_level::permissive:
            return index < max_idx;
            
        case mmap_security_level::strict:
            if (index >= max_idx) {
                return false;
            }
            // Additional alignment check
            return (index * sizeof(entry_type)) % alignof(entry_type) == 0;
            
        case mmap_security_level::paranoid:
            if (index >= max_idx) {
                return false;
            }
            // [Security Paranoid]: Full memory region validation
            const auto* entry_ptr = static_cast<const entry_type*>(mapped_memory_) + index;
            const auto* base = static_cast<const entry_type*>(mapped_memory_);
            const auto* end = base + max_idx;
            
            return entry_ptr >= base && entry_ptr < end;
        }
        
        return false;
    }
};

// [Enterprise Bean]: Memory-mapped cache with full enterprise features
template<typename T, size_t MaxEntries = 16384>
class mmap_cache_manager {
    static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable for mmap storage");
    static_assert(std::has_single_bit(MaxEntries), "MaxEntries must be power of 2");
    
public:
    using entry_type = mmap_cache_entry<T>;
    using file_manager_type = mmap_file_manager<T>;
    using key_type = std::string_view;
    
private:
    static constexpr size_t CAPACITY_MASK = MaxEntries - 1;
    static constexpr uint32_t DEFAULT_TTL = 3600;  // 1 hour for persistent cache
    static constexpr uint32_t MAX_PROBE_DISTANCE = 32;
    
    std::unique_ptr<file_manager_type> file_manager_;
    std::string cache_directory_;
    mmap_security_level security_level_;
    
    // [Performance Demon]: Atomic statistics for thread safety
    mutable std::atomic<uint64_t> stats_hits_{0};
    mutable std::atomic<uint64_t> stats_misses_{0};
    mutable std::atomic<uint64_t> stats_evictions_{0};
    
    [[nodiscard]] size_t hash_to_index(uint64_t hash) const noexcept {
        return hash & CAPACITY_MASK;
    }
    
    // [Hacktivist]: Simple FNV-1a hash because it's FAST and EFFECTIVE
    [[nodiscard]] uint64_t hash_key(key_type key) const noexcept {
        uint64_t hash = 14695981039346656037ULL;  // FNV offset basis
        for (char c : key) {
            hash ^= static_cast<uint64_t>(c);
            hash *= 1099511628211ULL;  // FNV prime
        }
        return hash;
    }
    
    // [Security Paranoid]: Secure slot finding with bounds checking
    [[nodiscard]] std::optional<size_t> find_slot(uint64_t key_hash, bool for_insertion = false) const noexcept {
        size_t index = hash_to_index(key_hash);
        size_t distance = 0;
        
        while (distance <= MAX_PROBE_DISTANCE) {
            const size_t current_index = (index + distance) & CAPACITY_MASK;
            const auto* entry = file_manager_->get_entry(current_index);
            
            if (!entry) {  // [Security Paranoid]: Bounds check failed
                return std::nullopt;
            }
            
            if (!entry->is_valid()) {
                // Invalid or empty slot
                return for_insertion ? std::optional{current_index} : std::nullopt;
            }
            
            if (entry->key_hash == key_hash && !entry->is_expired()) {
                return current_index;
            }
            
            ++distance;
        }
        
        return std::nullopt;
    }
    
    // [Minimalist Zen]: Simple LRU eviction
    [[nodiscard]] size_t find_eviction_candidate() const noexcept {
        uint64_t oldest_access = UINT64_MAX;
        size_t candidate = 0;
        
        for (size_t i = 0; i < MaxEntries; ++i) {
            const auto* entry = file_manager_->get_entry(i);
            if (!entry || !entry->is_valid()) {
                return i;  // Found invalid entry
            }
            
            if (entry->last_access_timestamp < oldest_access) {
                oldest_access = entry->last_access_timestamp;
                candidate = i;
            }
        }
        
        return candidate;
    }
    
public:
    // [Enterprise Bean]: Constructor with comprehensive configuration
    explicit mmap_cache_manager(const std::string& cache_dir = "/tmp/hsml_cache",
                               mmap_security_level security = mmap_security_level::strict)
        : cache_directory_(cache_dir), security_level_(security) {
        
        const std::string cache_file = cache_directory_ + "/spatial_cache.mmap";
        
        try {
            file_manager_ = std::make_unique<file_manager_type>(cache_file, MaxEntries, security);
            file_manager_->prefault_pages();  // [Performance Demon]: Warm up the cache
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to initialize mmap cache: " + std::string(e.what()));
        }
    }
    
    // [Hacktivist]: GET operation with zero-copy access
    template<typename Key>
    [[nodiscard]] std::optional<T> get(Key&& key) const noexcept {
        const uint64_t key_hash = hash_key(std::forward<Key>(key));
        const auto slot = find_slot(key_hash);
        
        if (!slot) {
            stats_misses_.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        auto* entry = file_manager_->get_entry(*slot);
        if (!entry || !entry->is_valid() || entry->is_expired()) {
            stats_misses_.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        // [Performance Demon]: Update access info directly in mapped memory
        const_cast<entry_type*>(entry)->update_access();
        stats_hits_.fetch_add(1, std::memory_order_relaxed);
        
        return entry->value;
    }
    
    // [Hacktivist]: SET operation with direct memory placement
    template<typename Key, typename Value>
    bool set(Key&& key, Value&& value, uint32_t ttl_seconds = DEFAULT_TTL) noexcept {
        const uint64_t key_hash = hash_key(std::forward<Key>(key));
        auto slot = find_slot(key_hash, true);
        
        if (!slot) {
            // Need to evict
            const size_t evict_index = find_eviction_candidate();
            slot = evict_index;
            stats_evictions_.fetch_add(1, std::memory_order_relaxed);
        }
        
        auto* entry = file_manager_->get_entry(*slot);
        if (!entry) {
            return false;  // [Security Paranoid]: Bounds check failed
        }
        
        // [Hacktivist]: Direct placement new into mapped memory
        new (entry) entry_type{key_hash, std::forward<Value>(value), ttl_seconds};
        
        return true;
    }
    
    // [Minimalist Zen]: Simple removal
    template<typename Key>
    bool remove(Key&& key) noexcept {
        const uint64_t key_hash = hash_key(std::forward<Key>(key));
        const auto slot = find_slot(key_hash);
        
        if (!slot) {
            return false;
        }
        
        auto* entry = file_manager_->get_entry(*slot);
        if (entry) {
            // [Security Paranoid]: Secure zeroing
            std::memset(entry, 0, sizeof(entry_type));
        }
        
        return true;
    }
    
    // [Performance Demon]: Batch cleanup of expired entries
    void cleanup_expired() noexcept {
        for (size_t i = 0; i < MaxEntries; ++i) {
            auto* entry = file_manager_->get_entry(i);
            if (entry && entry->is_expired()) {
                std::memset(entry, 0, sizeof(entry_type));
            }
        }
        
        // [Hacktivist]: Force sync to disk
        file_manager_->sync_to_disk();
    }
    
    // [Enterprise Bean]: Comprehensive health reporting
    struct health_report {
        uint64_t cache_hits;
        uint64_t cache_misses;
        uint64_t evictions;
        double hit_rate;
        size_t active_entries;
        double utilization;
        mmap_performance_metrics file_metrics;
        std::string cache_path;
    };
    
    [[nodiscard]] health_report health() const noexcept {
        const uint64_t hits = stats_hits_.load(std::memory_order_relaxed);
        const uint64_t misses = stats_misses_.load(std::memory_order_relaxed);
        const uint64_t total = hits + misses;
        
        size_t active = 0;
        for (size_t i = 0; i < MaxEntries; ++i) {
            const auto* entry = file_manager_->get_entry(i);
            if (entry && entry->is_valid() && !entry->is_expired()) {
                ++active;
            }
        }
        
        return health_report{
            .cache_hits = hits,
            .cache_misses = misses,
            .evictions = stats_evictions_.load(std::memory_order_relaxed),
            .hit_rate = total > 0 ? static_cast<double>(hits) / total : 0.0,
            .active_entries = active,
            .utilization = static_cast<double>(active) / MaxEntries,
            .file_metrics = file_manager_->get_metrics(),
            .cache_path = cache_directory_
        };
    }
    
    // [Security Paranoid]: Secure cache clearing
    void clear() noexcept {
        for (size_t i = 0; i < MaxEntries; ++i) {
            auto* entry = file_manager_->get_entry(i);
            if (entry) {
                secure_zero(entry, sizeof(entry_type));
            }
        }
        
        stats_hits_.store(0, std::memory_order_relaxed);
        stats_misses_.store(0, std::memory_order_relaxed);
        stats_evictions_.store(0, std::memory_order_relaxed);
        
        file_manager_->sync_to_disk();
    }
    
private:
    // [Security Paranoid]: Secure memory zeroing
    void secure_zero(void* ptr, size_t size) noexcept {
        volatile char* vptr = static_cast<volatile char*>(ptr);
        for (size_t i = 0; i < size; ++i) {
            vptr[i] = 0;
        }
    }
};

// [Enterprise Bean]: Factory for different mmap configurations
template<typename T>
class mmap_cache_factory {
public:
    enum class persistence_level {
        temporary,     // /tmp directory
        user_local,    // ~/.cache/hsml
        system_wide    // /var/cache/hsml
    };
    
    template<persistence_level Level, size_t Capacity = 16384>
    static auto create(mmap_security_level security = mmap_security_level::strict) {
        std::string cache_dir;
        
        if constexpr (Level == persistence_level::temporary) {
            cache_dir = "/tmp/hsml_cache";
        } else if constexpr (Level == persistence_level::user_local) {
            cache_dir = std::string(getenv("HOME")) + "/.cache/hsml";
        } else {
            cache_dir = "/var/cache/hsml";
        }
        
        return mmap_cache_manager<T, Capacity>{cache_dir, security};
    }
};

} // namespace compute
} // namespace hsml

/* [Hacktivist]: Memory-mapped files are the ULTIMATE performance hack! */
/* [Security Paranoid]: At least we have checksums and bounds checking... */
/* [Enterprise Bean]: The abstraction layers provide proper enterprise-grade functionality */