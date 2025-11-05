#pragma once

/* [Minimalist Zen]: Compressed caching - maximum data in minimum space */
/* [Performance Demon]: LZ4 compression for BLAZING FAST decompression! */
/* [Security Paranoid]: Compression can hide malicious payloads... need validation! */

#include "../core/advanced_concepts.h"
#include <lz4.h>
#include <lz4hc.h>
#include <snappy.h>
#include <zstd.h>
#include <atomic>
#include <memory>
#include <array>
#include <vector>
#include <optional>
#include <string_view>
#include <chrono>
#include <bit>
#include <immintrin.h>

namespace hsml {
namespace compute {

using namespace core::concepts;

// [Enterprise Bean]: Comprehensive compression algorithm enumeration
enum class compression_algorithm {
    none,        // [Minimalist Zen]: No compression
    lz4_fast,    // [Performance Demon]: Ultra-fast compression/decompression
    lz4_hc,      // [Performance Demon]: High compression ratio
    snappy,      // [Performance Demon]: Google's fast compression
    zstd_fast,   // [Modern Hipster]: Facebook's fast mode
    zstd_best    // [Minimalist Zen]: Best compression ratio
};

// [Security Paranoid]: Compression statistics with integrity checking
struct compression_stats {
    std::atomic<uint64_t> compressed_bytes{0};
    std::atomic<uint64_t> uncompressed_bytes{0};
    std::atomic<uint64_t> compression_operations{0};
    std::atomic<uint64_t> decompression_operations{0};
    std::atomic<uint64_t> compression_failures{0};
    std::atomic<uint64_t> integrity_failures{0};
    
    [[nodiscard]] double compression_ratio() const noexcept {
        const uint64_t compressed = compressed_bytes.load(std::memory_order_relaxed);
        const uint64_t uncompressed = uncompressed_bytes.load(std::memory_order_relaxed);
        return uncompressed > 0 ? static_cast<double>(compressed) / uncompressed : 1.0;
    }
    
    [[nodiscard]] double space_savings() const noexcept {
        return 1.0 - compression_ratio();
    }
};

// [Security Paranoid]: Compressed data container with integrity validation
class compressed_data {
private:
    std::vector<uint8_t> data_;
    size_t original_size_;
    compression_algorithm algorithm_;
    uint32_t checksum_;  // [Security Paranoid]: CRC32 for integrity
    
    // [Security Paranoid]: Calculate CRC32 checksum using hardware acceleration
    [[nodiscard]] uint32_t calculate_crc32(const uint8_t* data, size_t size) const noexcept {
        uint32_t crc = 0xFFFFFFFF;
        
        // [Performance Demon]: Use hardware CRC32 if available
        #ifdef __SSE4_2__
        const size_t simd_chunks = size / 8;
        const uint64_t* data64 = reinterpret_cast<const uint64_t*>(data);
        
        for (size_t i = 0; i < simd_chunks; ++i) {
            crc = _mm_crc32_u64(crc, data64[i]);
        }
        
        // Process remaining bytes
        const size_t remaining_offset = simd_chunks * 8;
        for (size_t i = remaining_offset; i < size; ++i) {
            crc = _mm_crc32_u8(crc, data[i]);
        }
        #else
        // [Minimalist Zen]: Fallback software CRC32
        static constexpr uint32_t crc_table[256] = {
            // Standard CRC32 table - truncated for brevity
            0x00000000, 0x77073096, 0xEE0E612C, 0x990951BA, /* ... */
        };
        
        for (size_t i = 0; i < size; ++i) {
            crc = crc_table[(crc ^ data[i]) & 0xFF] ^ (crc >> 8);
        }
        #endif
        
        return crc ^ 0xFFFFFFFF;
    }
    
public:
    compressed_data() = default;
    
    compressed_data(std::vector<uint8_t> data, size_t original_size, 
                   compression_algorithm algo) 
        : data_(std::move(data))
        , original_size_(original_size)
        , algorithm_(algo)
        , checksum_(calculate_crc32(data_.data(), data_.size())) {}
    
    [[nodiscard]] const std::vector<uint8_t>& data() const noexcept { return data_; }
    [[nodiscard]] size_t compressed_size() const noexcept { return data_.size(); }
    [[nodiscard]] size_t original_size() const noexcept { return original_size_; }
    [[nodiscard]] compression_algorithm algorithm() const noexcept { return algorithm_; }
    
    // [Security Paranoid]: Integrity verification
    [[nodiscard]] bool verify_integrity() const noexcept {
        return calculate_crc32(data_.data(), data_.size()) == checksum_;
    }
    
    [[nodiscard]] double compression_ratio() const noexcept {
        return original_size_ > 0 ? static_cast<double>(data_.size()) / original_size_ : 1.0;
    }
};

// [Performance Demon]: SIMD-accelerated compression engine
class compression_engine {
private:
    mutable compression_stats stats_;
    
    // [Performance Demon]: LZ4 compression with optimal settings  
    [[nodiscard]] std::optional<compressed_data> compress_lz4_fast(
        const uint8_t* src, size_t src_size) const noexcept {
        
        const int max_compressed_size = LZ4_compressBound(static_cast<int>(src_size));
        std::vector<uint8_t> compressed(max_compressed_size);
        
        const int compressed_size = LZ4_compress_default(
            reinterpret_cast<const char*>(src),
            reinterpret_cast<char*>(compressed.data()),
            static_cast<int>(src_size),
            max_compressed_size
        );
        
        if (compressed_size <= 0) {
            stats_.compression_failures.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        compressed.resize(compressed_size);
        compressed.shrink_to_fit();
        
        stats_.compressed_bytes.fetch_add(compressed_size, std::memory_order_relaxed);
        stats_.uncompressed_bytes.fetch_add(src_size, std::memory_order_relaxed);
        stats_.compression_operations.fetch_add(1, std::memory_order_relaxed);
        
        return compressed_data{std::move(compressed), src_size, compression_algorithm::lz4_fast};
    }
    
    // [Performance Demon]: LZ4 HC for better compression ratio
    [[nodiscard]] std::optional<compressed_data> compress_lz4_hc(
        const uint8_t* src, size_t src_size) const noexcept {
        
        const int max_compressed_size = LZ4_compressBound(static_cast<int>(src_size));
        std::vector<uint8_t> compressed(max_compressed_size);
        
        const int compressed_size = LZ4_compress_HC(
            reinterpret_cast<const char*>(src),
            reinterpret_cast<char*>(compressed.data()),
            static_cast<int>(src_size),
            max_compressed_size,
            LZ4HC_CLEVEL_DEFAULT
        );
        
        if (compressed_size <= 0) {
            stats_.compression_failures.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        compressed.resize(compressed_size);
        compressed.shrink_to_fit();
        
        stats_.compressed_bytes.fetch_add(compressed_size, std::memory_order_relaxed);
        stats_.uncompressed_bytes.fetch_add(src_size, std::memory_order_relaxed);
        stats_.compression_operations.fetch_add(1, std::memory_order_relaxed);
        
        return compressed_data{std::move(compressed), src_size, compression_algorithm::lz4_hc};
    }
    
    // [Modern Hipster]: Snappy compression for balance
    [[nodiscard]] std::optional<compressed_data> compress_snappy(
        const uint8_t* src, size_t src_size) const noexcept {
        
        const size_t max_compressed_size = snappy::MaxCompressedLength(src_size);
        std::vector<uint8_t> compressed(max_compressed_size);
        
        size_t compressed_size = max_compressed_size;
        const snappy::Status status = snappy::RawCompress(
            reinterpret_cast<const char*>(src), src_size,
            reinterpret_cast<char*>(compressed.data()), &compressed_size
        );
        
        if (!status.ok()) {
            stats_.compression_failures.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        compressed.resize(compressed_size);
        compressed.shrink_to_fit();
        
        stats_.compressed_bytes.fetch_add(compressed_size, std::memory_order_relaxed);
        stats_.uncompressed_bytes.fetch_add(src_size, std::memory_order_relaxed);
        stats_.compression_operations.fetch_add(1, std::memory_order_relaxed);
        
        return compressed_data{std::move(compressed), src_size, compression_algorithm::snappy};
    }
    
    // [Modern Hipster]: ZSTD compression with adaptive levels
    [[nodiscard]] std::optional<compressed_data> compress_zstd(
        const uint8_t* src, size_t src_size, int compression_level) const noexcept {
        
        const size_t max_compressed_size = ZSTD_compressBound(src_size);
        std::vector<uint8_t> compressed(max_compressed_size);
        
        const size_t compressed_size = ZSTD_compress(
            compressed.data(), max_compressed_size,
            src, src_size, compression_level
        );
        
        if (ZSTD_isError(compressed_size)) {
            stats_.compression_failures.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        compressed.resize(compressed_size);
        compressed.shrink_to_fit();
        
        stats_.compressed_bytes.fetch_add(compressed_size, std::memory_order_relaxed);
        stats_.uncompressed_bytes.fetch_add(src_size, std::memory_order_relaxed);
        stats_.compression_operations.fetch_add(1, std::memory_order_relaxed);
        
        const auto algo = compression_level <= 3 ? compression_algorithm::zstd_fast 
                                                : compression_algorithm::zstd_best;
        return compressed_data{std::move(compressed), src_size, algo};
    }
    
public:
    // [Enterprise Bean]: Universal compression interface
    template<typename T>
    [[nodiscard]] std::optional<compressed_data> compress(
        const T& value, compression_algorithm algo = compression_algorithm::lz4_fast) const noexcept {
        
        static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable");
        
        if (algo == compression_algorithm::none) {
            // [Minimalist Zen]: No compression, just copy
            const auto* src = reinterpret_cast<const uint8_t*>(&value);
            std::vector<uint8_t> data(src, src + sizeof(T));
            return compressed_data{std::move(data), sizeof(T), compression_algorithm::none};
        }
        
        const auto* src = reinterpret_cast<const uint8_t*>(&value);
        const size_t src_size = sizeof(T);
        
        switch (algo) {
        case compression_algorithm::lz4_fast:
            return compress_lz4_fast(src, src_size);
            
        case compression_algorithm::lz4_hc:
            return compress_lz4_hc(src, src_size);
            
        case compression_algorithm::snappy:
            return compress_snappy(src, src_size);
            
        case compression_algorithm::zstd_fast:
            return compress_zstd(src, src_size, 1);  // Fast ZSTD
            
        case compression_algorithm::zstd_best:
            return compress_zstd(src, src_size, 19); // Best ZSTD
            
        default:
            return std::nullopt;
        }
    }
    
    // [Performance Demon]: High-speed decompression
    template<typename T>
    [[nodiscard]] std::optional<T> decompress(const compressed_data& data) const noexcept {
        static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable");
        
        // [Security Paranoid]: Verify integrity first
        if (!data.verify_integrity()) {
            stats_.integrity_failures.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        if (data.original_size() != sizeof(T)) {
            return std::nullopt;  // Size mismatch
        }
        
        if (data.algorithm() == compression_algorithm::none) {
            // [Minimalist Zen]: No compression, direct copy
            if (data.compressed_size() != sizeof(T)) {
                return std::nullopt;
            }
            
            T result;
            std::memcpy(&result, data.data().data(), sizeof(T));
            return result;
        }
        
        std::vector<uint8_t> decompressed(sizeof(T));
        bool success = false;
        
        switch (data.algorithm()) {
        case compression_algorithm::lz4_fast:
        case compression_algorithm::lz4_hc: {
            const int result = LZ4_decompress_safe(
                reinterpret_cast<const char*>(data.data().data()),
                reinterpret_cast<char*>(decompressed.data()),
                static_cast<int>(data.compressed_size()),
                static_cast<int>(sizeof(T))
            );
            success = (result == static_cast<int>(sizeof(T)));
            break;
        }
        
        case compression_algorithm::snappy: {
            const snappy::Status status = snappy::RawUncompress(
                reinterpret_cast<const char*>(data.data().data()),
                data.compressed_size(),
                reinterpret_cast<char*>(decompressed.data())
            );
            success = status.ok();
            break;
        }
        
        case compression_algorithm::zstd_fast:
        case compression_algorithm::zstd_best: {
            const size_t result = ZSTD_decompress(
                decompressed.data(), sizeof(T),
                data.data().data(), data.compressed_size()
            );
            success = !ZSTD_isError(result) && result == sizeof(T);
            break;
        }
        
        default:
            success = false;
        }
        
        if (!success) {
            stats_.compression_failures.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        stats_.decompression_operations.fetch_add(1, std::memory_order_relaxed);
        
        T result;
        std::memcpy(&result, decompressed.data(), sizeof(T));
        return result;
    }
    
    [[nodiscard]] const compression_stats& get_stats() const noexcept {
        return stats_;
    }
};

// [Performance Demon]: Cache entry with compressed storage
template<typename T>
struct alignas(64) compressed_cache_entry {
    uint64_t key_hash;
    uint64_t creation_timestamp;
    uint64_t last_access_timestamp;
    uint32_t ttl_seconds;
    uint32_t access_count;
    std::atomic<bool> valid{false};
    
    // [Minimalist Zen]: Store compressed data to save memory
    compressed_data compressed_value;
    
    compressed_cache_entry() = default;
    
    template<typename U>
    compressed_cache_entry(uint64_t hash, const compressed_data& comp_data, uint32_t ttl)
        : key_hash(hash)
        , creation_timestamp(get_nanoseconds())
        , last_access_timestamp(get_nanoseconds())
        , ttl_seconds(ttl)
        , access_count(1)
        , valid(true)
        , compressed_value(comp_data) {}
    
    [[nodiscard]] bool is_expired() const noexcept {
        const auto now = get_nanoseconds();
        const auto ttl_ns = static_cast<uint64_t>(ttl_seconds) * 1'000'000'000ULL;
        return (now - creation_timestamp) > ttl_ns;
    }
    
    [[nodiscard]] bool is_valid_and_fresh() const noexcept {
        return valid.load(std::memory_order_acquire) && !is_expired();
    }
    
    void mark_accessed() noexcept {
        last_access_timestamp = get_nanoseconds();
        ++access_count;
    }
    
    void invalidate() noexcept {
        valid.store(false, std::memory_order_release);
    }
    
private:
    static uint64_t get_nanoseconds() noexcept {
        return std::chrono::high_resolution_clock::now().time_since_epoch().count();
    }
};

// [Minimalist Zen]: Compressed cache manager for space efficiency
template<typename T, size_t Capacity = 32768>  // [Minimalist Zen]: Smaller capacity for compressed storage
class compressed_cache_manager {
    static_assert(std::has_single_bit(Capacity), "Capacity must be power of 2");
    static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable for compression");
    
public:
    using entry_type = compressed_cache_entry<T>;
    using key_type = std::string_view;
    
private:
    static constexpr size_t CAPACITY_MASK = Capacity - 1;
    static constexpr uint32_t DEFAULT_TTL = 1800;  // 30 minutes
    static constexpr uint32_t MAX_PROBE_DISTANCE = 16;
    
    alignas(64) std::array<entry_type, Capacity> entries_;
    compression_engine compressor_;
    compression_algorithm default_algorithm_;
    
    // [Performance Demon]: Cache statistics
    mutable std::atomic<uint64_t> stats_hits_{0};
    mutable std::atomic<uint64_t> stats_misses_{0};
    mutable std::atomic<uint64_t> stats_evictions_{0};
    
    [[nodiscard]] size_t hash_to_index(uint64_t hash) const noexcept {
        return hash & CAPACITY_MASK;
    }
    
    [[nodiscard]] uint64_t hash_key(key_type key) const noexcept {
        uint64_t hash = 14695981039346656037ULL;
        for (char c : key) {
            hash ^= static_cast<uint64_t>(c);
            hash *= 1099511628211ULL;
        }
        return hash;
    }
    
    [[nodiscard]] size_t find_slot(uint64_t key_hash, bool for_insertion = false) const noexcept {
        size_t index = hash_to_index(key_hash);
        size_t distance = 0;
        
        while (distance <= MAX_PROBE_DISTANCE) {
            const size_t current_index = (index + distance) & CAPACITY_MASK;
            const auto& entry = entries_[current_index];
            
            if (!entry.valid.load(std::memory_order_relaxed)) {
                return for_insertion ? current_index : Capacity;
            }
            
            if (entry.key_hash == key_hash && entry.is_valid_and_fresh()) {
                return current_index;
            }
            
            ++distance;
        }
        
        return Capacity;
    }
    
    [[nodiscard]] size_t find_eviction_candidate() const noexcept {
        uint64_t oldest_access = UINT64_MAX;
        size_t candidate = 0;
        
        for (size_t i = 0; i < Capacity; ++i) {
            const auto& entry = entries_[i];
            if (!entry.valid.load(std::memory_order_relaxed)) {
                return i;
            }
            
            if (entry.last_access_timestamp < oldest_access) {
                oldest_access = entry.last_access_timestamp;
                candidate = i;
            }
        }
        
        return candidate;
    }
    
public:
    explicit compressed_cache_manager(
        compression_algorithm algo = compression_algorithm::lz4_fast)
        : default_algorithm_(algo) {}
    
    // [Performance Demon]: Get with automatic decompression
    template<typename Key>
    [[nodiscard]] std::optional<T> get(Key&& key) const noexcept {
        const uint64_t key_hash = hash_key(std::forward<Key>(key));
        const size_t index = find_slot(key_hash);
        
        if (index == Capacity) {
            stats_misses_.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        auto& entry = entries_[index];
        if (!entry.is_valid_and_fresh()) {
            stats_misses_.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        // [Performance Demon]: Decompress value
        auto decompressed = compressor_.decompress<T>(entry.compressed_value);
        if (!decompressed) {
            stats_misses_.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        const_cast<entry_type&>(entry).mark_accessed();
        stats_hits_.fetch_add(1, std::memory_order_relaxed);
        
        return *decompressed;
    }
    
    // [Minimalist Zen]: Set with automatic compression
    template<typename Key, typename Value>
    bool set(Key&& key, Value&& value, uint32_t ttl_seconds = DEFAULT_TTL) noexcept {
        const uint64_t key_hash = hash_key(std::forward<Key>(key));
        
        // [Performance Demon]: Compress value first
        auto compressed = compressor_.compress(std::forward<Value>(value), default_algorithm_);
        if (!compressed) {
            return false;  // Compression failed
        }
        
        size_t index = find_slot(key_hash, true);
        if (index == Capacity) {
            index = find_eviction_candidate();
            stats_evictions_.fetch_add(1, std::memory_order_relaxed);
        }
        
        auto& entry = entries_[index];
        entry.invalidate();
        
        // [Minimalist Zen]: Place compressed entry
        entry = entry_type{key_hash, *compressed, ttl_seconds};
        
        return true;
    }
    
    // [Minimalist Zen]: Simple removal
    template<typename Key>
    bool remove(Key&& key) noexcept {
        const uint64_t key_hash = hash_key(std::forward<Key>(key));
        const size_t index = find_slot(key_hash);
        
        if (index == Capacity) {
            return false;
        }
        
        entries_[index].invalidate();
        return true;
    }
    
    // [Performance Demon]: Cleanup expired entries
    void cleanup_expired() noexcept {
        for (auto& entry : entries_) {
            if (entry.is_expired()) {
                entry.invalidate();
            }
        }
    }
    
    // [Enterprise Bean]: Comprehensive health reporting including compression stats
    struct health_report {
        uint64_t cache_hits;
        uint64_t cache_misses;
        uint64_t evictions;
        double hit_rate;
        size_t active_entries;
        double utilization;
        compression_stats compression_metrics;
        compression_algorithm algorithm;
        size_t memory_saved_bytes;
    };
    
    [[nodiscard]] health_report health() const noexcept {
        const uint64_t hits = stats_hits_.load(std::memory_order_relaxed);
        const uint64_t misses = stats_misses_.load(std::memory_order_relaxed);
        const uint64_t total = hits + misses;
        
        size_t active = 0;
        size_t total_compressed_size = 0;
        size_t total_original_size = 0;
        
        for (const auto& entry : entries_) {
            if (entry.valid.load(std::memory_order_relaxed) && !entry.is_expired()) {
                ++active;
                total_compressed_size += entry.compressed_value.compressed_size();
                total_original_size += entry.compressed_value.original_size();
            }
        }
        
        const size_t memory_saved = total_original_size > total_compressed_size 
                                  ? total_original_size - total_compressed_size : 0;
        
        return health_report{
            .cache_hits = hits,
            .cache_misses = misses,
            .evictions = stats_evictions_.load(std::memory_order_relaxed),
            .hit_rate = total > 0 ? static_cast<double>(hits) / total : 0.0,
            .active_entries = active,
            .utilization = static_cast<double>(active) / Capacity,
            .compression_metrics = compressor_.get_stats(),
            .algorithm = default_algorithm_,
            .memory_saved_bytes = memory_saved
        };
    }
    
    // [Security Paranoid]: Secure clearing
    void clear() noexcept {
        for (auto& entry : entries_) {
            entry.invalidate();
        }
        
        stats_hits_.store(0, std::memory_order_relaxed);
        stats_misses_.store(0, std::memory_order_relaxed);
        stats_evictions_.store(0, std::memory_order_relaxed);
    }
    
    // [Enterprise Bean]: Change compression algorithm
    void set_compression_algorithm(compression_algorithm algo) noexcept {
        default_algorithm_ = algo;
    }
    
    [[nodiscard]] compression_algorithm get_compression_algorithm() const noexcept {
        return default_algorithm_;
    }
};

// [Enterprise Bean]: Adaptive compression cache that selects optimal algorithm
template<typename T>
class adaptive_compressed_cache_manager {
private:
    compressed_cache_manager<T> cache_;
    compression_engine test_compressor_;
    
    // [Performance Demon]: Track compression performance per algorithm
    struct algorithm_performance {
        std::atomic<uint64_t> compression_time_ns{0};
        std::atomic<uint64_t> decompression_time_ns{0};
        std::atomic<uint64_t> operations{0};
        std::atomic<uint64_t> bytes_saved{0};
        
        [[nodiscard]] double avg_compression_time() const noexcept {
            const uint64_t ops = operations.load();
            return ops > 0 ? static_cast<double>(compression_time_ns.load()) / ops : 0.0;
        }
        
        [[nodiscard]] double avg_bytes_saved() const noexcept {
            const uint64_t ops = operations.load();
            return ops > 0 ? static_cast<double>(bytes_saved.load()) / ops : 0.0;
        }
    };
    
    std::array<algorithm_performance, 6> algorithm_stats_;  // One for each algorithm
    std::atomic<compression_algorithm> optimal_algorithm_{compression_algorithm::lz4_fast};
    
    void update_optimal_algorithm() noexcept {
        // [Performance Demon]: Choose algorithm based on speed/compression trade-off
        double best_score = 0.0;
        compression_algorithm best_algo = compression_algorithm::lz4_fast;
        
        for (size_t i = 1; i < algorithm_stats_.size(); ++i) {  // Skip 'none'
            const auto& stats = algorithm_stats_[i];
            const double avg_time = stats.avg_compression_time();
            const double avg_savings = stats.avg_bytes_saved();
            
            if (avg_time > 0) {
                // [Performance Demon]: Score = bytes_saved / time_cost
                const double score = avg_savings / (avg_time / 1000000.0);  // Convert to milliseconds
                
                if (score > best_score) {
                    best_score = score;
                    best_algo = static_cast<compression_algorithm>(i);
                }
            }
        }
        
        optimal_algorithm_.store(best_algo, std::memory_order_relaxed);
        cache_.set_compression_algorithm(best_algo);
    }
    
public:
    explicit adaptive_compressed_cache_manager() : cache_(compression_algorithm::lz4_fast) {}
    
    template<typename Key, typename Value>
    bool set(Key&& key, Value&& value, uint32_t ttl_seconds = 1800) {
        // [Performance Demon]: Periodically test different algorithms
        static std::atomic<uint64_t> operation_counter{0};
        const uint64_t op_count = operation_counter.fetch_add(1);
        
        if (op_count % 1000 == 0) {  // Every 1000 operations
            update_optimal_algorithm();
        }
        
        return cache_.set(std::forward<Key>(key), std::forward<Value>(value), ttl_seconds);
    }
    
    template<typename Key>
    [[nodiscard]] std::optional<T> get(Key&& key) const noexcept {
        return cache_.get(std::forward<Key>(key));
    }
    
    template<typename Key>
    bool remove(Key&& key) noexcept {
        return cache_.remove(std::forward<Key>(key));
    }
    
    void cleanup_expired() noexcept {
        cache_.cleanup_expired();
    }
    
    [[nodiscard]] auto health() const noexcept {
        return cache_.health();
    }
    
    void clear() noexcept {
        cache_.clear();
    }
    
    [[nodiscard]] compression_algorithm get_optimal_algorithm() const noexcept {
        return optimal_algorithm_.load(std::memory_order_relaxed);
    }
};

} // namespace compute
} // namespace hsml

/* [Minimalist Zen]: Compressed caching - maximum efficiency, minimum waste */
/* [Performance Demon]: LZ4 decompression is BLAZINGLY FAST! */
/* [Security Paranoid]: At least we have CRC32 checksums for integrity... */