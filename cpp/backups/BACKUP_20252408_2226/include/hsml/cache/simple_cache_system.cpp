/**
 * HSML Simple Cache System
 * Essential caching without system-level complexity
 * 
 * Debloated from 736 lines to <150 lines while maintaining core functionality
 * 
 * @author Debloated by AspergicFullStackArchitect  
 * @version 8.0.0-DEBLOATED-EDITION
 */

#pragma once

#include <atomic>
#include <chrono>
#include <memory>
#include <shared_mutex>
#include <string>
#include <unordered_map>

namespace hsml::cache {

// Simple cache entry with TTL support
template<typename T>
struct CacheEntry {
    T value;
    std::chrono::steady_clock::time_point expiry;
    std::atomic<size_t> access_count{0};
    
    CacheEntry(T val, std::chrono::milliseconds ttl) 
        : value(std::move(val))
        , expiry(std::chrono::steady_clock::now() + ttl) {}
        
    bool is_expired() const {
        return std::chrono::steady_clock::now() > expiry;
    }
    
    void touch() {
        access_count.fetch_add(1, std::memory_order_relaxed);
    }
};

// Simple thread-safe cache
template<typename Key, typename Value>
class SimpleCache {
public:
    explicit SimpleCache(
        size_t max_size = 1000, 
        std::chrono::milliseconds default_ttl = std::chrono::minutes{10})
        : max_size_(max_size), default_ttl_(default_ttl) {}
    
    // Store value with TTL
    void put(const Key& key, Value value, 
             std::chrono::milliseconds ttl = std::chrono::milliseconds::zero()) {
        if (ttl == std::chrono::milliseconds::zero()) {
            ttl = default_ttl_;
        }
        
        std::unique_lock lock(mutex_);
        
        // Remove expired entries if cache is getting full
        if (cache_.size() >= max_size_) {
            cleanup_expired();
            
            // If still full, remove least recently used
            if (cache_.size() >= max_size_) {
                evict_lru();
            }
        }
        
        cache_[key] = std::make_unique<CacheEntry<Value>>(std::move(value), ttl);
        hit_count_.store(0, std::memory_order_relaxed);  // Reset stats on write
    }
    
    // Get value if exists and not expired
    std::optional<Value> get(const Key& key) {
        std::shared_lock lock(mutex_);
        
        auto it = cache_.find(key);
        if (it == cache_.end()) {
            miss_count_.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        if (it->second->is_expired()) {
            lock.unlock();
            std::unique_lock write_lock(mutex_);
            cache_.erase(it);
            miss_count_.fetch_add(1, std::memory_order_relaxed);
            return std::nullopt;
        }
        
        it->second->touch();
        hit_count_.fetch_add(1, std::memory_order_relaxed);
        return it->second->value;
    }
    
    // Remove specific key
    void remove(const Key& key) {
        std::unique_lock lock(mutex_);
        cache_.erase(key);
    }
    
    // Clear all entries
    void clear() {
        std::unique_lock lock(mutex_);
        cache_.clear();
        hit_count_.store(0, std::memory_order_relaxed);
        miss_count_.store(0, std::memory_order_relaxed);
    }
    
    // Get cache statistics
    struct Stats {
        size_t size;
        size_t hits;
        size_t misses;
        double hit_ratio;
    };
    
    Stats get_stats() const {
        std::shared_lock lock(mutex_);
        size_t hits = hit_count_.load(std::memory_order_relaxed);
        size_t misses = miss_count_.load(std::memory_order_relaxed);
        double ratio = (hits + misses > 0) ? static_cast<double>(hits) / (hits + misses) : 0.0;
        
        return {cache_.size(), hits, misses, ratio};
    }

private:
    mutable std::shared_mutex mutex_;
    std::unordered_map<Key, std::unique_ptr<CacheEntry<Value>>> cache_;
    const size_t max_size_;
    const std::chrono::milliseconds default_ttl_;
    
    mutable std::atomic<size_t> hit_count_{0};
    mutable std::atomic<size_t> miss_count_{0};
    
    // Remove expired entries (must hold unique lock)
    void cleanup_expired() {
        auto now = std::chrono::steady_clock::now();
        for (auto it = cache_.begin(); it != cache_.end();) {
            if (now > it->second->expiry) {
                it = cache_.erase(it);
            } else {
                ++it;
            }
        }
    }
    
    // Evict least recently used entry (must hold unique lock)
    void evict_lru() {
        if (cache_.empty()) return;
        
        auto lru_it = cache_.begin();
        size_t min_access = lru_it->second->access_count.load(std::memory_order_relaxed);
        
        for (auto it = cache_.begin(); it != cache_.end(); ++it) {
            size_t access = it->second->access_count.load(std::memory_order_relaxed);
            if (access < min_access) {
                min_access = access;
                lru_it = it;
            }
        }
        
        cache_.erase(lru_it);
    }
};

// Type aliases for common use cases
using StringCache = SimpleCache<std::string, std::string>;
using SpatialCache = SimpleCache<std::string, std::vector<double>>;

} // namespace hsml::cache