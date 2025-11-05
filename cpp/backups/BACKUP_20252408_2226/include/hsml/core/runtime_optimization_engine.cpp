/**
 * HSML Runtime Optimization Engine - C++20 Implementation
 * Revolutionary multi-threaded performance optimization for spherical coordinates
 * Death-cheating transposition with modern C++ paradigms
 */

#pragma once

#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <atomic>
#include <mutex>
#include <thread>
#include <chrono>
#include <concepts>
#include <ranges>
#include <algorithm>
#include <numeric>
#include <execution>
#include <optional>
#include <expected>
#include <span>
#include <coroutine>
#include <variant>

#include "spherical_coords.h"
#include "solid_angle.h"
#include "vector3.h"
#include "matrix4.h"

namespace hsml::core {

// Forward declarations
class SphericalCoordinateProcessor;
class SphericalPhysicsEngine;
class SphericalRenderer;

// Modern C++20 concepts for optimization system
template<typename T>
concept OptimizableObject = requires(T t) {
    { t.get_position() } -> std::convertible_to<SphericalCoords>;
    { t.get_distance_to(SphericalCoords{}) } -> std::convertible_to<double>;
    { t.is_visible() } -> std::convertible_to<bool>;
    { t.get_lod_level() } -> std::convertible_to<uint8_t>;
};

template<typename T>
concept PerformanceMonitor = requires(T t) {
    { t.get_frame_time() } -> std::convertible_to<double>;
    { t.get_memory_usage() } -> std::convertible_to<size_t>;
    { t.get_object_count() } -> std::convertible_to<size_t>;
};

// High-performance metrics structure with atomic operations
struct alignas(64) PerformanceMetrics {
    std::atomic<double> frame_time_ms{0.0};
    std::atomic<double> render_time_ms{0.0};
    std::atomic<double> physics_time_ms{0.0};
    std::atomic<size_t> memory_usage_bytes{0};
    std::atomic<size_t> object_count{0};
    std::atomic<size_t> solid_angle_calculations{0};
    std::atomic<double> cache_hit_rate{0.0};
    std::atomic<double> fps{0.0};
    
    enum class BottleneckType : uint8_t {
        NONE,
        RENDERING,
        PHYSICS,
        MEMORY,
        CALCULATION
    };
    
    std::atomic<BottleneckType> bottleneck_type{BottleneckType::NONE};
    
    // Copy constructor for atomic members
    PerformanceMetrics(const PerformanceMetrics& other) noexcept
        : frame_time_ms(other.frame_time_ms.load())
        , render_time_ms(other.render_time_ms.load())
        , physics_time_ms(other.physics_time_ms.load())
        , memory_usage_bytes(other.memory_usage_bytes.load())
        , object_count(other.object_count.load())
        , solid_angle_calculations(other.solid_angle_calculations.load())
        , cache_hit_rate(other.cache_hit_rate.load())
        , fps(other.fps.load())
        , bottleneck_type(other.bottleneck_type.load()) {}
    
    PerformanceMetrics& operator=(const PerformanceMetrics& other) noexcept {
        if (this != &other) {
            frame_time_ms.store(other.frame_time_ms.load());
            render_time_ms.store(other.render_time_ms.load());
            physics_time_ms.store(other.physics_time_ms.load());
            memory_usage_bytes.store(other.memory_usage_bytes.load());
            object_count.store(other.object_count.load());
            solid_angle_calculations.store(other.solid_angle_calculations.load());
            cache_hit_rate.store(other.cache_hit_rate.load());
            fps.store(other.fps.load());
            bottleneck_type.store(other.bottleneck_type.load());
        }
        return *this;
    }
    
    PerformanceMetrics() = default;
};

// Template-based optimization profile with compile-time validation
template<size_t MaxObjects = 10000>
struct OptimizationProfile {
    static constexpr size_t max_simultaneous_objects = MaxObjects;
    
    double lod_distance_threshold = 1000.0;
    double physics_time_step = 1.0 / 60.0;
    double culling_solid_angle_threshold = 0.001;
    std::chrono::milliseconds cache_expiry_time{5000};
    bool adaptive_quality_enabled = true;
    double performance_target_fps = 60.0;
    
    // Validation
    [[nodiscard]] constexpr bool is_valid() const noexcept {
        return lod_distance_threshold > 0.0 &&
               physics_time_step > 0.0 &&
               culling_solid_angle_threshold > 0.0 &&
               performance_target_fps > 0.0;
    }
};

// Level of Detail with constexpr calculations
struct LODLevel {
    double distance;
    uint32_t vertex_count;
    double material_complexity;
    double physics_accuracy;
    uint8_t render_priority;
    
    constexpr LODLevel(double d, uint32_t v, double m, double p, uint8_t r) noexcept
        : distance(d), vertex_count(v), material_complexity(m), 
          physics_accuracy(p), render_priority(r) {}
    
    [[nodiscard]] constexpr bool is_applicable(double object_distance) const noexcept {
        return object_distance >= distance;
    }
};

// Adaptive quality system with template metaprogramming
template<size_t NumLevels = 5>
struct AdaptiveQualitySystem {
    static constexpr size_t quality_levels = NumLevels;
    
    struct QualityLevel {
        std::string_view name;
        double multiplier;
        std::string_view description;
        
        constexpr QualityLevel(std::string_view n, double m, std::string_view d) noexcept
            : name(n), multiplier(m), description(d) {}
    };
    
    static constexpr std::array<QualityLevel, NumLevels> levels = {{
        {"potato", 0.25, "Minimum quality for low-end devices"},
        {"low", 0.5, "Reduced quality for better performance"},
        {"medium", 0.75, "Balanced quality and performance"},
        {"high", 1.0, "Full quality"},
        {"ultra", 1.25, "Enhanced quality for high-end devices"}
    }};
    
    std::atomic<size_t> current_level{3}; // Start at 'high'
    double target_frame_time_ms = 16.67; // 60 FPS
    double quality_adjustment_rate = 0.1;
    double performance_buffer = 2.0;
    std::atomic<std::chrono::steady_clock::time_point> last_adjustment;
    
    [[nodiscard]] constexpr const QualityLevel& get_current_level() const noexcept {
        return levels[current_level.load()];
    }
    
    [[nodiscard]] constexpr bool can_increase() const noexcept {
        return current_level.load() < quality_levels - 1;
    }
    
    [[nodiscard]] constexpr bool can_decrease() const noexcept {
        return current_level.load() > 0;
    }
};

// Lock-free spatial indexing system
template<size_t SectorCount = 1024>
class LockFreeSpatialIndex {
private:
    struct SpatialNode {
        std::string object_id;
        SphericalCoords position;
        std::atomic<SpatialNode*> next{nullptr};
        
        SpatialNode(std::string id, const SphericalCoords& pos)
            : object_id(std::move(id)), position(pos) {}
    };
    
    alignas(64) std::array<std::atomic<SpatialNode*>, SectorCount> sectors_;
    
public:
    void insert(const std::string& object_id, const SphericalCoords& position) {
        const size_t sector = calculate_sector_hash(position);
        auto* node = new SpatialNode(object_id, position);
        
        // Lock-free insertion
        auto* head = sectors_[sector].load();
        do {
            node->next.store(head);
        } while (!sectors_[sector].compare_exchange_weak(head, node));
    }
    
    [[nodiscard]] std::vector<std::string> query_radius(const SphericalCoords& center, double radius) const {
        std::vector<std::string> results;
        const size_t center_sector = calculate_sector_hash(center);
        
        // Search center sector and adjacent sectors
        for (size_t i = 0; i < SectorCount; ++i) {
            if (is_adjacent_sector(center_sector, i)) {
                auto* current = sectors_[i].load();
                while (current) {
                    const double distance = spherical_distance_fast(center, current->position);
                    if (distance <= radius) {
                        results.push_back(current->object_id);
                    }
                    current = current->next.load();
                }
            }
        }
        
        return results;
    }
    
    void clear() {
        for (auto& sector : sectors_) {
            auto* current = sector.load();
            while (current) {
                auto* next = current->next.load();
                delete current;
                current = next;
            }
            sector.store(nullptr);
        }
    }
    
    ~LockFreeSpatialIndex() { clear(); }

private:
    [[nodiscard]] constexpr size_t calculate_sector_hash(const SphericalCoords& pos) const noexcept {
        // Simple hash based on spherical coordinates
        const auto r_bucket = static_cast<size_t>(pos.r / 100.0);
        const auto theta_bucket = static_cast<size_t>((pos.theta / M_PI) * 16);
        const auto phi_bucket = static_cast<size_t>((pos.phi / (2.0 * M_PI)) * 32);
        
        return (r_bucket + theta_bucket * 17 + phi_bucket * 257) % SectorCount;
    }
    
    [[nodiscard]] constexpr bool is_adjacent_sector(size_t center, size_t candidate) const noexcept {
        // Simplified adjacency check
        const size_t diff = std::abs(static_cast<int64_t>(center) - static_cast<int64_t>(candidate));
        return diff <= 1 || diff == SectorCount - 1; // Handle wraparound
    }
    
    [[nodiscard]] static constexpr double spherical_distance_fast(const SphericalCoords& p1, 
                                                                 const SphericalCoords& p2) noexcept {
        // Fast approximation using haversine formula
        const double delta_theta = p2.theta - p1.theta;
        const double delta_phi = p2.phi - p1.phi;
        
        const double a = std::sin(delta_theta / 2.0) * std::sin(delta_theta / 2.0) +
                        std::cos(p1.theta) * std::cos(p2.theta) *
                        std::sin(delta_phi / 2.0) * std::sin(delta_phi / 2.0);
        
        return 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));
    }
};

// Memory pool for object reuse
template<typename T, size_t PoolSize = 1000>
class LockFreeMemoryPool {
private:
    alignas(64) std::array<T, PoolSize> pool_;
    std::atomic<size_t> next_free_{0};
    std::array<std::atomic<bool>, PoolSize> in_use_;
    
public:
    LockFreeMemoryPool() {
        for (auto& flag : in_use_) {
            flag.store(false);
        }
    }
    
    [[nodiscard]] std::optional<T*> acquire() {
        for (size_t attempts = 0; attempts < PoolSize; ++attempts) {
            const size_t index = next_free_.fetch_add(1) % PoolSize;
            
            bool expected = false;
            if (in_use_[index].compare_exchange_strong(expected, true)) {
                return &pool_[index];
            }
        }
        return std::nullopt; // Pool exhausted
    }
    
    void release(T* object) {
        const size_t index = object - pool_.data();
        if (index < PoolSize) {
            in_use_[index].store(false);
        }
    }
    
    [[nodiscard]] size_t available_count() const noexcept {
        return std::count_if(in_use_.begin(), in_use_.end(), 
                           [](const auto& flag) { return !flag.load(); });
    }
};

// Main runtime optimization engine
template<size_t MaxObjects = 10000>
class RuntimeOptimizationEngine {
private:
    // Singleton pattern with C++20 features
    static std::unique_ptr<RuntimeOptimizationEngine> instance_;
    static std::mutex instance_mutex_;
    
    // Core systems (dependency injection ready)
    std::shared_ptr<SphericalCoordinateProcessor> coordinate_processor_;
    std::shared_ptr<SphericalPhysicsEngine> physics_engine_;
    std::shared_ptr<SphericalRenderer> renderer_;
    
    // Performance monitoring
    PerformanceMetrics current_metrics_;
    std::vector<PerformanceMetrics> performance_history_;
    std::vector<double> frame_time_history_;
    mutable std::shared_mutex metrics_mutex_;
    
    // Optimization systems
    OptimizationProfile<MaxObjects> optimization_profile_;
    AdaptiveQualitySystem<5> adaptive_quality_;
    std::array<LODLevel, 5> lod_levels_;
    
    // Spatial optimization
    LockFreeSpatialIndex<1024> spatial_index_;
    std::unordered_map<std::string, bool> visibility_cache_;
    mutable std::shared_mutex visibility_cache_mutex_;
    
    // Memory management
    LockFreeMemoryPool<SphericalCoords, 10000> coord_pool_;
    size_t garbage_collection_threshold_ = 100 * 1024 * 1024; // 100MB
    std::atomic<std::chrono::steady_clock::time_point> last_gc_time_;
    
    // Thread management
    std::vector<std::thread> optimization_threads_;
    std::atomic<bool> should_terminate_{false};
    std::atomic<bool> optimization_enabled_{true};
    
    // Performance thresholds
    static constexpr struct {
        double critical_frame_time_ms = 33.33; // 30 FPS
        double target_frame_time_ms = 16.67;   // 60 FPS
        double excellent_frame_time_ms = 8.33; // 120 FPS
        size_t memory_warning_bytes = 150 * 1024 * 1024; // 150MB
        size_t cache_size_limit = 10000;
        std::chrono::milliseconds optimization_interval{1000};
    } PERFORMANCE_THRESHOLDS;
    
    // Private constructor for singleton
    RuntimeOptimizationEngine() {
        initialize_lod_levels();
        start_optimization_threads();
        last_gc_time_.store(std::chrono::steady_clock::now());
    }

public:
    // Singleton access with thread safety
    static auto get_instance() -> RuntimeOptimizationEngine& {
        std::lock_guard<std::mutex> lock(instance_mutex_);
        if (!instance_) {
            instance_ = std::unique_ptr<RuntimeOptimizationEngine>(
                new RuntimeOptimizationEngine()
            );
        }
        return *instance_;
    }
    
    // Dependency injection
    void set_coordinate_processor(std::shared_ptr<SphericalCoordinateProcessor> processor) {
        coordinate_processor_ = std::move(processor);
    }
    
    void set_physics_engine(std::shared_ptr<SphericalPhysicsEngine> engine) {
        physics_engine_ = std::move(engine);
    }
    
    void set_renderer(std::shared_ptr<SphericalRenderer> renderer) {
        renderer_ = std::move(renderer);
    }
    
    // Frame timing with RAII
    class FrameTimer {
    private:
        RuntimeOptimizationEngine& engine_;
        std::chrono::steady_clock::time_point start_time_;
        
    public:
        explicit FrameTimer(RuntimeOptimizationEngine& engine) 
            : engine_(engine), start_time_(std::chrono::steady_clock::now()) {}
        
        ~FrameTimer() {
            const auto end_time = std::chrono::steady_clock::now();
            const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
                end_time - start_time_).count() / 1000.0;
            
            engine_.end_frame(duration);
        }
    };
    
    [[nodiscard]] auto create_frame_timer() -> FrameTimer {
        return FrameTimer(*this);
    }
    
    // Performance monitoring
    void start_frame() {
        current_metrics_.frame_time_ms.store(
            std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::steady_clock::now().time_since_epoch()
            ).count() / 1000.0
        );
    }
    
    void end_frame(double frame_time_ms) {
        // Update frame time history
        frame_time_history_.push_back(frame_time_ms);
        if (frame_time_history_.size() > 60) {
            frame_time_history_.erase(frame_time_history_.begin());
        }
        
        // Calculate FPS
        const double fps = frame_time_ms > 0.0 ? 1000.0 / frame_time_ms : 0.0;
        current_metrics_.fps.store(fps);
        current_metrics_.frame_time_ms.store(frame_time_ms);
        
        // Update performance metrics
        update_performance_metrics();
        
        // Add to history
        {
            std::lock_guard<std::shared_mutex> lock(metrics_mutex_);
            performance_history_.push_back(current_metrics_);
            if (performance_history_.size() > 300) {
                performance_history_.erase(performance_history_.begin());
            }
        }
    }
    
    // LOD optimization
    template<OptimizableObject T>
    [[nodiscard]] auto calculate_lod_level(const T& object, const SphericalCoords& viewer_position) const -> const LODLevel& {
        const double distance = object.get_distance_to(viewer_position);
        
        // Find appropriate LOD level using ranges
        auto lod_iter = std::ranges::find_if(lod_levels_, 
            [distance](const LODLevel& lod) { return lod.is_applicable(distance); });
        
        return lod_iter != lod_levels_.end() ? *lod_iter : lod_levels_.back();
    }
    
    // Frustum culling with parallel execution
    template<std::ranges::range ObjectRange>
    [[nodiscard]] auto perform_frustum_culling(const ObjectRange& objects,
                                              const SphericalCoords& viewer_position,
                                              const SolidAngle& field_of_view) const 
        -> std::vector<std::string> {
        
        std::vector<std::string> visible_objects;
        
        // Parallel culling using C++20 execution policies
        std::mutex visible_mutex;
        std::for_each(std::execution::par_unseq, 
                     objects.begin(), objects.end(),
                     [&](const auto& object_pair) {
                         if (is_object_visible(object_pair.second, viewer_position, field_of_view)) {
                             std::lock_guard<std::mutex> lock(visible_mutex);
                             visible_objects.push_back(object_pair.first);
                         }
                     });
        
        return visible_objects;
    }
    
    // Spatial indexing
    void update_spatial_index(const std::unordered_map<std::string, SphericalCoords>& objects) {
        spatial_index_.clear();
        
        // Parallel insertion
        std::for_each(std::execution::par_unseq,
                     objects.begin(), objects.end(),
                     [this](const auto& obj_pair) {
                         spatial_index_.insert(obj_pair.first, obj_pair.second);
                     });
    }
    
    [[nodiscard]] auto get_objects_in_radius(const SphericalCoords& center, double radius) const 
        -> std::vector<std::string> {
        return spatial_index_.query_radius(center, radius);
    }
    
    // Memory management
    void perform_garbage_collection() {
        const auto now = std::chrono::steady_clock::now();
        const auto last_gc = last_gc_time_.load();
        
        if (now - last_gc < std::chrono::seconds(5)) {
            return; // Don't GC too frequently
        }
        
        // Clean up caches
        cleanup_caches();
        
        // Request coordinate processor cache cleanup
        if (coordinate_processor_) {
            // coordinate_processor_->cleanup_cache();
        }
        
        last_gc_time_.store(now);
    }
    
    // Coroutine-based async optimization
    struct OptimizationTask {
        struct promise_type {
            auto initial_suspend() { return std::suspend_never{}; }
            auto final_suspend() noexcept { return std::suspend_never{}; }
            auto get_return_object() { return OptimizationTask{}; }
            void return_void() {}
            void unhandled_exception() { std::terminate(); }
        };
    };
    
    auto optimize_async() -> OptimizationTask {
        run_optimizations();
        co_return;
    }
    
    // Public interface
    [[nodiscard]] auto get_performance_metrics() const -> PerformanceMetrics {
        return current_metrics_;
    }
    
    [[nodiscard]] auto get_optimization_profile() const -> const OptimizationProfile<MaxObjects>& {
        return optimization_profile_;
    }
    
    void update_optimization_profile(const OptimizationProfile<MaxObjects>& profile) {
        if (profile.is_valid()) {
            optimization_profile_ = profile;
        }
    }
    
    [[nodiscard]] auto get_current_quality_level() const -> std::tuple<size_t, std::string_view, std::string_view> {
        const auto& quality = adaptive_quality_.get_current_level();
        return {adaptive_quality_.current_level.load(), quality.name, quality.description};
    }
    
    void set_quality_level(size_t level) {
        if (level < AdaptiveQualitySystem<5>::quality_levels) {
            adaptive_quality_.current_level.store(level);
            apply_quality_level();
        }
    }
    
    void enable_optimization(bool enabled) {
        optimization_enabled_.store(enabled);
    }
    
    [[nodiscard]] std::string generate_performance_report() const {
        std::shared_lock<std::shared_mutex> lock(metrics_mutex_);
        
        if (performance_history_.empty()) {
            return "No performance data available";
        }
        
        // Calculate averages from recent history
        const size_t sample_size = std::min(size_t(10), performance_history_.size());
        const auto recent = std::span(performance_history_).last(sample_size);
        
        const double avg_fps = std::reduce(std::execution::par,
            recent.begin(), recent.end(), 0.0,
            [](double sum, const PerformanceMetrics& m) { 
                return sum + m.fps.load(); 
            }) / sample_size;
        
        const double avg_frame_time = std::reduce(std::execution::par,
            recent.begin(), recent.end(), 0.0,
            [](double sum, const PerformanceMetrics& m) { 
                return sum + m.frame_time_ms.load(); 
            }) / sample_size;
        
        const auto [level, name, description] = get_current_quality_level();
        
        return std::format(
            "Performance Report:\n"
            "  Average FPS: {:.1f}\n"
            "  Average Frame Time: {:.2f}ms\n"
            "  Quality Level: {}\n"
            "  Objects Rendered: {}\n"
            "  Cache Hit Rate: {:.1f}%\n"
            "  Memory Usage: {:.1f}MB\n"
            "  Current Bottleneck: {}\n",
            avg_fps,
            avg_frame_time,
            name,
            current_metrics_.object_count.load(),
            current_metrics_.cache_hit_rate.load() * 100.0,
            current_metrics_.memory_usage_bytes.load() / (1024.0 * 1024.0),
            bottleneck_to_string(current_metrics_.bottleneck_type.load())
        );
    }
    
    // Resource management
    ~RuntimeOptimizationEngine() {
        should_terminate_.store(true);
        
        for (auto& thread : optimization_threads_) {
            if (thread.joinable()) {
                thread.join();
            }
        }
    }
    
    // No copy/move semantics for singleton
    RuntimeOptimizationEngine(const RuntimeOptimizationEngine&) = delete;
    RuntimeOptimizationEngine& operator=(const RuntimeOptimizationEngine&) = delete;
    RuntimeOptimizationEngine(RuntimeOptimizationEngine&&) = delete;
    RuntimeOptimizationEngine& operator=(RuntimeOptimizationEngine&&) = delete;

private:
    void initialize_lod_levels() {
        lod_levels_ = {{
            {0.0,    1000, 1.0, 1.0, 1},
            {100.0,  500,  0.8, 0.9, 2},
            {500.0,  200,  0.6, 0.7, 3},
            {1000.0, 50,   0.4, 0.5, 4},
            {2000.0, 10,   0.2, 0.2, 5}
        }};
    }
    
    void start_optimization_threads() {
        const size_t thread_count = std::max(1u, std::thread::hardware_concurrency() / 4);
        
        for (size_t i = 0; i < thread_count; ++i) {
            optimization_threads_.emplace_back([this]() {
                optimization_worker_thread();
            });
        }
    }
    
    void optimization_worker_thread() {
        while (!should_terminate_.load()) {
            if (optimization_enabled_.load()) {
                run_optimizations();
            }
            
            std::this_thread::sleep_for(PERFORMANCE_THRESHOLDS.optimization_interval);
        }
    }
    
    void run_optimizations() {
        // Adjust quality based on performance
        adjust_quality_level();
        
        // Perform garbage collection if needed
        if (current_metrics_.memory_usage_bytes.load() > garbage_collection_threshold_) {
            perform_garbage_collection();
        }
        
        // Apply bottleneck-specific optimizations
        apply_bottleneck_optimizations();
    }
    
    void adjust_quality_level() {
        if (!optimization_profile_.adaptive_quality_enabled) {
            return;
        }
        
        const auto now = std::chrono::steady_clock::now();
        const auto last_adjustment = adaptive_quality_.last_adjustment.load();
        
        if (now - last_adjustment < std::chrono::milliseconds(500)) {
            return; // Debounce adjustments
        }
        
        if (frame_time_history_.empty()) {
            return;
        }
        
        const double avg_frame_time = std::reduce(frame_time_history_.begin(), 
                                                 frame_time_history_.end()) / frame_time_history_.size();
        const double target_frame_time = adaptive_quality_.target_frame_time_ms;
        const double buffer = adaptive_quality_.performance_buffer;
        
        if (avg_frame_time > target_frame_time * buffer && adaptive_quality_.can_decrease()) {
            // Performance is poor, reduce quality
            adaptive_quality_.current_level.fetch_sub(1);
            apply_quality_level();
            adaptive_quality_.last_adjustment.store(now);
        } else if (avg_frame_time < target_frame_time * 0.8 && adaptive_quality_.can_increase()) {
            // Performance is good, increase quality
            adaptive_quality_.current_level.fetch_add(1);
            apply_quality_level();
            adaptive_quality_.last_adjustment.store(now);
        }
    }
    
    void apply_quality_level() {
        const auto& quality_level = adaptive_quality_.get_current_level();
        const double multiplier = quality_level.multiplier;
        
        // Adjust LOD distances
        for (auto& lod : lod_levels_) {
            lod.distance *= multiplier;
            lod.vertex_count = static_cast<uint32_t>(lod.vertex_count * multiplier);
        }
        
        // Adjust optimization profile
        optimization_profile_.lod_distance_threshold *= multiplier;
        optimization_profile_.culling_solid_angle_threshold /= multiplier;
    }
    
    void apply_bottleneck_optimizations() {
        switch (current_metrics_.bottleneck_type.load()) {
            case PerformanceMetrics::BottleneckType::RENDERING:
                optimization_profile_.lod_distance_threshold *= 0.9;
                optimization_profile_.culling_solid_angle_threshold *= 1.1;
                break;
                
            case PerformanceMetrics::BottleneckType::PHYSICS:
                optimization_profile_.physics_time_step *= 1.1;
                break;
                
            case PerformanceMetrics::BottleneckType::MEMORY:
                perform_garbage_collection();
                optimization_profile_.cache_expiry_time = 
                    std::chrono::milliseconds(static_cast<int64_t>(
                        optimization_profile_.cache_expiry_time.count() * 0.8));
                break;
                
            case PerformanceMetrics::BottleneckType::CALCULATION:
                optimization_profile_.cache_expiry_time = 
                    std::chrono::milliseconds(static_cast<int64_t>(
                        optimization_profile_.cache_expiry_time.count() * 1.2));
                break;
                
            case PerformanceMetrics::BottleneckType::NONE:
                // No optimization needed
                break;
        }
    }
    
    void update_performance_metrics() {
        // Update memory usage
        // Note: This would need platform-specific implementation
        // current_metrics_.memory_usage_bytes.store(get_memory_usage());
        
        // Update other metrics from subsystems
        if (coordinate_processor_) {
            // auto stats = coordinate_processor_->get_performance_stats();
            // current_metrics_.cache_hit_rate.store(stats.cache_hit_rate);
            // current_metrics_.solid_angle_calculations.store(stats.total_computations);
        }
        
        if (renderer_) {
            // auto render_stats = renderer_->get_performance_stats();
            // current_metrics_.render_time_ms.store(render_stats.last_render_time);
            // current_metrics_.object_count.store(render_stats.total_objects);
        }
        
        if (physics_engine_) {
            // auto physics_stats = physics_engine_->get_simulation_stats();
            // current_metrics_.physics_time_ms.store(physics_stats.avg_frame_time);
        }
        
        // Identify bottleneck
        identify_bottleneck();
    }
    
    void identify_bottleneck() {
        const double frame_time = current_metrics_.frame_time_ms.load();
        const double render_time = current_metrics_.render_time_ms.load();
        const double physics_time = current_metrics_.physics_time_ms.load();
        const size_t memory_usage = current_metrics_.memory_usage_bytes.load();
        
        if (memory_usage > PERFORMANCE_THRESHOLDS.memory_warning_bytes) {
            current_metrics_.bottleneck_type.store(PerformanceMetrics::BottleneckType::MEMORY);
        } else if (render_time > frame_time * 0.6) {
            current_metrics_.bottleneck_type.store(PerformanceMetrics::BottleneckType::RENDERING);
        } else if (physics_time > frame_time * 0.3) {
            current_metrics_.bottleneck_type.store(PerformanceMetrics::BottleneckType::PHYSICS);
        } else if (current_metrics_.solid_angle_calculations.load() > 10000) {
            current_metrics_.bottleneck_type.store(PerformanceMetrics::BottleneckType::CALCULATION);
        } else {
            current_metrics_.bottleneck_type.store(PerformanceMetrics::BottleneckType::NONE);
        }
    }
    
    template<OptimizableObject T>
    [[nodiscard]] bool is_object_visible(const T& object,
                                        const SphericalCoords& viewer_position,
                                        const SolidAngle& field_of_view) const {
        
        const double distance = object.get_distance_to(viewer_position);
        const auto object_pos = object.get_position();
        
        // Calculate solid angle subtended by object
        const double angular_radius = std::atan(1.0 / distance); // Simplified radius
        const double object_solid_angle = 2.0 * M_PI * (1.0 - std::cos(angular_radius));
        
        // Cull objects that are too small
        if (object_solid_angle < optimization_profile_.culling_solid_angle_threshold) {
            return false;
        }
        
        // Check if object is within field of view
        const double theta_diff = std::abs(object_pos.theta - viewer_position.theta);
        const double phi_diff = std::abs(object_pos.phi - viewer_position.phi);
        
        return theta_diff < field_of_view.theta_extent() && 
               phi_diff < field_of_view.phi_extent();
    }
    
    void cleanup_caches() {
        std::lock_guard<std::shared_mutex> lock(visibility_cache_mutex_);
        visibility_cache_.clear();
    }
    
    [[nodiscard]] static constexpr std::string_view bottleneck_to_string(PerformanceMetrics::BottleneckType type) noexcept {
        switch (type) {
            case PerformanceMetrics::BottleneckType::NONE: return "none";
            case PerformanceMetrics::BottleneckType::RENDERING: return "rendering";
            case PerformanceMetrics::BottleneckType::PHYSICS: return "physics";
            case PerformanceMetrics::BottleneckType::MEMORY: return "memory";
            case PerformanceMetrics::BottleneckType::CALCULATION: return "calculation";
            default: return "unknown";
        }
    }
};

// Static member definitions
template<size_t MaxObjects>
std::unique_ptr<RuntimeOptimizationEngine<MaxObjects>> RuntimeOptimizationEngine<MaxObjects>::instance_;

template<size_t MaxObjects>
std::mutex RuntimeOptimizationEngine<MaxObjects>::instance_mutex_;

// Type aliases for common configurations
using StandardOptimizationEngine = RuntimeOptimizationEngine<10000>;
using HighPerformanceOptimizationEngine = RuntimeOptimizationEngine<50000>;
using LowEndOptimizationEngine = RuntimeOptimizationEngine<1000>;

} // namespace hsml::core