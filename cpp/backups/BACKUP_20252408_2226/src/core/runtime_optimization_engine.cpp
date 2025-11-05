/**
 * HSML Runtime Optimization Engine Implementation - C++20
 * Death-cheating transposition with advanced optimization algorithms
 */

#include "hsml/core/runtime_optimization_engine.h"
#include <cmath>
#include <random>
#include <filesystem>
#include <fstream>

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#elif defined(__linux__)
#include <sys/resource.h>
#include <fstream>
#elif defined(__APPLE__)
#include <mach/mach.h>
#include <sys/resource.h>
#endif

namespace hsml::core {

// Platform-specific memory usage detection
namespace platform_memory {
    
    [[nodiscard]] size_t get_current_memory_usage() noexcept {
#ifdef _WIN32
        PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
            return pmc.WorkingSetSize;
        }
#elif defined(__linux__)
        std::ifstream status_file("/proc/self/status");
        std::string line;
        while (std::getline(status_file, line)) {
            if (line.substr(0, 6) == "VmRSS:") {
                // Extract memory value in KB
                const size_t pos = line.find_first_of("0123456789");
                if (pos != std::string::npos) {
                    const size_t memory_kb = std::stoull(line.substr(pos));
                    return memory_kb * 1024; // Convert to bytes
                }
            }
        }
#elif defined(__APPLE__)
        task_basic_info info;
        mach_msg_type_number_t info_count = TASK_BASIC_INFO_COUNT;
        if (task_info(mach_task_self(), TASK_BASIC_INFO, 
                     reinterpret_cast<task_info_t>(&info), &info_count) == KERN_SUCCESS) {
            return info.resident_size;
        }
#endif
        return 0; // Fallback
    }
    
    [[nodiscard]] size_t get_available_memory() noexcept {
#ifdef _WIN32
        MEMORYSTATUSEX memory_status;
        memory_status.dwLength = sizeof(memory_status);
        if (GlobalMemoryStatusEx(&memory_status)) {
            return memory_status.ullAvailPhys;
        }
#elif defined(__linux__)
        std::ifstream meminfo("/proc/meminfo");
        std::string line;
        size_t available_memory = 0;
        
        while (std::getline(meminfo, line)) {
            if (line.substr(0, 12) == "MemAvailable") {
                const size_t pos = line.find_first_of("0123456789");
                if (pos != std::string::npos) {
                    available_memory = std::stoull(line.substr(pos)) * 1024;
                    break;
                }
            }
        }
        return available_memory;
#elif defined(__APPLE__)
        vm_size_t page_size;
        vm_statistics64_data_t vm_stat;
        mach_msg_type_number_t host_size = sizeof(vm_stat) / sizeof(natural_t);
        
        if (host_page_size(mach_host_self(), &page_size) == KERN_SUCCESS &&
            host_statistics64(mach_host_self(), HOST_VM_INFO64, 
                            reinterpret_cast<host_info64_t>(&vm_stat), &host_size) == KERN_SUCCESS) {
            return (vm_stat.free_count + vm_stat.inactive_count) * page_size;
        }
#endif
        return 4ULL * 1024 * 1024 * 1024; // 4GB fallback
    }
}

// SIMD-optimized distance calculations
namespace simd_math {
    
#ifdef __AVX2__
    #include <immintrin.h>
    
    // Vectorized spherical distance calculation for multiple points
    void calculate_distances_avx2(const SphericalCoords* positions,
                                 const SphericalCoords& viewer,
                                 double* distances,
                                 size_t count) {
        const __m256d viewer_theta = _mm256_set1_pd(viewer.theta);
        const __m256d viewer_phi = _mm256_set1_pd(viewer.phi);
        const __m256d viewer_r = _mm256_set1_pd(viewer.r);
        
        for (size_t i = 0; i < count; i += 4) {
            // Load 4 coordinates at once
            __m256d theta = _mm256_set_pd(
                (i + 3 < count) ? positions[i + 3].theta : 0.0,
                (i + 2 < count) ? positions[i + 2].theta : 0.0,
                (i + 1 < count) ? positions[i + 1].theta : 0.0,
                positions[i].theta
            );
            
            __m256d phi = _mm256_set_pd(
                (i + 3 < count) ? positions[i + 3].phi : 0.0,
                (i + 2 < count) ? positions[i + 2].phi : 0.0,
                (i + 1 < count) ? positions[i + 1].phi : 0.0,
                positions[i].phi
            );
            
            __m256d r = _mm256_set_pd(
                (i + 3 < count) ? positions[i + 3].r : 0.0,
                (i + 2 < count) ? positions[i + 2].r : 0.0,
                (i + 1 < count) ? positions[i + 1].r : 0.0,
                positions[i].r
            );
            
            // Calculate spherical law of cosines
            __m256d delta_theta = _mm256_sub_pd(theta, viewer_theta);
            __m256d delta_phi = _mm256_sub_pd(phi, viewer_phi);
            
            __m256d cos_delta_theta = _mm256_set_pd(
                std::cos(_mm256_extract_pd(delta_theta, 3)),
                std::cos(_mm256_extract_pd(delta_theta, 2)),
                std::cos(_mm256_extract_pd(delta_theta, 1)),
                std::cos(_mm256_extract_pd(delta_theta, 0))
            );
            
            __m256d cos_delta_phi = _mm256_set_pd(
                std::cos(_mm256_extract_pd(delta_phi, 3)),
                std::cos(_mm256_extract_pd(delta_phi, 2)),
                std::cos(_mm256_extract_pd(delta_phi, 1)),
                std::cos(_mm256_extract_pd(delta_phi, 0))
            );
            
            // Great circle distance calculation
            __m256d cos_distance = _mm256_add_pd(
                _mm256_mul_pd(
                    _mm256_set_pd(std::sin(viewer.theta), std::sin(viewer.theta), 
                                 std::sin(viewer.theta), std::sin(viewer.theta)),
                    _mm256_set_pd(
                        std::sin(_mm256_extract_pd(theta, 3)),
                        std::sin(_mm256_extract_pd(theta, 2)),
                        std::sin(_mm256_extract_pd(theta, 1)),
                        std::sin(_mm256_extract_pd(theta, 0))
                    )
                ),
                _mm256_mul_pd(
                    _mm256_mul_pd(
                        _mm256_set_pd(std::cos(viewer.theta), std::cos(viewer.theta),
                                     std::cos(viewer.theta), std::cos(viewer.theta)),
                        _mm256_set_pd(
                            std::cos(_mm256_extract_pd(theta, 3)),
                            std::cos(_mm256_extract_pd(theta, 2)),
                            std::cos(_mm256_extract_pd(theta, 1)),
                            std::cos(_mm256_extract_pd(theta, 0))
                        )
                    ),
                    cos_delta_phi
                )
            );
            
            // Store results
            alignas(32) double temp_distances[4];
            temp_distances[0] = std::acos(std::clamp(_mm256_extract_pd(cos_distance, 0), -1.0, 1.0));
            temp_distances[1] = std::acos(std::clamp(_mm256_extract_pd(cos_distance, 1), -1.0, 1.0));
            temp_distances[2] = std::acos(std::clamp(_mm256_extract_pd(cos_distance, 2), -1.0, 1.0));
            temp_distances[3] = std::acos(std::clamp(_mm256_extract_pd(cos_distance, 3), -1.0, 1.0));
            
            for (size_t j = 0; j < 4 && i + j < count; ++j) {
                distances[i + j] = temp_distances[j];
            }
        }
    }
#endif
    
    // Fallback scalar implementation
    void calculate_distances_scalar(const SphericalCoords* positions,
                                   const SphericalCoords& viewer,
                                   double* distances,
                                   size_t count) {
        const double viewer_sin_theta = std::sin(viewer.theta);
        const double viewer_cos_theta = std::cos(viewer.theta);
        
        for (size_t i = 0; i < count; ++i) {
            const double delta_phi = positions[i].phi - viewer.phi;
            const double cos_distance = 
                viewer_sin_theta * std::sin(positions[i].theta) * std::cos(delta_phi) +
                viewer_cos_theta * std::cos(positions[i].theta);
            
            distances[i] = std::acos(std::clamp(cos_distance, -1.0, 1.0));
        }
    }
}

// Advanced caching system with LRU eviction
template<typename Key, typename Value, size_t MaxSize = 10000>
class LRUCache {
private:
    struct CacheNode {
        Key key;
        Value value;
        std::chrono::steady_clock::time_point access_time;
        std::shared_ptr<CacheNode> prev;
        std::shared_ptr<CacheNode> next;
        
        CacheNode(Key k, Value v) 
            : key(std::move(k)), value(std::move(v)), 
              access_time(std::chrono::steady_clock::now()) {}
    };
    
    std::unordered_map<Key, std::shared_ptr<CacheNode>> cache_map_;
    std::shared_ptr<CacheNode> head_;
    std::shared_ptr<CacheNode> tail_;
    mutable std::shared_mutex cache_mutex_;
    std::atomic<size_t> hit_count_{0};
    std::atomic<size_t> miss_count_{0};
    
public:
    LRUCache() {
        head_ = std::make_shared<CacheNode>(Key{}, Value{});
        tail_ = std::make_shared<CacheNode>(Key{}, Value{});
        head_->next = tail_;
        tail_->prev = head_;
    }
    
    [[nodiscard]] std::optional<Value> get(const Key& key) {
        std::shared_lock<std::shared_mutex> lock(cache_mutex_);
        
        auto it = cache_map_.find(key);
        if (it != cache_map_.end()) {
            // Move to front (most recently used)
            move_to_front(it->second);
            hit_count_.fetch_add(1, std::memory_order_relaxed);
            return it->second->value;
        }
        
        miss_count_.fetch_add(1, std::memory_order_relaxed);
        return std::nullopt;
    }
    
    void put(const Key& key, Value value) {
        std::unique_lock<std::shared_mutex> lock(cache_mutex_);
        
        auto it = cache_map_.find(key);
        if (it != cache_map_.end()) {
            // Update existing value
            it->second->value = std::move(value);
            it->second->access_time = std::chrono::steady_clock::now();
            move_to_front(it->second);
            return;
        }
        
        // Add new entry
        auto new_node = std::make_shared<CacheNode>(key, std::move(value));
        cache_map_[key] = new_node;
        add_to_front(new_node);
        
        // Evict if necessary
        if (cache_map_.size() > MaxSize) {
            auto last_node = tail_->prev;
            remove_node(last_node);
            cache_map_.erase(last_node->key);
        }
    }
    
    void clear() {
        std::unique_lock<std::shared_mutex> lock(cache_mutex_);
        cache_map_.clear();
        head_->next = tail_;
        tail_->prev = head_;
        hit_count_.store(0);
        miss_count_.store(0);
    }
    
    [[nodiscard]] double get_hit_rate() const noexcept {
        const size_t hits = hit_count_.load(std::memory_order_relaxed);
        const size_t misses = miss_count_.load(std::memory_order_relaxed);
        const size_t total = hits + misses;
        return total > 0 ? static_cast<double>(hits) / total : 0.0;
    }
    
    [[nodiscard]] size_t size() const noexcept {
        std::shared_lock<std::shared_mutex> lock(cache_mutex_);
        return cache_map_.size();
    }

private:
    void move_to_front(std::shared_ptr<CacheNode> node) {
        remove_node(node);
        add_to_front(node);
    }
    
    void add_to_front(std::shared_ptr<CacheNode> node) {
        node->prev = head_;
        node->next = head_->next;
        head_->next->prev = node;
        head_->next = node;
    }
    
    void remove_node(std::shared_ptr<CacheNode> node) {
        node->prev->next = node->next;
        node->next->prev = node->prev;
    }
};

// Global optimization state management
struct GlobalOptimizationState {
    std::atomic<bool> global_optimization_enabled{true};
    std::atomic<double> global_performance_target{60.0}; // Target FPS
    std::atomic<size_t> total_objects_rendered{0};
    std::atomic<double> total_render_time{0.0};
    std::chrono::steady_clock::time_point startup_time;
    
    // Performance profiling
    LRUCache<std::string, double> performance_cache;
    
    GlobalOptimizationState() : startup_time(std::chrono::steady_clock::now()) {}
    
    [[nodiscard]] std::chrono::milliseconds get_uptime() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - startup_time
        );
    }
};

static GlobalOptimizationState g_optimization_state;

// Advanced optimization algorithms
namespace optimization_algorithms {
    
    // Genetic algorithm for optimization profile tuning
    struct OptimizationGene {
        double lod_threshold_factor = 1.0;
        double culling_threshold_factor = 1.0;
        double physics_step_factor = 1.0;
        double cache_expiry_factor = 1.0;
        double fitness = 0.0;
    };
    
    class GeneticOptimizer {
    private:
        static constexpr size_t POPULATION_SIZE = 20;
        static constexpr size_t GENERATIONS = 10;
        static constexpr double MUTATION_RATE = 0.1;
        static constexpr double CROSSOVER_RATE = 0.7;
        
        std::vector<OptimizationGene> population_;
        std::mt19937 rng_;
        std::uniform_real_distribution<double> dist_{0.5, 2.0};
        std::uniform_real_distribution<double> mutation_dist_{0.9, 1.1};
        
    public:
        GeneticOptimizer() : rng_(std::random_device{}()) {
            initialize_population();
        }
        
        [[nodiscard]] OptimizationGene evolve(const PerformanceMetrics& target_metrics) {
            for (size_t generation = 0; generation < GENERATIONS; ++generation) {
                // Evaluate fitness
                for (auto& gene : population_) {
                    gene.fitness = calculate_fitness(gene, target_metrics);
                }
                
                // Sort by fitness
                std::sort(population_.begin(), population_.end(),
                         [](const auto& a, const auto& b) { return a.fitness > b.fitness; });
                
                // Create next generation
                std::vector<OptimizationGene> next_generation;
                next_generation.reserve(POPULATION_SIZE);
                
                // Keep best individuals (elitism)
                for (size_t i = 0; i < POPULATION_SIZE / 4; ++i) {
                    next_generation.push_back(population_[i]);
                }
                
                // Crossover and mutation
                while (next_generation.size() < POPULATION_SIZE) {
                    const auto& parent1 = select_parent();
                    const auto& parent2 = select_parent();
                    
                    auto offspring = crossover(parent1, parent2);
                    mutate(offspring);
                    next_generation.push_back(offspring);
                }
                
                population_ = std::move(next_generation);
            }
            
            // Return best gene
            return *std::max_element(population_.begin(), population_.end(),
                                   [](const auto& a, const auto& b) { return a.fitness < b.fitness; });
        }
        
    private:
        void initialize_population() {
            population_.reserve(POPULATION_SIZE);
            for (size_t i = 0; i < POPULATION_SIZE; ++i) {
                OptimizationGene gene;
                gene.lod_threshold_factor = dist_(rng_);
                gene.culling_threshold_factor = dist_(rng_);
                gene.physics_step_factor = dist_(rng_);
                gene.cache_expiry_factor = dist_(rng_);
                population_.push_back(gene);
            }
        }
        
        [[nodiscard]] double calculate_fitness(const OptimizationGene& gene, 
                                              const PerformanceMetrics& target) const {
            // Fitness based on how well the gene would perform
            double fitness = 0.0;
            
            // Reward good FPS
            const double fps_score = std::min(1.0, target.fps.load() / 60.0);
            fitness += fps_score * 0.4;
            
            // Reward low memory usage
            const double memory_score = 1.0 - std::min(1.0, 
                target.memory_usage_bytes.load() / (200.0 * 1024 * 1024));
            fitness += memory_score * 0.3;
            
            // Reward high cache hit rate
            fitness += target.cache_hit_rate.load() * 0.3;
            
            return fitness;
        }
        
        [[nodiscard]] const OptimizationGene& select_parent() const {
            // Tournament selection
            const size_t tournament_size = 3;
            size_t best_idx = rng_() % population_.size();
            double best_fitness = population_[best_idx].fitness;
            
            for (size_t i = 1; i < tournament_size; ++i) {
                const size_t idx = rng_() % population_.size();
                if (population_[idx].fitness > best_fitness) {
                    best_idx = idx;
                    best_fitness = population_[idx].fitness;
                }
            }
            
            return population_[best_idx];
        }
        
        [[nodiscard]] OptimizationGene crossover(const OptimizationGene& parent1, 
                                                const OptimizationGene& parent2) const {
            if (std::uniform_real_distribution<double>{0.0, 1.0}(rng_) > CROSSOVER_RATE) {
                return parent1; // No crossover
            }
            
            OptimizationGene offspring;
            offspring.lod_threshold_factor = rng_() % 2 ? parent1.lod_threshold_factor : parent2.lod_threshold_factor;
            offspring.culling_threshold_factor = rng_() % 2 ? parent1.culling_threshold_factor : parent2.culling_threshold_factor;
            offspring.physics_step_factor = rng_() % 2 ? parent1.physics_step_factor : parent2.physics_step_factor;
            offspring.cache_expiry_factor = rng_() % 2 ? parent1.cache_expiry_factor : parent2.cache_expiry_factor;
            
            return offspring;
        }
        
        void mutate(OptimizationGene& gene) const {
            if (std::uniform_real_distribution<double>{0.0, 1.0}(rng_) < MUTATION_RATE) {
                gene.lod_threshold_factor *= mutation_dist_(rng_);
            }
            if (std::uniform_real_distribution<double>{0.0, 1.0}(rng_) < MUTATION_RATE) {
                gene.culling_threshold_factor *= mutation_dist_(rng_);
            }
            if (std::uniform_real_distribution<double>{0.0, 1.0}(rng_) < MUTATION_RATE) {
                gene.physics_step_factor *= mutation_dist_(rng_);
            }
            if (std::uniform_real_distribution<double>{0.0, 1.0}(rng_) < MUTATION_RATE) {
                gene.cache_expiry_factor *= mutation_dist_(rng_);
            }
        }
    };
    
    // Particle Swarm Optimization for real-time parameter tuning
    class ParticleSwarmOptimizer {
    private:
        struct Particle {
            std::vector<double> position;
            std::vector<double> velocity;
            std::vector<double> best_position;
            double best_fitness = -std::numeric_limits<double>::max();
        };
        
        static constexpr size_t SWARM_SIZE = 15;
        static constexpr size_t DIMENSIONS = 4; // lod, culling, physics, cache
        static constexpr double INERTIA = 0.729;
        static constexpr double COGNITIVE = 1.49445;
        static constexpr double SOCIAL = 1.49445;
        
        std::vector<Particle> swarm_;
        std::vector<double> global_best_position_;
        double global_best_fitness_ = -std::numeric_limits<double>::max();
        std::mt19937 rng_;
        
    public:
        ParticleSwarmOptimizer() : rng_(std::random_device{}()) {
            initialize_swarm();
        }
        
        [[nodiscard]] std::vector<double> optimize(std::function<double(const std::vector<double>&)> fitness_func,
                                                  size_t iterations = 50) {
            for (size_t iter = 0; iter < iterations; ++iter) {
                for (auto& particle : swarm_) {
                    // Evaluate fitness
                    const double fitness = fitness_func(particle.position);
                    
                    // Update personal best
                    if (fitness > particle.best_fitness) {
                        particle.best_fitness = fitness;
                        particle.best_position = particle.position;
                    }
                    
                    // Update global best
                    if (fitness > global_best_fitness_) {
                        global_best_fitness_ = fitness;
                        global_best_position_ = particle.position;
                    }
                }
                
                // Update velocities and positions
                for (auto& particle : swarm_) {
                    for (size_t d = 0; d < DIMENSIONS; ++d) {
                        const double r1 = std::uniform_real_distribution<double>{0.0, 1.0}(rng_);
                        const double r2 = std::uniform_real_distribution<double>{0.0, 1.0}(rng_);
                        
                        particle.velocity[d] = INERTIA * particle.velocity[d] +
                                             COGNITIVE * r1 * (particle.best_position[d] - particle.position[d]) +
                                             SOCIAL * r2 * (global_best_position_[d] - particle.position[d]);
                        
                        particle.position[d] += particle.velocity[d];
                        
                        // Clamp to bounds
                        particle.position[d] = std::clamp(particle.position[d], 0.1, 5.0);
                    }
                }
            }
            
            return global_best_position_;
        }
        
    private:
        void initialize_swarm() {
            swarm_.resize(SWARM_SIZE);
            global_best_position_.resize(DIMENSIONS);
            
            std::uniform_real_distribution<double> pos_dist(0.5, 2.0);
            std::uniform_real_distribution<double> vel_dist(-0.1, 0.1);
            
            for (auto& particle : swarm_) {
                particle.position.resize(DIMENSIONS);
                particle.velocity.resize(DIMENSIONS);
                particle.best_position.resize(DIMENSIONS);
                
                for (size_t d = 0; d < DIMENSIONS; ++d) {
                    particle.position[d] = pos_dist(rng_);
                    particle.velocity[d] = vel_dist(rng_);
                    particle.best_position[d] = particle.position[d];
                }
            }
            
            // Initialize global best
            for (size_t d = 0; d < DIMENSIONS; ++d) {
                global_best_position_[d] = 1.0; // Default values
            }
        }
    };
}

// Template instantiations for common configurations
template class RuntimeOptimizationEngine<1000>;   // LowEndOptimizationEngine
template class RuntimeOptimizationEngine<10000>;  // StandardOptimizationEngine
template class RuntimeOptimizationEngine<50000>;  // HighPerformanceOptimizationEngine

// Performance logging and analytics
class PerformanceLogger {
private:
    std::ofstream log_file_;
    std::mutex log_mutex_;
    bool logging_enabled_ = false;
    
public:
    void enable_logging(const std::string& filename) {
        std::lock_guard<std::mutex> lock(log_mutex_);
        log_file_.open(filename, std::ios::app);
        logging_enabled_ = log_file_.is_open();
        
        if (logging_enabled_) {
            log_file_ << "timestamp,fps,frame_time_ms,render_time_ms,physics_time_ms,"
                      << "memory_mb,object_count,cache_hit_rate,bottleneck_type\n";
        }
    }
    
    void log_metrics(const PerformanceMetrics& metrics) {
        if (!logging_enabled_) return;
        
        std::lock_guard<std::mutex> lock(log_mutex_);
        const auto now = std::chrono::system_clock::now();
        const auto timestamp = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()).count();
        
        log_file_ << timestamp << ","
                  << metrics.fps.load() << ","
                  << metrics.frame_time_ms.load() << ","
                  << metrics.render_time_ms.load() << ","
                  << metrics.physics_time_ms.load() << ","
                  << (metrics.memory_usage_bytes.load() / (1024.0 * 1024.0)) << ","
                  << metrics.object_count.load() << ","
                  << metrics.cache_hit_rate.load() << ","
                  << static_cast<int>(metrics.bottleneck_type.load()) << "\n";
        
        log_file_.flush();
    }
    
    void disable_logging() {
        std::lock_guard<std::mutex> lock(log_mutex_);
        if (log_file_.is_open()) {
            log_file_.close();
        }
        logging_enabled_ = false;
    }
    
    ~PerformanceLogger() {
        disable_logging();
    }
};

static PerformanceLogger g_performance_logger;

// Export C interface for FFI compatibility
extern "C" {
    struct RuntimeOptimizerC;
    
    RuntimeOptimizerC* runtime_optimizer_create() {
        try {
            auto& engine = StandardOptimizationEngine::get_instance();
            return reinterpret_cast<RuntimeOptimizerC*>(&engine);
        } catch (...) {
            return nullptr;
        }
    }
    
    void runtime_optimizer_start_frame(RuntimeOptimizerC* optimizer) {
        if (optimizer) {
            auto* engine = reinterpret_cast<StandardOptimizationEngine*>(optimizer);
            engine->start_frame();
        }
    }
    
    void runtime_optimizer_end_frame(RuntimeOptimizerC* optimizer, double frame_time_ms) {
        if (optimizer) {
            auto* engine = reinterpret_cast<StandardOptimizationEngine*>(optimizer);
            engine->end_frame(frame_time_ms);
        }
    }
    
    void runtime_optimizer_get_metrics(RuntimeOptimizerC* optimizer, 
                                      double* fps, double* frame_time, 
                                      size_t* memory_usage, size_t* object_count) {
        if (optimizer && fps && frame_time && memory_usage && object_count) {
            auto* engine = reinterpret_cast<StandardOptimizationEngine*>(optimizer);
            const auto metrics = engine->get_performance_metrics();
            
            *fps = metrics.fps.load();
            *frame_time = metrics.frame_time_ms.load();
            *memory_usage = metrics.memory_usage_bytes.load();
            *object_count = metrics.object_count.load();
        }
    }
    
    void runtime_optimizer_set_quality_level(RuntimeOptimizerC* optimizer, size_t level) {
        if (optimizer) {
            auto* engine = reinterpret_cast<StandardOptimizationEngine*>(optimizer);
            engine->set_quality_level(level);
        }
    }
    
    void runtime_optimizer_enable_optimization(RuntimeOptimizerC* optimizer, bool enabled) {
        if (optimizer) {
            auto* engine = reinterpret_cast<StandardOptimizationEngine*>(optimizer);
            engine->enable_optimization(enabled);
        }
    }
    
    const char* runtime_optimizer_get_report(RuntimeOptimizerC* optimizer) {
        if (optimizer) {
            auto* engine = reinterpret_cast<StandardOptimizationEngine*>(optimizer);
            static thread_local std::string report;
            report = engine->generate_performance_report();
            return report.c_str();
        }
        return nullptr;
    }
}

} // namespace hsml::core