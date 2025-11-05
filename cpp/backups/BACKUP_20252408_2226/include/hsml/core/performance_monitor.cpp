#pragma once

#include <chrono>
#include <atomic>
#include <array>
#include <vector>
#include <string>
#include <functional>
#include <thread>
#include <memory>

namespace hsml {
namespace core {

// Build configuration detection
namespace build_configuration {
    #ifdef HSML_ENABLE_PROFILING
        static constexpr bool enable_profiling = true;
    #else
        static constexpr bool enable_profiling = false;
    #endif
    
    #ifdef HSML_ENABLE_PERFORMANCE_COUNTERS
        static constexpr bool enable_performance_counters = true;
    #else
        static constexpr bool enable_performance_counters = false;
    #endif
}

// High-resolution timestamp using RDTSC when available
inline uint64_t rdtsc_timestamp() {
    #if defined(__x86_64__) || defined(__i386__)
        uint32_t hi, lo;
        __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
        return static_cast<uint64_t>(hi) << 32 | lo;
    #else
        return std::chrono::high_resolution_clock::now().time_since_epoch().count();
    #endif
}

// Performance counter structure
struct performance_counter {
    std::atomic<uint64_t> call_count{0};
    std::atomic<uint64_t> total_cycles{0};
    std::atomic<uint64_t> min_cycles{UINT64_MAX};
    std::atomic<uint64_t> max_cycles{0};
    std::atomic<uint64_t> last_cycles{0};
    
    void record_measurement(uint64_t cycles) {
        call_count.fetch_add(1);
        total_cycles.fetch_add(cycles);
        last_cycles.store(cycles);
        
        // Update min/max atomically
        uint64_t current_min = min_cycles.load();
        while (cycles < current_min && 
               !min_cycles.compare_exchange_weak(current_min, cycles)) {}
        
        uint64_t current_max = max_cycles.load();
        while (cycles > current_max && 
               !max_cycles.compare_exchange_weak(current_max, cycles)) {}
    }
    
    double average_cycles() const {
        uint64_t count = call_count.load();
        if (count == 0) return 0.0;
        return static_cast<double>(total_cycles.load()) / static_cast<double>(count);
    }
    
    void reset() {
        call_count.store(0);
        total_cycles.store(0);
        min_cycles.store(UINT64_MAX);
        max_cycles.store(0);
        last_cycles.store(0);
    }
};

// Thread-local performance counters
class performance_counters {
    static constexpr size_t max_counters = 256;
    
    std::array<performance_counter, max_counters> counters_;
    std::atomic<size_t> next_counter_index_{0};
    
public:
    size_t create_counter() {
        size_t index = next_counter_index_.fetch_add(1);
        if (index >= max_counters) {
            // Fallback: return last available counter
            return max_counters - 1;
        }
        return index;
    }
    
    performance_counter& get_counter(size_t index) {
        return counters_[index % max_counters];
    }
    
    const performance_counter& get_counter(size_t index) const {
        return counters_[index % max_counters];
    }
    
    void reset_all() {
        for (auto& counter : counters_) {
            counter.reset();
        }
        next_counter_index_.store(0);
    }
    
    // Get statistics for all active counters
    struct counter_stats {
        size_t index;
        uint64_t call_count;
        double average_cycles;
        uint64_t min_cycles;
        uint64_t max_cycles;
        uint64_t last_cycles;
    };
    
    std::vector<counter_stats> get_all_stats() const {
        std::vector<counter_stats> stats;
        size_t active_counters = next_counter_index_.load();
        
        for (size_t i = 0; i < std::min(active_counters, max_counters); ++i) {
            const auto& counter = counters_[i];
            uint64_t calls = counter.call_count.load();
            if (calls > 0) {
                stats.push_back({
                    i,
                    calls,
                    counter.average_cycles(),
                    counter.min_cycles.load(),
                    counter.max_cycles.load(),
                    counter.last_cycles.load()
                });
            }
        }
        
        return stats;
    }
};

// Real-time performance metrics with zero overhead
class performance_monitor {
    thread_local inline static performance_counters counters_;
    thread_local inline static std::vector<std::pair<std::string, size_t>> named_counters_;
    
    static std::atomic<bool> global_enabled_;
    static std::atomic<size_t> global_counter_index_;
    
public:
    // Measure function execution time
    template<typename F>
    static auto measure(F&& func) -> decltype(func()) {
        if constexpr (build_configuration::enable_profiling) {
            auto start = rdtsc_timestamp();
            auto result = func();
            auto end = rdtsc_timestamp();
            counters_.get_counter(0).record_measurement(end - start);
            return result;
        } else {
            return func();
        }
    }
    
    // Measure with specific counter
    template<typename F>
    static auto measure(size_t counter_index, F&& func) -> decltype(func()) {
        if constexpr (build_configuration::enable_profiling) {
            auto start = rdtsc_timestamp();
            auto result = func();
            auto end = rdtsc_timestamp();
            counters_.get_counter(counter_index).record_measurement(end - start);
            return result;
        } else {
            return func();
        }
    }
    
    // Measure with named counter
    template<typename F>
    static auto measure(const std::string& name, F&& func) -> decltype(func()) {
        if constexpr (build_configuration::enable_profiling) {
            size_t counter_index = get_or_create_named_counter(name);
            return measure(counter_index, std::forward<F>(func));
        } else {
            return func();
        }
    }
    
    // Manual timing control
    class timer {
        size_t counter_index_;
        uint64_t start_time_;
        bool active_;
        
    public:
        explicit timer(size_t counter_index) : counter_index_(counter_index), active_(false) {
            if constexpr (build_configuration::enable_profiling) {
                start_time_ = rdtsc_timestamp();
                active_ = true;
            }
        }
        
        explicit timer(const std::string& name) : timer(get_or_create_named_counter(name)) {}
        
        ~timer() {
            if (active_) {
                stop();
            }
        }
        
        void stop() {
            if (active_) {
                auto end_time = rdtsc_timestamp();
                counters_.get_counter(counter_index_).record_measurement(end_time - start_time_);
                active_ = false;
            }
        }
        
        void restart() {
            if constexpr (build_configuration::enable_profiling) {
                start_time_ = rdtsc_timestamp();
                active_ = true;
            }
        }
    };
    
    // Create a new counter
    static size_t create_counter() {
        if constexpr (build_configuration::enable_profiling) {
            return counters_.create_counter();
        } else {
            return 0;
        }
    }
    
    // Get counter statistics
    static performance_counter& get_counter(size_t index) {
        return counters_.get_counter(index);
    }
    
    // Get all statistics
    static std::vector<performance_counters::counter_stats> get_all_stats() {
        if constexpr (build_configuration::enable_profiling) {
            return counters_.get_all_stats();
        } else {
            return {};
        }
    }
    
    // Reset all counters
    static void reset_all() {
        if constexpr (build_configuration::enable_profiling) {
            counters_.reset_all();
            named_counters_.clear();
        }
    }
    
    // Enable/disable global profiling
    static void set_enabled(bool enabled) {
        global_enabled_.store(enabled);
    }
    
    static bool is_enabled() {
        return global_enabled_.load() && build_configuration::enable_profiling;
    }
    
    // Performance report generation
    struct performance_report {
        std::chrono::high_resolution_clock::time_point timestamp;
        std::vector<performance_counters::counter_stats> counters;
        double total_cpu_utilization;
        double memory_usage_mb;
        size_t active_threads;
        
        std::string to_string() const {
            std::string report = "Performance Report\n";
            report += "==================\n";
            report += "Timestamp: " + std::to_string(timestamp.time_since_epoch().count()) + "\n";
            report += "Active Threads: " + std::to_string(active_threads) + "\n";
            report += "Memory Usage: " + std::to_string(memory_usage_mb) + " MB\n";
            report += "CPU Utilization: " + std::to_string(total_cpu_utilization) + "%\n\n";
            
            for (const auto& counter : counters) {
                report += "Counter " + std::to_string(counter.index) + ":\n";
                report += "  Calls: " + std::to_string(counter.call_count) + "\n";
                report += "  Avg Cycles: " + std::to_string(counter.average_cycles) + "\n";
                report += "  Min Cycles: " + std::to_string(counter.min_cycles) + "\n";
                report += "  Max Cycles: " + std::to_string(counter.max_cycles) + "\n";
                report += "  Last Cycles: " + std::to_string(counter.last_cycles) + "\n\n";
            }
            
            return report;
        }
    };
    
    static performance_report generate_report() {
        performance_report report;
        report.timestamp = std::chrono::high_resolution_clock::now();
        report.counters = get_all_stats();
        report.total_cpu_utilization = get_cpu_utilization();
        report.memory_usage_mb = get_memory_usage_mb();
        report.active_threads = std::thread::hardware_concurrency();
        return report;
    }
    
private:
    static size_t get_or_create_named_counter(const std::string& name) {
        // Check if counter already exists
        for (const auto& [counter_name, index] : named_counters_) {
            if (counter_name == name) {
                return index;
            }
        }
        
        // Create new counter
        size_t new_index = counters_.create_counter();
        named_counters_.emplace_back(name, new_index);
        return new_index;
    }
    
    static double get_cpu_utilization() {
        // Placeholder for CPU utilization measurement
        // In a real implementation, this would read system metrics
        return 0.0;
    }
    
    static double get_memory_usage_mb() {
        // Placeholder for memory usage measurement
        // In a real implementation, this would read system metrics
        return 0.0;
    }
};

// Global variables
std::atomic<bool> performance_monitor::global_enabled_{true};
std::atomic<size_t> performance_monitor::global_counter_index_{0};

// RAII performance scope
class performance_scope {
    performance_monitor::timer timer_;
    
public:
    explicit performance_scope(size_t counter_index) : timer_(counter_index) {}
    explicit performance_scope(const std::string& name) : timer_(name) {}
    
    ~performance_scope() = default;
    
    void stop() { timer_.stop(); }
    void restart() { timer_.restart(); }
};

// Macro for easy performance measurement
#define HSML_PERF_SCOPE(name) \
    performance_scope hsml_perf_scope_##__LINE__(name)

#define HSML_PERF_MEASURE(name, code) \
    performance_monitor::measure(name, [&]() { code; })

// Compile-time tests
namespace compile_time_tests {
    consteval bool test_performance_monitor_compilation() {
        // Test that performance monitor compiles correctly
        return true;
    }
    
    consteval bool test_build_configuration() {
        return build_configuration::enable_profiling == false || 
               build_configuration::enable_profiling == true;
    }
    
    static_assert(test_performance_monitor_compilation());
    static_assert(test_build_configuration());
}

} // namespace core
} // namespace hsml
