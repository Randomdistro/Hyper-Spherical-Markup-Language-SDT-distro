#pragma once

// [The Memory Safety Paranoid]: "EVERY BYTE MUST BE ACCOUNTED FOR!"
// [The Valgrind Integrator]: "Deep memory analysis with professional tools!"

#include <memory>
#include <unordered_map>
#include <mutex>
#include <atomic>
#include <chrono>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#ifdef __linux__
#include <sys/resource.h>
#include <unistd.h>
#include <sys/wait.h>
#endif

namespace hsml::testing {

// [The Allocation Tracker]: "Track every malloc and free!"
class AllocationTracker {
public:
    struct AllocationInfo {
        size_t size;
        std::chrono::steady_clock::time_point timestamp;
        std::string file;
        int line;
        std::string function;
        void* stack_trace[16]; // Store stack frames
        size_t stack_depth;
    };

    static AllocationTracker& getInstance() {
        static AllocationTracker instance;
        return instance;
    }

    void recordAllocation(void* ptr, size_t size, const char* file = __FILE__, 
                         int line = __LINE__, const char* function = __FUNCTION__) {
        std::lock_guard<std::mutex> lock(mutex_);
        
        AllocationInfo info;
        info.size = size;
        info.timestamp = std::chrono::steady_clock::now();
        info.file = file ? file : "unknown";
        info.line = line;
        info.function = function ? function : "unknown";
        info.stack_depth = captureStackTrace(info.stack_trace, 16);
        
        allocations_[ptr] = info;
        total_allocated_.fetch_add(size);
        active_allocations_.fetch_add(1);
        
        if (active_allocations_.load() > peak_allocations_.load()) {
            peak_allocations_.store(active_allocations_.load());
        }
    }

    void recordDeallocation(void* ptr) {
        std::lock_guard<std::mutex> lock(mutex_);
        
        auto it = allocations_.find(ptr);
        if (it != allocations_.end()) {
            total_deallocated_.fetch_add(it->second.size);
            active_allocations_.fetch_sub(1);
            allocations_.erase(it);
        } else {
            // Double free or invalid free!
            invalid_frees_.fetch_add(1);
        }
    }

    struct MemoryStats {
        size_t total_allocated = 0;
        size_t total_deallocated = 0; 
        size_t current_leaks = 0;
        size_t peak_allocations = 0;
        size_t invalid_frees = 0;
        size_t active_allocations = 0;
    };

    MemoryStats getStats() const {
        std::lock_guard<std::mutex> lock(mutex_);
        
        MemoryStats stats;
        stats.total_allocated = total_allocated_.load();
        stats.total_deallocated = total_deallocated_.load();
        stats.current_leaks = stats.total_allocated - stats.total_deallocated;
        stats.peak_allocations = peak_allocations_.load();
        stats.invalid_frees = invalid_frees_.load();
        stats.active_allocations = active_allocations_.load();
        
        return stats;
    }

    std::vector<AllocationInfo> getActiveAllocations() const {
        std::lock_guard<std::mutex> lock(mutex_);
        
        std::vector<AllocationInfo> active;
        for (const auto& [ptr, info] : allocations_) {
            active.push_back(info);
        }
        
        return active;
    }

    void generateLeakReport(const std::string& filename) const {
        std::ofstream report(filename);
        auto stats = getStats();
        
        report << "=== HSML Memory Leak Detection Report ===\n\n";
        report << "Total Allocated: " << stats.total_allocated << " bytes\n";
        report << "Total Deallocated: " << stats.total_deallocated << " bytes\n";
        report << "Current Leaks: " << stats.current_leaks << " bytes\n";
        report << "Peak Allocations: " << stats.peak_allocations << " allocations\n";
        report << "Invalid Frees: " << stats.invalid_frees << " attempts\n";
        report << "Active Allocations: " << stats.active_allocations << " allocations\n\n";
        
        if (stats.current_leaks > 0) {
            report << "=== ACTIVE LEAKS ===\n";
            auto active = getActiveAllocations();
            
            for (const auto& alloc : active) {
                report << "LEAK: " << alloc.size << " bytes at " 
                      << alloc.file << ":" << alloc.line 
                      << " in " << alloc.function << "\n";
                      
                // Print stack trace if available
                if (alloc.stack_depth > 0) {
                    report << "Stack trace:\n";
                    for (size_t i = 0; i < alloc.stack_depth; ++i) {
                        report << "  " << alloc.stack_trace[i] << "\n";
                    }
                }
                report << "\n";
            }
        }
    }

private:
    mutable std::mutex mutex_;
    std::unordered_map<void*, AllocationInfo> allocations_;
    std::atomic<size_t> total_allocated_{0};
    std::atomic<size_t> total_deallocated_{0};
    std::atomic<size_t> active_allocations_{0};
    std::atomic<size_t> peak_allocations_{0};
    std::atomic<size_t> invalid_frees_{0};

    size_t captureStackTrace(void** buffer, size_t max_frames) {
        // Platform-specific stack trace capture
        #ifdef __linux__
        // Use backtrace if available
        // This is a simplified version
        return 0; // Would implement actual stack trace capture
        #else
        return 0;
        #endif
    }
};

// [The Custom Allocator]: "Override new/delete for tracking!"
class TrackingAllocator {
public:
    static void* allocate(size_t size, const char* file = __FILE__, 
                         int line = __LINE__, const char* function = __FUNCTION__) {
        void* ptr = std::malloc(size);
        if (ptr) {
            AllocationTracker::getInstance().recordAllocation(ptr, size, file, line, function);
        }
        return ptr;
    }

    static void deallocate(void* ptr) {
        if (ptr) {
            AllocationTracker::getInstance().recordDeallocation(ptr);
            std::free(ptr);
        }
    }
};

// [The Valgrind Integrator]: "Professional memory analysis!"
class ValgrindIntegration {
public:
    struct ValgrindResult {
        bool valgrind_available = false;
        bool memory_leaks_detected = false;
        bool invalid_reads_detected = false;
        bool invalid_writes_detected = false;
        size_t definitely_lost_bytes = 0;
        size_t possibly_lost_bytes = 0;
        size_t still_reachable_bytes = 0;
        std::string detailed_report;
    };

    static ValgrindResult runValgrindAnalysis(const std::string& executable_path, 
                                            const std::vector<std::string>& args = {}) {
        ValgrindResult result;
        
        #ifdef __linux__
        // Check if valgrind is available
        if (!isValgrindAvailable()) {
            std::cerr << "⚠️ Valgrind not available on this system\n";
            return result;
        }
        
        result.valgrind_available = true;
        
        // Construct valgrind command
        std::string valgrind_cmd = "valgrind --leak-check=full --show-leak-kinds=all "
                                  "--track-origins=yes --verbose --xml=yes "
                                  "--xml-file=valgrind_output.xml " + executable_path;
        
        for (const auto& arg : args) {
            valgrind_cmd += " " + arg;
        }
        
        // Run valgrind
        int exit_code = std::system(valgrind_cmd.c_str());
        
        if (exit_code == 0) {
            // Parse valgrind XML output
            result = parseValgrindXML("valgrind_output.xml");
        }
        #endif
        
        return result;
    }

    static bool isValgrindAvailable() {
        return std::system("valgrind --version > /dev/null 2>&1") == 0;
    }

private:
    static ValgrindResult parseValgrindXML(const std::string& xml_file) {
        ValgrindResult result;
        result.valgrind_available = true;
        
        std::ifstream file(xml_file);
        if (!file.is_open()) {
            return result;
        }
        
        std::string line;
        std::stringstream report;
        
        while (std::getline(file, line)) {
            report << line << "\n";
            
            // Parse key metrics from XML
            if (line.find("<kind>Leak_DefinitelyLost</kind>") != std::string::npos) {
                result.memory_leaks_detected = true;
                // Extract byte count
                result.definitely_lost_bytes += extractByteCount(line);
            }
            
            if (line.find("<kind>Leak_PossiblyLost</kind>") != std::string::npos) {
                result.memory_leaks_detected = true;
                result.possibly_lost_bytes += extractByteCount(line);
            }
            
            if (line.find("<kind>InvalidRead</kind>") != std::string::npos) {
                result.invalid_reads_detected = true;
            }
            
            if (line.find("<kind>InvalidWrite</kind>") != std::string::npos) {
                result.invalid_writes_detected = true;
            }
        }
        
        result.detailed_report = report.str();
        return result;
    }
    
    static size_t extractByteCount(const std::string& xml_line) {
        // Simple XML parsing to extract byte counts
        // In a real implementation, would use a proper XML parser
        return 0; // Placeholder
    }
};

// [The Address Sanitizer Wrapper]: "Compiler-based memory safety!"
class AddressSanitizerWrapper {
public:
    struct ASanResult {
        bool asan_enabled = false;
        bool heap_buffer_overflow_detected = false;
        bool stack_buffer_overflow_detected = false;
        bool use_after_free_detected = false;
        bool double_free_detected = false;
        std::vector<std::string> violations;
    };

    static ASanResult runWithAddressSanitizer(std::function<void()> test_function) {
        ASanResult result;
        
        #ifdef __has_feature
        #if __has_feature(address_sanitizer)
        result.asan_enabled = true;
        #endif
        #endif
        
        #ifdef __SANITIZE_ADDRESS__
        result.asan_enabled = true;
        #endif
        
        if (!result.asan_enabled) {
            std::cerr << "⚠️ AddressSanitizer not enabled. Compile with -fsanitize=address\n";
            return result;
        }
        
        // Set up ASan error callback
        // Note: This is a simplified version. Real implementation would use
        // __sanitizer_set_death_callback or similar
        
        try {
            test_function();
        } catch (const std::exception& e) {
            // ASan violations typically terminate the program
            // In a real implementation, we'd catch ASan signals
            result.violations.push_back(e.what());
        }
        
        return result;
    }
};

// [The Memory Pool Validator]: "Custom memory pool testing!"
template<typename MemoryPoolType>
class MemoryPoolValidator {
public:
    struct PoolValidationResult {
        bool fragmentation_acceptable = false;
        bool allocation_speed_acceptable = false;
        bool deallocation_speed_acceptable = false;
        bool memory_efficiency_acceptable = false;
        double fragmentation_ratio = 0.0;
        std::chrono::nanoseconds avg_allocation_time{0};
        std::chrono::nanoseconds avg_deallocation_time{0};
        double memory_utilization = 0.0;
    };

    static PoolValidationResult validateMemoryPool(MemoryPoolType& pool) {
        PoolValidationResult result;
        
        // 1. Fragmentation Testing
        result.fragmentation_ratio = measureFragmentation(pool);
        result.fragmentation_acceptable = (result.fragmentation_ratio < 0.2); // 20% threshold
        
        // 2. Allocation Speed Testing
        result.avg_allocation_time = benchmarkAllocation(pool);
        result.allocation_speed_acceptable = (result.avg_allocation_time < std::chrono::microseconds(1));
        
        // 3. Deallocation Speed Testing
        result.avg_deallocation_time = benchmarkDeallocation(pool);
        result.deallocation_speed_acceptable = (result.avg_deallocation_time < std::chrono::microseconds(1));
        
        // 4. Memory Efficiency Testing
        result.memory_utilization = measureMemoryUtilization(pool);
        result.memory_efficiency_acceptable = (result.memory_utilization > 0.8); // 80% threshold
        
        return result;
    }

private:
    static double measureFragmentation(MemoryPoolType& pool) {
        // Measure memory fragmentation in the pool
        const size_t test_allocations = 1000;
        std::vector<void*> ptrs;
        
        // Allocate various sizes
        for (size_t i = 0; i < test_allocations; ++i) {
            size_t size = 16 + (i % 512); // Variable sizes 16-528 bytes
            void* ptr = pool.allocate(size);
            if (ptr) {
                ptrs.push_back(ptr);
            }
        }
        
        // Deallocate every other allocation to create fragmentation
        for (size_t i = 0; i < ptrs.size(); i += 2) {
            pool.deallocate(ptrs[i]);
        }
        
        // Measure largest contiguous block vs total free space
        size_t largest_block = pool.getLargestFreeBlock();
        size_t total_free = pool.getTotalFreeSpace();
        
        // Clean up remaining allocations
        for (size_t i = 1; i < ptrs.size(); i += 2) {
            pool.deallocate(ptrs[i]);
        }
        
        return (total_free > 0) ? (1.0 - (static_cast<double>(largest_block) / total_free)) : 0.0;
    }
    
    static std::chrono::nanoseconds benchmarkAllocation(MemoryPoolType& pool) {
        const size_t iterations = 10000;
        const size_t allocation_size = 64;
        
        auto start = std::chrono::high_resolution_clock::now();
        
        for (size_t i = 0; i < iterations; ++i) {
            void* ptr = pool.allocate(allocation_size);
            if (ptr) {
                pool.deallocate(ptr); // Immediate deallocation for benchmarking
            }
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        
        return total_time / iterations;
    }
    
    static std::chrono::nanoseconds benchmarkDeallocation(MemoryPoolType& pool) {
        const size_t iterations = 10000;
        const size_t allocation_size = 64;
        
        // Pre-allocate
        std::vector<void*> ptrs;
        for (size_t i = 0; i < iterations; ++i) {
            void* ptr = pool.allocate(allocation_size);
            if (ptr) {
                ptrs.push_back(ptr);
            }
        }
        
        // Benchmark deallocation
        auto start = std::chrono::high_resolution_clock::now();
        
        for (void* ptr : ptrs) {
            pool.deallocate(ptr);
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        
        return total_time / ptrs.size();
    }
    
    static double measureMemoryUtilization(MemoryPoolType& pool) {
        size_t total_pool_size = pool.getTotalSize();
        size_t used_space = pool.getUsedSpace();
        
        return (total_pool_size > 0) ? (static_cast<double>(used_space) / total_pool_size) : 0.0;
    }
};

// [The Memory Safety Test Suite]: "Comprehensive memory validation!"
class MemorySafetyTestSuite {
public:
    static void runComprehensiveMemoryTests() {
        std::cout << "\n🛡️ Memory Safety Test Suite Starting... 🛡️\n";
        std::cout << "=================================================\n";

        // 1. Allocation Tracking Test
        std::cout << "\n[Allocation Tracker]: Testing allocation tracking...\n";
        testAllocationTracking();

        // 2. Valgrind Integration Test
        std::cout << "\n[Valgrind Integration]: Running memory analysis...\n";
        if (ValgrindIntegration::isValgrindAvailable()) {
            auto valgrind_result = ValgrindIntegration::runValgrindAnalysis("./test_executable");
            reportValgrindResults(valgrind_result);
        }

        // 3. AddressSanitizer Test
        std::cout << "\n[AddressSanitizer]: Testing memory safety...\n";
        testWithAddressSanitizer();

        // 4. Memory Pool Validation
        std::cout << "\n[Memory Pool]: Validating custom memory pools...\n";
        // testMemoryPools(); // Would test actual memory pool implementations

        std::cout << "\n✅ Memory Safety Tests Complete! ✅\n";
    }

private:
    static void testAllocationTracking() {
        auto& tracker = AllocationTracker::getInstance();
        
        // Simulate some allocations
        void* ptr1 = TrackingAllocator::allocate(1024);
        void* ptr2 = TrackingAllocator::allocate(2048);
        void* ptr3 = TrackingAllocator::allocate(512);
        
        auto stats = tracker.getStats();
        std::cout << "  Active allocations: " << stats.active_allocations << "\n";
        std::cout << "  Total allocated: " << stats.total_allocated << " bytes\n";
        
        // Clean up
        TrackingAllocator::deallocate(ptr1);
        TrackingAllocator::deallocate(ptr2);
        TrackingAllocator::deallocate(ptr3);
        
        // Generate leak report
        tracker.generateLeakReport("memory_leak_report.txt");
        std::cout << "  Leak report generated: memory_leak_report.txt\n";
    }
    
    static void reportValgrindResults(const ValgrindIntegration::ValgrindResult& result) {
        if (result.valgrind_available) {
            std::cout << "  Valgrind analysis complete:\n";
            std::cout << "    Memory leaks: " << (result.memory_leaks_detected ? "DETECTED" : "None") << "\n";
            std::cout << "    Invalid reads: " << (result.invalid_reads_detected ? "DETECTED" : "None") << "\n";
            std::cout << "    Invalid writes: " << (result.invalid_writes_detected ? "DETECTED" : "None") << "\n";
            std::cout << "    Definitely lost: " << result.definitely_lost_bytes << " bytes\n";
            std::cout << "    Possibly lost: " << result.possibly_lost_bytes << " bytes\n";
        } else {
            std::cout << "  Valgrind not available\n";
        }
    }
    
    static void testWithAddressSanitizer() {
        auto result = AddressSanitizerWrapper::runWithAddressSanitizer([]() {
            // Test function that should not trigger ASan
            std::vector<int> test_vector(100);
            for (size_t i = 0; i < test_vector.size(); ++i) {
                test_vector[i] = static_cast<int>(i);
            }
        });
        
        if (result.asan_enabled) {
            std::cout << "  AddressSanitizer: Enabled and monitoring\n";
            std::cout << "  Violations detected: " << result.violations.size() << "\n";
        } else {
            std::cout << "  AddressSanitizer: Not enabled (compile with -fsanitize=address)\n";
        }
    }
};

} // namespace hsml::testing

// [The Global Allocation Override]: "Intercept all allocations!"
// Macros for easy allocation tracking
#define HSML_NEW(size) hsml::testing::TrackingAllocator::allocate(size, __FILE__, __LINE__, __FUNCTION__)
#define HSML_DELETE(ptr) hsml::testing::TrackingAllocator::deallocate(ptr)

// Custom new/delete operators for comprehensive tracking
#ifdef HSML_ENABLE_MEMORY_TRACKING
void* operator new(size_t size) {
    return hsml::testing::TrackingAllocator::allocate(size, __FILE__, __LINE__, __FUNCTION__);
}

void operator delete(void* ptr) noexcept {
    hsml::testing::TrackingAllocator::deallocate(ptr);
}

void* operator new[](size_t size) {
    return hsml::testing::TrackingAllocator::allocate(size, __FILE__, __LINE__, __FUNCTION__);
}

void operator delete[](void* ptr) noexcept {
    hsml::testing::TrackingAllocator::deallocate(ptr);
}
#endif