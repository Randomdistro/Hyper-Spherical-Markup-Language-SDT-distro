#pragma once

#include "../core/advanced_concepts.h"
#include "../core/constexpr_spherical_coords.h"
#include "gpu_memory_manager.h"
#include <future>
#include <span>
#include <vector>
#include <thread>
#include <atomic>
#include <memory>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <type_traits>

namespace hsml {
namespace compute {

using namespace core;
using namespace core::concepts;

// Compute device capabilities
struct compute_capabilities {
    bool has_simd = false;
    bool has_avx2 = false;
    bool has_avx512 = false;
    bool has_gpu = false;
    bool has_opencl = false;
    bool has_cuda = false;
    bool has_vulkan_compute = false;
    
    size_t cpu_cores = 0;
    size_t gpu_compute_units = 0;
    size_t gpu_memory_mb = 0;
    size_t cpu_cache_size_kb = 0;
    
    double cpu_base_frequency_ghz = 0.0;
    double gpu_base_frequency_ghz = 0.0;
    
    [[nodiscard]] double get_cpu_performance_score() const noexcept {
        double score = cpu_cores * cpu_base_frequency_ghz;
        if (has_avx512) score *= 2.0;
        else if (has_avx2) score *= 1.5;
        else if (has_simd) score *= 1.2;
        return score;
    }
    
    [[nodiscard]] double get_gpu_performance_score() const noexcept {
        if (!has_gpu) return 0.0;
        return gpu_compute_units * gpu_base_frequency_ghz * 0.1;  // Scaled for comparison
    }
    
    [[nodiscard]] bool is_gpu_preferred() const noexcept {
        return has_gpu && (get_gpu_performance_score() > get_cpu_performance_score());
    }
};

// Runtime device detection
class compute_device_detector {
private:
    static compute_capabilities detect_cpu_capabilities() {
        compute_capabilities caps;
        
        caps.cpu_cores = std::thread::hardware_concurrency();
        caps.cpu_base_frequency_ghz = 2.5;  // Default estimate
        caps.cpu_cache_size_kb = 32768;     // Default estimate
        
#if defined(__x86_64__) || defined(_M_X64)
        // Check CPU features using CPUID
        caps.has_simd = true;  // SSE2 is guaranteed on x86_64
        
#ifdef __AVX2__
        caps.has_avx2 = true;
#endif

#ifdef __AVX512F__
        caps.has_avx512 = true;
#endif
        
#elif defined(__ARM_NEON) || defined(__ARM_NEON__)
        caps.has_simd = true;
#endif

        return caps;
    }
    
    static compute_capabilities detect_gpu_capabilities() {
        compute_capabilities caps;
        
        // Simplified GPU detection - in real implementation would use APIs
        caps.has_gpu = false;           // Conservative default
        caps.has_opencl = false;        // Would check for OpenCL runtime
        caps.has_cuda = false;          // Would check for CUDA runtime
        caps.has_vulkan_compute = false; // Would check for Vulkan
        
        caps.gpu_compute_units = 0;
        caps.gpu_memory_mb = 0;
        caps.gpu_base_frequency_ghz = 0.0;
        
        return caps;
    }
    
public:
    [[nodiscard]] static compute_capabilities detect_all_capabilities() {
        auto cpu_caps = detect_cpu_capabilities();
        auto gpu_caps = detect_gpu_capabilities();
        
        // Merge capabilities
        compute_capabilities merged = cpu_caps;
        merged.has_gpu = gpu_caps.has_gpu;
        merged.has_opencl = gpu_caps.has_opencl;
        merged.has_cuda = gpu_caps.has_cuda;
        merged.has_vulkan_compute = gpu_caps.has_vulkan_compute;
        merged.gpu_compute_units = gpu_caps.gpu_compute_units;
        merged.gpu_memory_mb = gpu_caps.gpu_memory_mb;
        merged.gpu_base_frequency_ghz = gpu_caps.gpu_base_frequency_ghz;
        
        return merged;
    }
};

// Task types for different compute patterns
enum class compute_task_type {
    spherical_transform,    // Coordinate transformations
    solid_angle_calc,      // Solid angle calculations
    matrix_multiply,       // Matrix operations
    vector_operations,     // Vector math
    pixel_processing,      // Pixel shader equivalent
    reduction_ops,         // Sum, min, max, etc.
    convolution,          // Filtering operations
    custom_kernel         // User-defined compute
};

// Compute task descriptor
template<typename InputType, typename OutputType>
struct compute_task {
    compute_task_type type;
    std::span<const InputType> inputs;
    std::span<OutputType> outputs;
    std::function<void(std::span<const InputType>, std::span<OutputType>)> cpu_kernel;
    std::string gpu_kernel_source;  // Shader/kernel source code
    size_t local_work_size = 64;    // GPU workgroup size
    
    compute_task(compute_task_type t, 
                std::span<const InputType> in, 
                std::span<OutputType> out)
        : type(t), inputs(in), outputs(out) {}
};

// Performance measurement utilities
struct performance_measurement {
    std::chrono::nanoseconds execution_time{0};
    size_t elements_processed = 0;
    compute_backend backend_used = compute_backend::auto_detect;
    
    [[nodiscard]] double get_throughput_gops() const noexcept {
        if (execution_time.count() == 0) return 0.0;
        
        double operations_per_sec = static_cast<double>(elements_processed) * 1e9 / execution_time.count();
        return operations_per_sec / 1e9;  // Convert to GOPS
    }
    
    [[nodiscard]] double get_latency_us() const noexcept {
        return static_cast<double>(execution_time.count()) / 1000.0;
    }
};

// Automatic work distribution engine
template<compute_backend Backend = compute_backend::auto_detect>
class heterogeneous_compute_engine {
private:
    compute_capabilities capabilities_;
    std::unique_ptr<vulkan_memory_manager> gpu_memory_manager_;
    
    // Performance history for adaptive scheduling
    mutable std::vector<performance_measurement> cpu_history_;
    mutable std::vector<performance_measurement> gpu_history_;
    mutable std::mutex history_mutex_;
    
    // Thresholds for backend selection
    static constexpr size_t MIN_GPU_WORKLOAD = 1024;
    static constexpr double GPU_OVERHEAD_FACTOR = 1.5;
    
    template<typename InputType, typename OutputType>
    [[nodiscard]] compute_backend select_optimal_backend(
        const compute_task<InputType, OutputType>& task) const {
        
        if constexpr (Backend != compute_backend::auto_detect) {
            return Backend;
        }
        
        size_t workload_size = task.inputs.size();
        
        // Small workloads prefer CPU
        if (workload_size < MIN_GPU_WORKLOAD) {
            return capabilities_.has_simd ? compute_backend::simd : compute_backend::cpu;
        }
        
        // Large workloads prefer GPU if available
        if (capabilities_.has_gpu && workload_size > MIN_GPU_WORKLOAD * 4) {
            return compute_backend::gpu;
        }
        
        // Use performance history for adaptive selection
        std::lock_guard<std::mutex> lock(history_mutex_);
        
        if (!cpu_history_.empty() && !gpu_history_.empty()) {
            auto avg_cpu_throughput = std::accumulate(cpu_history_.begin(), cpu_history_.end(), 0.0,
                [](double sum, const performance_measurement& pm) {
                    return sum + pm.get_throughput_gops();
                }) / cpu_history_.size();
            
            auto avg_gpu_throughput = std::accumulate(gpu_history_.begin(), gpu_history_.end(), 0.0,
                [](double sum, const performance_measurement& pm) {
                    return sum + pm.get_throughput_gops();
                }) / gpu_history_.size();
            
            // Account for GPU transfer overhead
            if (avg_gpu_throughput > avg_cpu_throughput * GPU_OVERHEAD_FACTOR) {
                return compute_backend::gpu;
            }
        }
        
        // Default to best available CPU backend
        if (capabilities_.has_avx512) return compute_backend::simd;
        if (capabilities_.has_avx2) return compute_backend::simd;
        return compute_backend::cpu;
    }
    
    template<typename InputType, typename OutputType>
    performance_measurement execute_cpu_task(const compute_task<InputType, OutputType>& task) const {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        if (task.cpu_kernel) {
            task.cpu_kernel(task.inputs, task.outputs);
        } else {
            // Default implementations for common task types
            execute_default_cpu_kernel(task);
        }
        
        auto end_time = std::chrono::high_resolution_clock::now();
        
        performance_measurement measurement;
        measurement.execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
        measurement.elements_processed = task.inputs.size();
        measurement.backend_used = compute_backend::cpu;
        
        return measurement;
    }
    
    template<typename InputType, typename OutputType>
    performance_measurement execute_simd_task(const compute_task<InputType, OutputType>& task) const {
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // SIMD implementation (simplified)
        execute_simd_kernel(task);
        
        auto end_time = std::chrono::high_resolution_clock::now();
        
        performance_measurement measurement;
        measurement.execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
        measurement.elements_processed = task.inputs.size();
        measurement.backend_used = compute_backend::simd;
        
        return measurement;
    }
    
    template<typename InputType, typename OutputType>
    performance_measurement execute_gpu_task(const compute_task<InputType, OutputType>& task) const {
        if (!capabilities_.has_gpu || !gpu_memory_manager_) {
            return execute_cpu_task(task);
        }
        
        auto start_time = std::chrono::high_resolution_clock::now();
        
        // Allocate GPU buffers
        auto input_buffer = gpu_memory_manager_->allocate_buffer<InputType>(
            task.inputs.size(), memory_usage::static_read);
        auto output_buffer = gpu_memory_manager_->allocate_buffer<OutputType>(
            task.outputs.size(), memory_usage::dynamic_write);
        
        // Copy input data to GPU
        std::copy(task.inputs.begin(), task.inputs.end(), input_buffer.begin());
        input_buffer.sync_to_gpu();
        
        // Execute GPU kernel (simplified)
        execute_gpu_kernel(task, input_buffer, output_buffer);
        
        // Copy results back
        output_buffer.sync_from_gpu();
        std::copy(output_buffer.begin(), output_buffer.end(), task.outputs.begin());
        
        auto end_time = std::chrono::high_resolution_clock::now();
        
        performance_measurement measurement;
        measurement.execution_time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time);
        measurement.elements_processed = task.inputs.size();
        measurement.backend_used = compute_backend::gpu;
        
        return measurement;
    }
    
    template<typename InputType, typename OutputType>
    void execute_default_cpu_kernel(const compute_task<InputType, OutputType>& task) const {
        switch (task.type) {
            case compute_task_type::spherical_transform:
                execute_spherical_transform_cpu(task);
                break;
            case compute_task_type::solid_angle_calc:
                execute_solid_angle_calc_cpu(task);
                break;
            case compute_task_type::vector_operations:
                execute_vector_ops_cpu(task);
                break;
            default:
                // Generic element-wise processing
                std::transform(task.inputs.begin(), task.inputs.end(), task.outputs.begin(),
                    [](const InputType& input) -> OutputType {
                        return static_cast<OutputType>(input);
                    });
        }
    }
    
    template<typename InputType, typename OutputType>
    void execute_spherical_transform_cpu(const compute_task<InputType, OutputType>& task) const {
        if constexpr (spherical_coords_like<InputType> && vector3_like<OutputType>) {
            std::transform(task.inputs.begin(), task.inputs.end(), task.outputs.begin(),
                [](const InputType& spherical) -> OutputType {
                    return spherical.template to_cartesian<OutputType>();
                });
        }
    }
    
    template<typename InputType, typename OutputType>
    void execute_solid_angle_calc_cpu(const compute_task<InputType, OutputType>& task) const {
        using T = typename InputType::value_type;
        
        if constexpr (std::is_same_v<InputType, pixel_coordinate<T>>) {
            display_geometry<T> default_geometry{1920, 1080, 650};
            
            std::transform(task.inputs.begin(), task.inputs.end(), task.outputs.begin(),
                [&](const InputType& pixel) -> OutputType {
                    return static_cast<OutputType>(
                        constexpr_solid_angle_engine<T>::calculate_pixel_steradian(pixel, default_geometry)
                    );
                });
        }
    }
    
    template<typename InputType, typename OutputType>
    void execute_vector_ops_cpu(const compute_task<InputType, OutputType>& task) const {
        if constexpr (vector3_like<InputType> && vector3_like<OutputType>) {
            std::transform(task.inputs.begin(), task.inputs.end(), task.outputs.begin(),
                [](const InputType& vec) -> OutputType {
                    return OutputType{vec.normalized()};
                });
        }
    }
    
    template<typename InputType, typename OutputType>
    void execute_simd_kernel(const compute_task<InputType, OutputType>& task) const {
        // SIMD optimized versions of the CPU kernels
        // For demonstration, fall back to CPU
        execute_default_cpu_kernel(task);
    }
    
    template<typename InputType, typename OutputType>
    void execute_gpu_kernel(const compute_task<InputType, OutputType>& task,
                           const gpu_buffer_view<InputType>& input_buffer,
                           gpu_buffer_view<OutputType>& output_buffer) const {
        // GPU kernel execution (simplified - would use actual GPU APIs)
        // For demonstration, copy data and process on CPU
        execute_default_cpu_kernel(task);
    }
    
    void update_performance_history(const performance_measurement& measurement) const {
        std::lock_guard<std::mutex> lock(history_mutex_);
        
        constexpr size_t MAX_HISTORY_SIZE = 100;
        
        if (measurement.backend_used == compute_backend::cpu || 
            measurement.backend_used == compute_backend::simd) {
            cpu_history_.push_back(measurement);
            if (cpu_history_.size() > MAX_HISTORY_SIZE) {
                cpu_history_.erase(cpu_history_.begin());
            }
        } else if (measurement.backend_used == compute_backend::gpu) {
            gpu_history_.push_back(measurement);
            if (gpu_history_.size() > MAX_HISTORY_SIZE) {
                gpu_history_.erase(gpu_history_.begin());
            }
        }
    }
    
public:
    heterogeneous_compute_engine() 
        : capabilities_(compute_device_detector::detect_all_capabilities()) {
        
        if (capabilities_.has_gpu) {
            gpu_memory_manager_ = std::make_unique<vulkan_memory_manager>();
        }
    }
    
    // Execute single compute task with automatic backend selection
    template<spherical_operation_concept Op>
    auto execute_batch(std::span<const typename Op::input_type> inputs) 
        -> std::future<std::vector<typename Op::result_type>> {
        
        using InputType = typename Op::input_type;
        using OutputType = typename Op::result_type;
        
        return std::async(std::launch::async, [this, inputs]() -> std::vector<OutputType> {
            std::vector<OutputType> outputs(inputs.size());
            
            compute_task<InputType, OutputType> task{
                compute_task_type::custom_kernel,
                inputs,
                std::span<OutputType>{outputs}
            };
            
            // Set up operation-specific kernel
            task.cpu_kernel = [](std::span<const InputType> in, std::span<OutputType> out) {
                Op operation;
                std::transform(in.begin(), in.end(), out.begin(),
                    [&operation](const InputType& input) -> OutputType {
                        return operation.execute(input);
                    });
            };
            
            auto backend = select_optimal_backend(task);
            performance_measurement measurement;
            
            switch (backend) {
                case compute_backend::cpu:
                    measurement = execute_cpu_task(task);
                    break;
                case compute_backend::simd:
                    measurement = execute_simd_task(task);
                    break;
                case compute_backend::gpu:
                    measurement = execute_gpu_task(task);
                    break;
                default:
                    measurement = execute_cpu_task(task);
                    break;
            }
            
            update_performance_history(measurement);
            return outputs;
        });
    }
    
    // Execute generic compute task
    template<typename InputType, typename OutputType>
    performance_measurement execute_task(compute_task<InputType, OutputType>& task) {
        auto backend = select_optimal_backend(task);
        performance_measurement measurement;
        
        switch (backend) {
            case compute_backend::cpu:
                measurement = execute_cpu_task(task);
                break;
            case compute_backend::simd:
                measurement = execute_simd_task(task);
                break;
            case compute_backend::gpu:
                measurement = execute_gpu_task(task);
                break;
            default:
                measurement = execute_cpu_task(task);
                break;
        }
        
        update_performance_history(measurement);
        return measurement;
    }
    
    // Batch execution with load balancing
    template<typename InputType, typename OutputType>
    std::vector<performance_measurement> execute_tasks_parallel(
        std::vector<compute_task<InputType, OutputType>>& tasks) {
        
        std::vector<std::future<performance_measurement>> futures;
        std::vector<performance_measurement> measurements;
        
        futures.reserve(tasks.size());
        measurements.reserve(tasks.size());
        
        // Submit all tasks asynchronously
        for (auto& task : tasks) {
            futures.push_back(std::async(std::launch::async, 
                [this, &task]() -> performance_measurement {
                    return execute_task(task);
                }));
        }
        
        // Collect results
        for (auto& future : futures) {
            measurements.push_back(future.get());
        }
        
        return measurements;
    }
    
    // System information
    [[nodiscard]] const compute_capabilities& get_capabilities() const noexcept {
        return capabilities_;
    }
    
    [[nodiscard]] std::string get_device_info() const {
        std::string info = "Heterogeneous Compute Engine\n";
        info += "CPU Cores: " + std::to_string(capabilities_.cpu_cores) + "\n";
        info += "CPU Performance Score: " + std::to_string(capabilities_.get_cpu_performance_score()) + "\n";
        info += "SIMD Support: " + (capabilities_.has_simd ? "Yes" : "No") + "\n";
        info += "AVX2 Support: " + (capabilities_.has_avx2 ? "Yes" : "No") + "\n";
        info += "AVX512 Support: " + (capabilities_.has_avx512 ? "Yes" : "No") + "\n";
        info += "GPU Available: " + (capabilities_.has_gpu ? "Yes" : "No") + "\n";
        
        if (capabilities_.has_gpu) {
            info += "GPU Compute Units: " + std::to_string(capabilities_.gpu_compute_units) + "\n";
            info += "GPU Memory: " + std::to_string(capabilities_.gpu_memory_mb) + " MB\n";
            info += "GPU Performance Score: " + std::to_string(capabilities_.get_gpu_performance_score()) + "\n";
        }
        
        return info;
    }
    
    [[nodiscard]] std::vector<performance_measurement> get_performance_history() const {
        std::lock_guard<std::mutex> lock(history_mutex_);
        
        std::vector<performance_measurement> combined;
        combined.reserve(cpu_history_.size() + gpu_history_.size());
        
        combined.insert(combined.end(), cpu_history_.begin(), cpu_history_.end());
        combined.insert(combined.end(), gpu_history_.begin(), gpu_history_.end());
        
        return combined;
    }
    
    // Force specific backend for testing/benchmarking
    void set_preferred_backend(compute_backend backend) {
        // Would modify backend selection logic
    }
};

// Type aliases for common configurations
using auto_compute_engine = heterogeneous_compute_engine<compute_backend::auto_detect>;
using cpu_compute_engine = heterogeneous_compute_engine<compute_backend::cpu>;
using gpu_compute_engine = heterogeneous_compute_engine<compute_backend::gpu>;
using simd_compute_engine = heterogeneous_compute_engine<compute_backend::simd>;

// Example spherical operation for concept verification
struct spherical_to_cartesian_operation {
    using input_type = constexpr_spherical_coords<double>;
    using result_type = struct { double x, y, z; double x() const { return x; } double y() const { return y; } double z() const { return z; } };
    
    result_type execute(const input_type& input) const {
        return input.template to_cartesian<result_type>();
    }
};

// Concept verification
static_assert(spherical_operation_concept<spherical_to_cartesian_operation>);

} // namespace compute
} // namespace hsml